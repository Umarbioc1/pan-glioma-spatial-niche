"""24_neurooncology_build.py -- renumber citations by first appearance, count words, assemble Neuro-Oncology manuscript markdown."""
import re
src=open('ms_no_src.md').read()
refs_old=open('refs_ama.txt').read().strip().split('\n')
def expand(s):
    out=[]
    for p in s.split(','):
        if '-' in p: a,b=map(int,p.split('-')); out+=range(a,b+1)
        else: out.append(int(p))
    return out
order=[]
for m in re.finditer(r'\[\[([\d,\-]+)\]\]',src):
    for n in expand(m.group(1)):
        if n not in order: order.append(n)
print('cited',len(order),'missing',set(range(1,38))-set(order))
mp={o:i+1 for i,o in enumerate(order)}
def fmt(nums):
    nums=sorted(mp[n] for n in nums); out=[]; i=0
    while i<len(nums):
        j=i
        while j+1<len(nums) and nums[j+1]==nums[j]+1: j+=1
        out.append(f"{nums[i]}–{nums[j]}" if j-i>=2 else ",".join(map(str,nums[i:j+1])))
        i=j+1
    return ",".join(out)
body=re.sub(r'\[\[([\d,\-]+)\]\]',lambda m:'^'+fmt(expand(m.group(1)))+'^',src)
newrefs=[]
for o in order:
    r=refs_old[o-1].split('. ',1)[1].replace('–','-')
    newrefs.append(f"{mp[o]}. {r}")
def wc(s):
    s=re.sub(r'\^[\d,–]+\^','',s); s=re.sub(r'^#.*$','',s,flags=re.M)
    return len(re.findall(r"\S+",s))
abs_=body.split('## Abstract')[1].split('## Key Points')[0]
txt=body.split('## Introduction')[1].split('## Discussion')[1]  # placeholder
txt=body.split('## Introduction')[0]  # unused
main='## Introduction'+body.split('## Introduction')[1]
A=wc(abs_.replace('**Background.**','').replace('**Methods.**','').replace('**Results.**','').replace('**Conclusions.**',''))
T=wc(main)
kp=re.findall(r'^- (.+)$',body.split('## Key Points')[1].split('## Importance')[0],re.M)
imp=wc(body.split('## Importance of the Study')[1].split('## Introduction')[0])
print('abstract',A,'text',T,'importance',imp,'keypoints',[len(k) for k in kp],sum(len(k) for k in kp))
for s in ['Introduction','Materials and Methods','Results','Discussion']:
    print(s, wc(body.split('## '+s)[1].split('\n## ')[0]))
title=re.search(r'\| Title \| (.+) \|',body).group(1); rt=re.search(r'\| Running title \| (.+) \|',body).group(1)
print('title chars',len(title),'running',len(rt))
body=body.replace('{ABS}',str(A)).replace('{TXT}',f"{T:,}")
open('ms_no_body.md','w').write(body); open('refs_no.txt','w').write('\n'.join(newrefs))
import json; json.dump(mp,open('refmap_no.json','w'))
