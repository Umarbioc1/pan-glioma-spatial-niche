"""25_supp_methods.py -- Supplementary Methods + supplementary figure legends (Neuro-Oncology numbering)."""
import re, json
v2=open('ms_v2.md').read()
mp={int(k):v for k,v in json.load(open('refmap_no.json')).items()}
US=[('tumours','tumors'),('tumour','tumor'),('Tumour','Tumor'),('signalling','signaling'),('Signalling','Signaling'),('neighbourhood','neighborhood'),('neighbours','neighbors'),
('normalised','normalized'),('Normalised','Normalized'),('normalisation','normalization'),('initialisations','initializations'),('summarised','summarized'),('analysed','analyzed'),
('prioritisation','prioritization'),('Prioritisation','Prioritization'),('colour','color'),('Colour','Color'),('labelled','labeled'),('modelled','modeled'),('centre','center'),('grey','gray'),('Standardised','Standardized'),('standardised','standardized'),('characterisation','characterization')]
def us(s):
    for a,b in US: s=s.replace(a,b)
    return s
def refs(s):
    def rep(m):
        nums=[]
        for p in re.split(r',\s*',m.group(1)):
            if '–' in p: a,b=map(int,p.split('–')); nums+=range(a,b+1)
            else: nums.append(int(p))
        return '^'+','.join(str(mp[n]) for n in sorted(nums,key=lambda n:mp[n]))+'^'
    return re.sub(r'\[(\d+(?:[,–]\s?\d+)*)\]',rep,s)
def pv(s):
    s=re.sub(r'\bp ([=<>]) 0\.(\d+)',r'P \1 .\2',s); s=re.sub(r'\*p < 0\.05','*P < .05',s); s=re.sub(r'\*\*p < 0\.01','**P < .01',s); s=re.sub(r'\*\*\*p < 0\.001','***P < .001',s)
    return s
meth=v2.split('## Methods')[1].split('## Declarations')[0].strip()
meth=pv(refs(us(meth)))
leg=v2.split('### Supplementary figures')[1].split('### Supplementary tables')[0].strip()
leg=pv(refs(us(leg)))
f7=re.search(r'\*\*Figure 7\. ([^*]+)\*\*\s*(.+)',v2.split('## Figure legends')[1])
s10=f"**Figure S10. {us(f7.group(1))}** {pv(us(f7.group(2)))}"
tabs=v2.split('### Supplementary tables')[1].split('## Before submission')[0].strip()
tabs=us(tabs).replace('(Supplementary\\_Tables\\_S1-S16.xlsx)','')
md=f"""---
title: "Supplementary Material"
subtitle: "Pan-glioma spatial transcriptomics separates niche from cell state and hypoxia and nominates ANGPTL4 in the glioblastoma peri-necrotic niche"
---

# Supplementary Methods

Reference numbers refer to the reference list of the main article.

{meth.replace('### ','## ')}

# Supplementary Figure Legends

Supplementary Figures S1–S10 are provided as separate files.

{leg}

{s10}

# Supplementary Tables

Supplementary Tables S1–S16 are provided in Supplementary_Tables_S1-S16.xlsx.

{tabs}
"""
open('no_sub/Supplementary_Material.md','w').write(md)
print(len(md)); print(re.findall(r'\[\d',md)[:5])
