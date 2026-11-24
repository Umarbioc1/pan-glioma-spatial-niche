"""24a_convert_style.py -- convert v2 Results/Discussion to Neuro-Oncology style (US spelling, AMA P values, figure labels, Fig 7 -> Fig S10)."""
import re
v2=open('ms_v2.md').read()
src=open('ms_no_src.md').read()
US=[('tumours','tumors'),('tumour','tumor'),('Tumours','Tumors'),('Tumour','Tumor'),('signalling','signaling'),('Signalling','Signaling'),
('neighbourhood','neighborhood'),('organise','organize'),('organis','organiz'),('oedema','edema'),('colour','color'),('analysed','analyzed'),('analyse','analyze'),
('prioritisation','prioritization'),('Prioritisation','Prioritization'),('prioritis','prioritiz'),('normalis','normaliz'),('localis','localiz'),('antagonise','antagonize'),
('labelled','labeled'),('modelled','modeled'),('centre','center'),('towards','toward'),('recognis','recogniz'),('characteris','characteriz'),('haem','hem'),('favour','favor'),('behaviour','behavior')]
def us(s):
    for a,b in US: s=s.replace(a,b)
    return s
def pval(s):
    s=re.sub(r'\bp ([=<>]) 0\.(\d+)',r'P \1 .\2',s)
    s=re.sub(r'\bp ([=<>]) (\d)',r'P \1 \2',s)
    s=re.sub(r'\bp-values?\b','P values',s)
    s=re.sub(r'(\d)×10',r'\1 × 10',s)
    s=re.sub(r'all p > ','all P > ',s)
    return s
def figs(s):
    s=re.sub(r'\(Fig 7; Table S16\)','(Supplementary Figure S10; Supplementary Table S16)',s)
    s=re.sub(r'\bFig 7\b','Supplementary Figure S10',s)
    s=re.sub(r'\bFig (S\d+[A-C]?)(, (S\d+))?',lambda m:'Supplementary Figure'+('s ' if m.group(2) else ' ')+m.group(1)+(' and '+m.group(3) if m.group(2) else ''),s)
    s=re.sub(r'\bFig (\d)',r'Figure \1',s)
    s=re.sub(r'\bTables (S\d+), (S\d+)',r'Supplementary Tables \1 and \2',s)
    s=re.sub(r'\bTables (S\d+–S\d+)',r'Supplementary Tables \1',s)
    s=re.sub(r'(?<!Supplementary )\bTable (S\d+)',r'Supplementary Table \1',s)
    return s
def refs(s):
    return re.sub(r'\[(\d+(?:[,–]\s?\d+)*)\]',lambda m:'[['+m.group(1).replace('–','-').replace(' ','')+']]',s)
def clean(sec):
    sec='\n'.join(l for l in sec.split('\n') if not l.startswith('![') and not re.match(r'\*\*Figure \d\.',l))
    sec=re.sub(r'\n{3,}','\n\n',sec)
    return refs(figs(pval(us(sec)))).strip()
R=clean(v2.split('## Results')[1].split('## Discussion')[0])
D=clean(v2.split('## Discussion')[1].split('## Methods')[0])
# drop the 'Four points deserve comment' enumerations' style? keep.
pre=src.split('## Results')[0]
out=pre+'## Results\n\n'+R+'\n\n## Discussion\n\n'+D+'\n'
open('ms_no_src.md','w').write(out)
print(re.findall(r'Supplementary Figures? S\d+[^;)]*',R)[:30])
