"""Collects figures + source tables into deck/assets and deck/deck_data.json for the editable PPTX record."""
import pandas as pd, numpy as np, json, os
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
A="deck/assets"; os.makedirs(A,exist_ok=True)
def img(src,name,split=False,maxw=2400):
    im=Image.open(src).convert("RGB")
    outs=[]
    parts=[im] if not split else [im.crop((0,0,im.width,im.height//2)),im.crop((0,im.height//2,im.width,im.height))]
    for i,p in enumerate(parts):
        if p.width>maxw: p=p.resize((maxw,int(p.height*maxw/p.width)),Image.LANCZOS)
        f=f"{A}/{name}{'_'+str(i+1) if split else ''}.png"; p.save(f,optimize=True); outs.append({"file":os.path.abspath(f),"w":p.width,"h":p.height})
    return outs
d={}
d["img"]={
 "qc_g":img("qc_out/QC_spatial_UMI_Greenwald2024.png","qc_greenwald"),
 "qc_r":img("qc_ravi/QC_spatial_UMI.png","qc_ravi"),
 "prog_g":img("step2a_out/program_maps_all_samples.png","programs_greenwald",split=True),
 "prog_r":img("step2a_ravi/program_maps_all_samples.png","programs_ravi",split=True),
 "cnv":img("cnv_out/compartment_maps.png","cnv_compartments"),
 "niche_maps":img("niche_deconv/niche_maps.png","niche_maps_deconv"),
 "niche_prof":img("niche_deconv/niche_profiles.png","niche_profiles_deconv"),
 "lr_dot":img("lr_out/lr_niche_dotplot.png","lr_dotplot"),
 "old_prog":img("niche_out/niche_profiles.png","superseded_niches_programscore"),
 "old_harm":img("niche_out_h/niche_profiles.png","superseded_niches_harmony"),
}
P=pd.read_csv("niche_deconv/niche_program_profile.csv",index_col=0)
order=sorted(P.index,key=lambda s:int(s.split("_")[1]))
P=P.loc[order]
d["profile"]={"niches":[n.replace("Niche_","N") for n in order],"types":list(P.columns),"z":P.T.round(2).values.tolist()}
F=pd.read_csv("niche_deconv/niche_frequency_by_subtype.csv",index_col=0)[order]
d["freq"]={"niches":[n.replace("Niche_","Niche ") for n in order],"series":{k:(F.loc[k]*100).round(1).tolist() for k in F.index}}
S=pd.read_csv("lr_out/lr_niche_significant.csv"); B=pd.read_csv("lr_out/lr_significant_by_subtype.csv")
S=S.merge(B,on=["niche","pair"],how="left")
S["nnum"]=S.niche.str.split("_").str[1].astype(int)
S=S.sort_values(["nnum","median_auc"],ascending=[True,False])
def rows(df):
    out=[]
    for _,r in df.iterrows():
        f=lambda v:"n/a" if pd.isna(v) else f"{v:.2f}"
        out.append([r.niche.replace("Niche_","N"),r.pair.replace("->"," → ").replace("_","–"),f(r.median_auc),f"{r['frac_auc_gt0.6']*100:.0f}%",
                    f"{int(r.n_sections)} ({int(r.n_greenwald)}/{int(r.n_ravi)})",f(r.median_auc_GBM),f(r["median_auc_IDH-mut"]),
                    str(r.sender_top) if pd.notna(r.sender_top) else "",str(r.receiver_top) if pd.notna(r.receiver_top) else "",f"{r.fdr:.1e}"])
    return out
d["lr6"]=rows(S[S.nnum==6]); d["lr_other"]=rows(S[S.nnum!=6])
d["deconv_valid"]={"types":["Neuron","Oligodendrocyte","All myeloid","MES-like"],"r":[0.77,0.72,0.34,0.50]}
if os.path.exists("ncdeg_out/niche6_malignant_unadjusted_vs_stateadjusted.csv"):
    d["img"]["ncdeg"]=img("ncdeg_out/ncdeg_niche6_malignant.png","ncdeg_niche6_malignant")
    J=pd.read_csv("ncdeg_out/niche6_malignant_unadjusted_vs_stateadjusted.csv"); J=J[J.NCDEG]
    up=J[J.mean_beta>0].sort_values("mean_beta",ascending=False).head(14); dn=J[J.mean_beta<0].sort_values("mean_beta").head(8)
    d["ncdeg_tab"]=[[r.gene,"Up" if r.mean_beta>0 else "Down",f"{r.mean_beta:.2f}",f"{r.fdr:.3f}",f"{r.mean_beta_adj:.2f}",f"{r.p_adj:.4f}",f"{r.sign_consistency_adj*100:.0f}%"] for _,r in pd.concat([up,dn]).iterrows()]
    E=pd.read_csv("ncdeg_out/ncdeg_metaprogram_enrichment_sig.csv")
    E=E[E.meta_program.str.startswith(("Cancer|"))].head(8)
    d["mp_enr"]={"labels":[f"{m.split('|')[1].split(' ',1)[1].strip()} ({dd})" for m,dd in zip(E.meta_program,E.direction)],"values":(-np.log10(E.fdr)).round(1).tolist()}
if os.path.exists("bulk_out/cox_results.csv"):
    d["img"]["bulk_box"]=img("bulk_out/scores_by_type_grade.png","bulk_scores_by_type_grade")
    d["img"]["bulk_km"]=img("bulk_out/KM_HypoxicMES_N6.png","bulk_KM_hypoxicMES")
    C=pd.read_csv("bulk_out/cox_results.csv"); C=C[C.endpoint=="OS"]
    keep=[]
    for g in ["All gliomas","GBM, IDH-wt","Astrocytoma, IDH-mut","Oligodendroglioma","IDH-wt lower-grade"]:
        for v,adj in [("HypoxicMES_N6",None),("ANGPTL4",None),("HypoxicMES_N6","MP_Hypoxia"),("ANGPTL4","MP_Hypoxia")]:
            x=C[(C.group==g)&(C.variable==v)]
            x=x[x.adjusted_for.str.endswith(adj)] if adj else x[~x.adjusted_for.str.contains("MP_|BUFFA")]
            if len(x): r=x.iloc[0]; keep.append([g,v.replace("_N6"," (N6)"),"+ hypoxia program" if adj else "age (+grade, type)",str(int(r.n)),str(int(r.events)),f"{r.HR_per_SD:.2f} ({r.CI_low:.2f}–{r.CI_high:.2f})",f"{r.p:.3f}"])
    d["cox_tab"]=keep
if os.path.exists("targets_out/niche6_target_prioritisation.csv"):
    T=pd.read_csv("targets_out/niche6_target_prioritisation.csv")
    sel=T[T.tier.astype(str).str.startswith("1") | T.gene.isin(["VEGFA","LGALS3","MIF","CD44","ITGA5","TNC","ADM"])]
    f=lambda v,fmt: "" if pd.isna(v) else fmt.format(v)
    d["target_tab"]=[[r.gene,str(r.tier).split(" ")[0],("" if pd.isna(r.niche6_role) else r.niche6_role),f(r.stateadj_beta,"{:.2f}"),f(r.stateadj_p,"{:.3f}"),f(r.OT_GBM_assoc,"{:.2f}"),f(r.DepMap_glioma_median_effect,"{:.2f}"),
        ("" if pd.isna(r.example_drugs) else str(r.example_drugs).split("; ")[0])] for _,r in sel.iterrows()]
if os.path.exists("cgga_out/meta_TCGA_CGGA.csv"):
    d["img"]["cgga_box"]=img("cgga_out/cgga_scores_by_type_grade.png","cgga_scores_by_type_grade")
    d["img"]["cgga_km"]=img("cgga_out/cgga_KM_HypoxicMES_N6.png","cgga_KM")
    d["img"]["meta_forest"]=img("cgga_out/meta_forest.png","meta_forest")
    M=pd.read_csv("cgga_out/meta_TCGA_CGGA.csv"); rows=[]
    for v,adj,lab in [("HypoxicMES_N6","age (+grade, type)","N6 score"),("HypoxicMES_N6","+MP_Hypoxia","N6 beyond hypoxia"),("HypoxicMES_N6","+MP_MES","N6 beyond MES"),("Neuronal_N1","age (+grade, type)","N1 score"),("ANGPTL4","age (+grade, type)","ANGPTL4")]:
        for g in ["All gliomas","Astrocytoma, IDH-mut","GBM, IDH-wt","IDH-wt lower-grade","Oligodendroglioma"]:
            x=M[(M.variable==v)&(M.extra_adjustment==adj)&(M.group==g)]
            if len(x): r=x.iloc[0]; rows.append([lab,g,str(int(r.n_total)),str(int(r.events)),f"{r.HR_meta:.2f} ({r.CI_low:.2f}–{r.CI_high:.2f})",f"{r.p_meta:.1e}" if r.p_meta<0.001 else f"{r.p_meta:.3f}",f"{r.I2:.2f}",r.per_cohort_HR.replace("CGGA","C")])
    d["meta_tab"]=[r for r in rows if r[1] in ("All gliomas","Astrocytoma, IDH-mut","GBM, IDH-wt","Oligodendroglioma") ]
if os.path.exists("prog/progression_niche_score_maps.png"):
    d["img"]["prog"]=img("prog/progression_niche_score_maps.png","progression_case_maps")
if os.path.exists("icb_out/gse121810_arm_scores.png"):
    d["img"]["icb"]=img("icb_out/gse121810_arm_scores.png","icb_gse121810")
if os.path.exists("manuscript_figs/Figure_S5.png"):
    d["img"]["ivy"]=img("manuscript_figs/Figure_S5.png","ivygap_regions")
for k in ["S6","S7","S8","S9"]:
    if os.path.exists(f"manuscript_figs/Figure_{k}.png"): d["img"][k]=img(f"manuscript_figs/Figure_{k}.png",f"fig_{k}")
if os.path.exists("figs_v2/Figure_5.png"):
    d["img"]["fig5v2"]=img("figs_v2/Figure_5.png","fig5_v2")
if os.path.exists("figs_v2/renumbering_map.json"):
    d["renum"]=[[k.replace("_"," ").replace(" new"," (new ANGPTL4 axis)"),v.replace("_"," ")] for k,v in json.load(open("figs_v2/renumbering_map.json")).items()]
json.dump(d,open("deck/deck_data.json","w"),indent=1)
print(len(d["lr6"]),len(d["lr_other"]))
