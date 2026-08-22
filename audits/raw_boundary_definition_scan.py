#!/usr/bin/env python3
"""Compare plausible raw-SEG boundary definitions against historical window stability."""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.feature_selection import chi2

SOURCE = {3:["bp_1p_bin19","bp_1p_bin20","bp_6q_bin19","bp_6q_bin20","bp_4p_bin20"],
          5:["bp_1p_bin20","bp_6q_bin19","bp_19q_bin11","bp_19q_bin10","bp_2p_bin16"],
          7:["bp_1p_bin19","bp_1p_bin20","bp_6q_bin20","bp_19q_bin11","bp_3p_bin11"]}

def endpoint(d,h):
    t=pd.to_numeric(d.rfs_time_years,errors="coerce").to_numpy();e=pd.to_numeric(d.rfs_event,errors="coerce").to_numpy()
    return np.where((e==1)&(t<=h),1,np.where(t>=h,0,np.nan))

def make_boundaries(seg, arms, which, exclude_centromere=False):
    parts=[]
    choices=("Start","End") if which=="both" else (which.capitalize(),)
    for col in choices:
        z=seg[["cg","Chr",col,"Event"]].rename(columns={col:"pos"}).copy();z["side"]=col.lower();parts.append(z)
    b=pd.concat(parts,ignore_index=True);b["Chr"]=b.Chr.astype(str).str.replace("^chr","",regex=True)
    b=b.merge(arms,on="Chr",how="inner");b=b[(b.pos>=b.start)&(b.pos<=b.end)].copy()
    # Current feature builder excludes only chromosome-arm telomeres.
    b=b[~((b.Arm=="p")&(b.pos==b.start)) & ~((b.Arm=="q")&(b.pos==b.end))]
    if exclude_centromere:
        b=b[~((b.Arm=="p")&(b.pos==b.end)) & ~((b.Arm=="q")&(b.pos==b.start))]
    b["rel"]=np.where(b.Arm=="p",(b.pos-b.start)/b.size_arm,(b.end-b.pos)/b.size_arm).clip(0,.999999)
    b["bin"]=(np.floor(b.rel*20)+1).astype(int)
    b["feature"]="bp_"+b.Chr+b.Arm+"_bin"+b["bin"].astype(str).str.zfill(2)
    return b.drop_duplicates(["cg","feature"])

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segments", required=True, type=Path)
    parser.add_argument("--arm-reference", required=True, type=Path)
    parser.add_argument("--model-data", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(args.model_data,low_memory=False);ids=set(d.cg.astype(str))
    seg=pd.read_csv(args.segments,sep="\t",usecols=["cg","Chr","Start","End","Event"],low_memory=False)
    seg=seg[seg.cg.astype(str).isin(ids)&(seg.Event!="LOH")].copy()
    arms=pd.read_csv(args.arm_reference);arms["Chr"]=arms.Chr.astype(str).str.replace("^chr","",regex=True)
    for c in ["size_arm","start","end"]:arms[c]=pd.to_numeric(arms[c].astype(str).str.replace(",",""),errors="coerce")
    allowed=set([f"{n}{a}" for n in range(1,23) for a in "pq"])-{"13p","14p","15p","21p","22p"}
    arms=arms[(arms.Chr+arms.Arm).isin(allowed)]
    rows=[]
    for variant,which,exc in [("current_both","both",False),("both_internal_only","both",True),
                              ("starts_only","start",False),("ends_only","end",False),
                              ("starts_internal_only","start",True),("ends_internal_only","end",True)]:
        b=make_boundaries(seg,arms,which,exc)
        wide=pd.crosstab(b.cg,b.feature).reindex(index=d.cg.astype(str),fill_value=0)
        x=(wide.to_numpy()>0).astype(np.uint8);names=np.asarray(wide.columns)
        support=x.sum(0);keep=support>=20;x=x[:,keep];names=names[keep]
        for h in (3,5,7):
            y=endpoint(d,h);ok=np.isfinite(y);s=chi2(x[ok],y[ok].astype(int))[0]
            order=np.argsort(-np.nan_to_num(s,nan=-np.inf),kind="stable")
            rank=np.empty(len(order),int);rank[order]=np.arange(1,len(order)+1)
            for f in SOURCE[h]:
                hit=np.where(names==f)[0]
                rows.append({"variant":variant,"horizon":h,"n_candidates":len(names),"feature":f,
                             "rank":int(rank[hit[0]]) if len(hit) else np.nan,
                             "support":int(support[keep][hit[0]]) if len(hit) else np.nan,
                             "chi2":float(s[hit[0]]) if len(hit) else np.nan})
    out=pd.DataFrame(rows);out.to_csv(args.output_dir/"source_window_ranks.csv",index=False)
    print(out.to_string(index=False))

if __name__=="__main__":main()
