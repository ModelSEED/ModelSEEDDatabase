#!/usr/bin/env python3
"""Compare dGPredictor-ModelSEED (retrained) reaction energies to TECRDB.

Matching is done on a "SMILES key" derived from each compound's SMILES via
RDKit, at TWO tiers:
  * stereo_exact : full InChIKey of the charge-NEUTRALISED, largest-fragment
                   molecule (distinguishes anomers/stereoisomers; charge/pH
                   independent -- ATP == ATP4-).
  * skeleton     : InChIKey connectivity block (first 14 chars); unifies
                   stereoisomers & protonation states (e.g. all aldohexopyranoses).
Both sides (ModelSEED cpd SMILES and TECRDB KEGG SMILES) go through the SAME
RDKit pipeline, so the keys are directly comparable.

Reaction key = (reactant multiset, product multiset) of (key, coeff), protons
dropped, matched forward or reverse (sign flipped for reverse).
TECRDB tested energy: dG'o = -R*T*ln(K'), aggregated (median) per reaction.
dGPredictor energy: staged modelseed_retrained_dG.json dG_mean (kJ/mol),
for the reaction as written in ModelSEED. Disparity = dGpred - TECRDB.
"""
import json, glob, re, math, collections, csv, sys, os
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
# PORTED 2026-09-07 from Cooper Taylor's build_comparison.py
# (/scratch/ctaylor/dgpredictor_tecrdb/scripts/). The matching logic below is
# HIS and is unchanged; only the paths are rewritten to derive from this
# repository so the file it produces is regenerable here.
#
# Why: tecrdb_vs_dgpredictor_modelseed.csv is the ground truth behind every
# accuracy number in the paper -- 802 stereo-exact anchors -- and it was read
# by two committed scripts and produced by none. A release whose headline
# figure cannot be regenerated from the repository is not reproducible, and
# data_availability.tex claims otherwise.
#
# TWO DELIBERATE DIFFERENCES FROM HIS RUN, both of which move the numbers:
#   * the ModelSEED snapshot is THIS repository's live Biochemistry/, not his
#     pinned June copy;
#   * TECRDB comes from our tecrdb_v11 (1.66 MB) because his data/TECRDB.csv is
#     mode 600 and unreadable. Same four columns are used -- reaction,
#     temperature, p_h, K_prime -- so the schema is compatible, but it is a
#     different vintage.
# Expect the match counts to differ from his 1,550 / 802. Compare before
# adopting the output as the anchor set.
# ---------------------------------------------------------------------------
REPO   = Path(os.environ.get("MSDB_ROOT", Path(__file__).resolve().parents[3]))
SNAP   = Path(os.environ.get("MSDB_BIOCHEM", REPO / "Biochemistry"))
MSDB   = SNAP.parent
STAGED = Path(os.environ.get(
    "DGPREDICTOR_JSON", SNAP / "Thermodynamics" / "dGPredictor" / "retrained_dG.json"))
TECRDB = Path(os.environ.get(
    "TECRDB_SOURCE",
    "/scratch/seaver/Claude_Projects/eQuilibrator/data/tecrdb_v11/TECRDB.csv"))
OUTDIR = Path(os.environ.get("TECRDB_COMPARISON_OUT", SNAP / "Thermodynamics" / "SourceGrading"))
R_KJ   = 8.314462618e-3
PROTON = "GPRLSGONYQIRFK"      # InChIKey14 of [H+]
_UN    = rdMolStandardize.Uncharger()
_cache = {}

def keys_from_smiles(smiles):
    """-> (full_neutral_inchikey, skeleton_block1) or (None,None)."""
    if not smiles or smiles in ("null",""): return (None,None)
    if smiles in _cache: return _cache[smiles]
    res=(None,None)
    m=Chem.MolFromSmiles(smiles)
    if m is not None:
        try:
            skel=Chem.MolToInchiKey(m).split("-")[0]
            try:
                mp=rdMolStandardize.FragmentParent(m); mp=_UN.uncharge(mp)
                full=Chem.MolToInchiKey(mp)
            except Exception:
                full=Chem.MolToInchiKey(m)
            res=(full,skel)
        except Exception:
            res=(None,None)
    _cache[smiles]=res
    return res

_canon={}
def canonical_smiles(smiles):
    """RDKit-canonical SMILES of the served structure (as stored). None if unparseable."""
    if not smiles or smiles in ("null",""): return None
    if smiles in _canon: return _canon[smiles]
    m=Chem.MolFromSmiles(smiles)
    out=Chem.MolToSmiles(m) if m is not None else None
    _canon[smiles]=out
    return out

def reaction_smiles(stoich, cpd_smiles):
    """Build (reactants_str, products_str, reaction_smiles, complete) in the
    ModelSEED-written direction. reactants/products_str are coeff-annotated and
    ';'-joined; reaction_smiles is a RDKit-style 'r.r.p>>...' with each species
    repeated by its (integer) coefficient. complete=False if any compound lacks a
    parseable SMILES (reaction_smiles then left blank)."""
    react=[]; prod=[]; r_expand=[]; p_expand=[]; complete=True
    for s in stoich:
        c=s["coefficient"]; smi=canonical_smiles(cpd_smiles.get(s["compound"]))
        if smi is None:
            complete=False; smi_disp=f"?{s['compound']}"
        else:
            smi_disp=smi
        n=abs(c)
        side_ann = react if c<0 else prod
        side_ann.append(f"{n:g} {smi_disp}")
        if smi is not None:
            reps = int(round(n)) if abs(n-round(n))<1e-6 else 1
            (r_expand if c<0 else p_expand).extend([smi]*max(reps,1))
    rs = ".".join(r_expand)+">>"+".".join(p_expand) if complete else ""
    return "; ".join(react), "; ".join(prod), rs, complete

# ---- compounds ----
print("loading compounds ...",file=sys.stderr)
cpd_name={};cpd_smiles={};cpd_inchikey={};cpd_formula={};cpd_charge={}
for f in glob.glob(str(SNAP / "compound_*.json")):
    for c in json.load(open(f)):
        cid=c["id"]
        cpd_name[cid]=c.get("name");cpd_smiles[cid]=c.get("smiles")
        cpd_inchikey[cid]=c.get("inchikey");cpd_formula[cid]=c.get("formula")
        cpd_charge[cid]=c.get("charge")

# ---- reactions ----
print("loading reactions ...",file=sys.stderr)
rxn=collections.OrderedDict()
for f in sorted(glob.glob(str(SNAP / "reaction_*.json"))):
    for r in json.load(open(f)): rxn[r["id"]]=r
staged=json.load(open(STAGED))

# ---- cpd -> keys (only cpds in predicted reactions) ----
needed=set()
for rid in staged:
    r=rxn.get(rid)
    if r:
        for s in r.get("stoichiometry",[]): needed.add(s["compound"])
cpd_full={};cpd_skel={}
for cid in needed:
    full,skel=keys_from_smiles(cpd_smiles.get(cid))
    if skel is None and cpd_inchikey.get(cid):
        skel=cpd_inchikey[cid].split("-")[0]; full=full or cpd_inchikey[cid]
    if skel: cpd_skel[cid]=skel
    if full: cpd_full[cid]=full
skel2cpds=collections.defaultdict(set)
for cid,k in cpd_skel.items(): skel2cpds[k].add(cid)
print(f"  cpds needed={len(needed)}  with skeleton={len(cpd_skel)}  with full={len(cpd_full)}",file=sys.stderr)

# ---- KEGG -> keys ----
kegg_smi_charged=collections.defaultdict(set);kegg_smi_orig=collections.defaultdict(set)
with open(MSDB / "Biochemistry" / "Structures" / "All_ModelSEED_Structures.txt") as fh:
    for line in fh:
        p=line.rstrip("\n").split("\t")
        if len(p)<8: continue
        cpd,typ,charged,alias,source,formula,charge,struct=p[:8]
        if source!="KEGG" or not alias.startswith("C") or struct in ("","null"): continue
        if typ=="SMILE":
            (kegg_smi_charged if charged=="Charged" else kegg_smi_orig)[alias].add(struct)
kegg2cpd=collections.defaultdict(set)
with open(MSDB / "Biochemistry" / "Aliases" / "Unique_ModelSEED_Compound_Aliases.txt") as fh:
    next(fh,None)
    for line in fh:
        p=line.rstrip("\n").split("\t")
        if len(p)>=3 and p[2]=="KEGG" and p[1].startswith("C"): kegg2cpd[p[1]].add(p[0])

def kegg_keys(cid):
    """(full,skel) preferring a microspecies matching a ModelSEED reaction cpd."""
    cands=[]
    for smi in list(kegg_smi_charged.get(cid,()))+list(kegg_smi_orig.get(cid,())):
        f,s=keys_from_smiles(smi)
        if s: cands.append((f,s))
    for f,s in cands:
        if s in skel2cpds: return (f,s)          # skeleton matches a reaction cpd
    if cands: return cands[0]
    for c in kegg2cpd.get(cid,()):               # fallback: alias cpd
        if c in cpd_skel: return (cpd_full.get(c),cpd_skel[c])
    return (None,None)

# ---- parse TECRDB ----
df=pd.read_csv(TECRDB)
def parse_side(side, tier):
    d=collections.defaultdict(float); unresolved=[]
    for part in side.split("+"):
        part=part.strip()
        if not part: continue
        m=re.match(r'^(\d+(?:\.\d+)?)\s+(.*)$',part)
        coeff=float(m.group(1)) if m else 1.0
        tok=(m.group(2) if m else part).strip()
        if not tok.startswith("kegg:"): unresolved.append(tok); continue
        full,skel=kegg_keys(tok.split("kegg:")[1])
        k = full if tier=="full" else skel
        if k is None: unresolved.append(tok); continue
        if skel==PROTON: continue
        d[k]+=coeff
    return d,unresolved
def norm(d): return tuple(sorted((k,round(v,4)) for k,v in d.items()))

def ms_side(stoich, want_react, tier):
    d=collections.defaultdict(float)
    for s in stoich:
        c=s["coefficient"]
        if (c<0)!=want_react: continue
        skel=cpd_skel.get(s["compound"]); full=cpd_full.get(s["compound"])
        k = full if tier=="full" else skel
        if k is None: return None
        if skel==PROTON: continue
        d[k]+=abs(c)
    return d

def build_tecr(tier):
    groups={}
    for _,row in df.iterrows():
        rs=row["reaction"]
        if not isinstance(rs,str) or "=" not in rs: continue
        left,right=rs.split("=",1)
        L,uL=parse_side(left,tier); Rr,uR=parse_side(right,tier)
        if uL or uR or not L or not Rr: continue
        Ln,Rn=norm(L),norm(Rr)
        canon=(Ln,Rn,True) if Ln<=Rn else (Rn,Ln,False)
        key=(canon[0],canon[1]); fwd=canon[2]
        g=groups.setdefault(key,dict(react=canon[0],prod=canon[1],dGs=[],phs=[],Ts=[],
              ECs=set(),enz=set(),kegg_rxns=set()))
        g["kegg_rxns"].add(rs);g["ECs"].add(str(row.get("EC")));g["enz"].add(str(row.get("enzyme_name")))
        Kp=row.get("K_prime")
        if pd.notna(Kp) and Kp>0:
            T=float(row["temperature"]) if pd.notna(row.get("temperature")) else 298.15
            dg=-R_KJ*T*math.log(float(Kp)); g["dGs"].append(dg if fwd else -dg)
            if pd.notna(row.get("p_h")): g["phs"].append(float(row["p_h"]))
            g["Ts"].append(T)
    index={}
    for key,g in groups.items():
        index[(g["react"],g["prod"])]=(key,True)
        index[(g["prod"],g["react"])]=(key,False)
    return groups,index

def med(xs):
    xs=sorted(xs); n=len(xs)
    return None if n==0 else (xs[n//2] if n%2 else (xs[n//2-1]+xs[n//2])/2)

tecr={};idx={}
for tier in ("full","skel"):
    tecr[tier],idx[tier]=build_tecr(tier)
    print(f"  TECRDB groups[{tier}]={len(tecr[tier])}",file=sys.stderr)

# ---- match MS reactions ----
matches={}   # rid -> row (prefer stereo_exact)
for rid in staged:
    r=rxn.get(rid)
    if not r: continue
    for tier,label in (("full","stereo_exact"),("skel","skeleton")):
        L=ms_side(r["stoichiometry"],True,tier); Rr=ms_side(r["stoichiometry"],False,tier)
        if L is None or Rr is None or not L or not Rr: continue
        hit=idx[tier].get((norm(L),norm(Rr)))
        if not hit: continue
        key,ms_same=hit; g=tecr[tier][key]
        if not g["dGs"]: continue
        tecr_canon=med(g["dGs"]); tecr_ms=tecr_canon if ms_same else -tecr_canon
        dgs_ms=g["dGs"] if ms_same else [-x for x in g["dGs"]]
        dgpred=staged[rid]["dG_mean"];dgerr=staged[rid].get("dG_uncer")
        diff=dgpred-tecr_ms
        # combined uncertainty = prediction self-error (+) experimental scatter (quadrature)
        _mean=sum(dgs_ms)/len(dgs_ms)
        exp_sd=(sum((x-_mean)**2 for x in dgs_ms)/(len(dgs_ms)-1))**0.5 if len(dgs_ms)>=2 else 0.0
        pred_err=dgerr if dgerr is not None else 0.0
        combined_err=math.sqrt(pred_err**2+exp_sd**2)
        significant = abs(diff) > combined_err   # disparity exceeds the combined error
        # other thermo methods (stored kcal, ModelSEED written direction) -> kJ
        th=r.get("thermodynamics") or {}
        def _kj(name):
            v=th.get(name)
            return round(v[0]*4.184,3) if isinstance(v,list) and v and isinstance(v[0],(int,float)) else None
        gc_kj=_kj("Group contribution"); dgp_orig_kj=_kj("dGPredictor"); eq_kj=_kj("eQuilibrator")
        rsmi_r, rsmi_p, rsmi, rsmi_ok = reaction_smiles(r["stoichiometry"], cpd_smiles)
        row=dict(modelseed_rxn=rid,name=r.get("name"),
            reactants_smiles=rsmi_r, products_smiles=rsmi_p,
            reaction_smiles=rsmi, reaction_smiles_complete=rsmi_ok,
            ec=";".join(sorted(x for x in g["ECs"] if x and x!="nan")),
            enzyme_name=";".join(sorted(x for x in g["enz"] if x and x!="nan"))[:200],
            equation_definition=r.get("definition"),equation_ids=r.get("equation"),
            tecrdb_reaction=" | ".join(sorted(g["kegg_rxns"]))[:300],
            n_measurements=len(g["dGs"]),
            dGpredictor_modelseed_dG_kJ=round(dgpred,3),
            tecrdb_dG_kJ=round(tecr_ms,3),
            diff_kJ=round(diff,3),abs_diff_kJ=round(abs(diff),3),
            dGpredictor_modelseed_dG_kcal=round(dgpred/4.184,3),
            tecrdb_dG_kcal=round(tecr_ms/4.184,3),
            diff_kcal=round(diff/4.184,3),abs_diff_kcal=round(abs(diff)/4.184,3),
            dGpredictor_modelseed_err_kJ=round(dgerr,3) if dgerr is not None else None,
            tecrdb_dG_sd_kJ=round(exp_sd,3),combined_err_kJ=round(combined_err,3),significant=significant,
            tecrdb_dG_min_kJ=round(min(dgs_ms),3),tecrdb_dG_max_kJ=round(max(dgs_ms),3),
            other_GroupContribution_dG_kJ=gc_kj,
            other_dGPredictor_original_dG_kJ=dgp_orig_kj,
            other_eQuilibrator_dG_kJ=eq_kj,
            pH_min=round(min(g["phs"]),2) if g["phs"] else None,
            pH_max=round(max(g["phs"]),2) if g["phs"] else None,
            T_min=round(min(g["Ts"]),2) if g["Ts"] else None,
            T_max=round(max(g["Ts"]),2) if g["Ts"] else None,
            ms_orientation_vs_canonical=("same" if ms_same else "reversed"),
            match_tier=label,
            structure_sig=norm(ms_side(r["stoichiometry"],True,"full") or {})+norm(ms_side(r["stoichiometry"],False,"full") or {}))
        # keep the stereo_exact match if available; else skeleton
        if rid not in matches or label=="stereo_exact":
            matches[rid]=row
        break_needed = (label=="stereo_exact")
        if break_needed: break
rows=list(matches.values())
rows.sort(key=lambda x:-x["abs_diff_kJ"])
print(f"matched MS reactions: {len(rows)}  "
      f"(stereo_exact={sum(1 for r in rows if r['match_tier']=='stereo_exact')}, "
      f"skeleton={sum(1 for r in rows if r['match_tier']=='skeleton')})",file=sys.stderr)

cols=["modelseed_rxn","name","ec","enzyme_name","equation_definition","equation_ids",
      "tecrdb_reaction","n_measurements","match_tier",
      "dGpredictor_modelseed_dG_kJ","tecrdb_dG_kJ","diff_kJ","abs_diff_kJ",
      "dGpredictor_modelseed_dG_kcal","tecrdb_dG_kcal","diff_kcal","abs_diff_kcal",
      "dGpredictor_modelseed_err_kJ","tecrdb_dG_sd_kJ","combined_err_kJ","significant",
      "other_GroupContribution_dG_kJ","other_dGPredictor_original_dG_kJ","other_eQuilibrator_dG_kJ",
      "tecrdb_dG_min_kJ","tecrdb_dG_max_kJ",
      "pH_min","pH_max","T_min","T_max","ms_orientation_vs_canonical",
      "reactants_smiles","products_smiles","reaction_smiles","reaction_smiles_complete"]
OUTDIR.mkdir(parents=True, exist_ok=True)
out=str(OUTDIR / "tecrdb_vs_dgpredictor_modelseed.csv")
with open(out,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=cols); w.writeheader()
    for m in rows: w.writerow({k:m.get(k) for k in cols})
print("WROTE",out)

# ---- dedup to distinct chemistries, save structure_sig map for top-10 step ----
json.dump({m["modelseed_rxn"]:list(m["structure_sig"]) for m in rows},
          open(OUTDIR / "_structure_sig.json","w"))
diag=dict(tecrdb_rows_total=int(len(df)),
    predicted_reactions=len(staged),matched_reactions=len(rows),
    stereo_exact=sum(1 for r in rows if r["match_tier"]=="stereo_exact"),
    skeleton_only=sum(1 for r in rows if r["match_tier"]=="skeleton"),
    tecr_groups_full=len(tecr["full"]),tecr_groups_skel=len(tecr["skel"]))
json.dump(diag,open(OUTDIR / "tecrdb_comparison_diagnostics.json","w"),indent=2)
print(json.dumps(diag,indent=2))
