#!/usr/bin/env python
import argparse, json, os, warnings
import scanpy as sc
import pandas as pd
from scib.metrics import ilisi_graph, kbett, silhouette_batch, silhouette_label, graph_connectivity

warnings.filterwarnings("ignore")

def run_harmony(adata, key_batch):
    import harmonypy as hm
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, key_batch)
    adata.obsm["X_harmony"] = ho.Z_corr.T
    sc.pp.neighbors(adata, use_rep="X_harmony")
    sc.tl.umap(adata)
    return adata

def run_scanorama(adata, key_batch):
    import scanorama, numpy as np
    adatas = [adata[adata.obs[key_batch]==b].copy() for b in adata.obs[key_batch].unique()]
    corrected, _ = scanorama.correct_scanpy(adatas, return_dimred=True)
    order = pd.Index([])
    for a in corrected: order = order.append(pd.Index(a.obs_names))
    X = pd.concat([pd.DataFrame(a.obsm["X_scanorama"], index=a.obs_names) for a in corrected]).loc[adata.obs_names]
    adata.obsm["X_scanorama"] = X.values
    sc.pp.neighbors(adata, use_rep="X_scanorama")
    sc.tl.umap(adata)
    return adata

def compute_scib(adata, batch_key, label_key, rep):
    sc.pp.neighbors(adata, use_rep=rep)
    scores = {
        "iLISI": ilisi_graph(adata, batch_key=batch_key),
        "kBET": kbett(adata, batch_key=batch_key, label_key=label_key),
        "Silhouette_batch": silhouette_batch(adata, batch_key=batch_key),
        "Silhouette_label": silhouette_label(adata, label_key=label_key),
        "Graph_connectivity": graph_connectivity(adata, label_key=label_key),
    }
    return scores

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="JSON mapping name->h5ad path")
    ap.add_argument("--outdir", default="results/integration_benchmark")
    ap.add_argument("--batch_key", default="batch")
    ap.add_argument("--label_key", default="celltype")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    datasets = json.load(open(args.input))
    all_scores = []
    for name, path in datasets.items():
        adata = sc.read_h5ad(path)
        # Harmony
        ah = adata.copy(); ah = run_harmony(ah, args.batch_key)
        hs = compute_scib(ah, args.batch_key, args.label_key, "X_harmony"); hs["method"]="Harmony"; hs["dataset"]=name; all_scores.append(hs)
        ah.write(os.path.join(args.outdir, f"{name}_harmony.h5ad"))
        # Scanorama
        asn = adata.copy(); asn = run_scanorama(asn, args.batch_key)
        ss = compute_scib(asn, args.batch_key, args.label_key, "X_scanorama"); ss["method"]="Scanorama"; ss["dataset"]=name; all_scores.append(ss)
        asn.write(os.path.join(args.outdir, f"{name}_scanorama.h5ad"))
    df = pd.DataFrame(all_scores)
    df.to_csv(os.path.join(args.outdir, "scib_scores.csv"), index=False)
    print(df.groupby(["dataset","method"]).mean(numeric_only=True))

if __name__ == "__main__":
    main()
