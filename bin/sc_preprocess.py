#!/usr/bin/env python3
"""Single-cell preprocessing for AnnData (.h5ad) and MuData (.h5mu).

In multimodal MuData, GEX and protein are preprocessed separately, and
weighted nearest neighbors (WNN) integrates modality-specific PCA embeddings.
Only explicitly requested output files are written.
"""
from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import sparse

LOG = logging.getLogger("preprocess")
GEX_NAMES = {"gex", "rna", "gene_expression", "gene-expression"}
PRO_NAMES = {"pro", "prot", "protein", "adt", "citeseq", "cite-seq"}


def read_dataset(path):
    import scanpy as sc
    import mudata as md

    if path.suffix.lower() == ".h5ad":
        return sc.read_h5ad(path)
    if path.suffix.lower() == ".h5mu":
        return md.read_h5mu(path)
    raise ValueError("Input must have a .h5ad or .h5mu extension")


def is_mudata(data):
    from mudata import MuData
    return isinstance(data, MuData)


def modality_kind(name, adata, standalone=False):
    lower = name.lower()
    if lower in GEX_NAMES:
        return "gex"
    if lower in PRO_NAMES:
        return "protein"
    if "feature_types" in adata.var:
        types = set(adata.var["feature_types"].dropna().astype(str).unique())
        if types == {"Antibody Capture"}:
            return "protein"
        if types == {"Gene Expression"}:
            return "gex"
    if standalone:
        LOG.warning("Unrecognized AnnData modality: treating it as GEX (use --data-type protein for ADT-only input)")
        return "gex"
    raise ValueError(f"Unsupported MuData modality {name!r}; expected GEX/RNA or protein/ADT")


def modalities(data, data_type="auto"):
    if is_mudata(data):
        if data_type != "auto":
            raise ValueError("--data-type applies only to .h5ad inputs")
        return [(name, adata, modality_kind(name, adata)) for name, adata in data.mod.items()]
    if "feature_types" in data.var and data.var["feature_types"].eq("Antibody Capture").any():
        raise ValueError("Input .h5ad contains Antibody Capture features; use .h5mu with a separate protein modality")
    kind = modality_kind("unknown", data, standalone=True) if data_type == "auto" else data_type
    return [(kind, data, kind)]


def ensure_matrix(adata, label):
    if adata.X is None or adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"{label}: non-empty .X matrix is required")


def validate_counts(adata, label):
    ensure_matrix(adata, label)
    if "log1p" in adata.uns:
        raise ValueError(f"{label}: expected raw counts in .X, but .uns['log1p'] indicates log-transformed data")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError(f"{label}: counts contain negative or non-finite entries")
    if not np.allclose(values, np.round(values), atol=1e-5, rtol=0):
        raise ValueError(f"{label}: .X contains non-integer values; supply raw counts in .X")


def extract_raw_protein(path):
    import scanpy as sc
    raw = sc.read_h5ad(path)
    if "feature_types" in raw.var:
        mask = (raw.var["feature_types"] == "Antibody Capture").to_numpy()
        raw = raw[:, mask].copy()
    if raw.n_vars == 0:
        raise ValueError("--protein-raw contains no antibody features")
    validate_counts(raw, "--protein-raw")
    return raw


def normalize(data, args, mods):
    import scanpy as sc
    import muon as mu

    raw = extract_raw_protein(args.protein_raw) if args.protein_raw else None
    for name, adata, kind in mods:
        validate_counts(adata, name)
        adata.layers["count"] = adata.X.copy()
        if kind == "gex":
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            adata.layers["normalized_gex"] = adata.X.copy()
            LOG.info("%s: total-count normalization + log1p", name)
        else:
            if raw is not None:
                # DSB requires a background/raw droplet matrix, not merely filtered counts.
                mu.prot.pp.dsb(adata, raw)
                method = "dsb"
            else:
                mu.prot.pp.clr(adata, flavor="seurat")
                method = "clr"
            adata.layers["normalized_pro"] = adata.X.copy()
            LOG.info("%s: %s protein normalization", name, method)


def select_hvg(data, args, mods):
    import scanpy as sc

    for name, adata, kind in mods:
        if kind == "protein":
            LOG.info("%s: HVG not applicable to protein; all protein features kept for PCA", name)
            continue
        ensure_matrix(adata, name)
        if "log1p" not in adata.uns:
            raise ValueError(f"{name}: expected log-normalized expression in .X with .uns['log1p']; run normalize first")
        top = args.n_top_genes or max(int(0.1 * adata.n_vars), 10)
        top = min(top, adata.n_vars)
        batch = args.batch_key if args.batch_key in adata.obs else None
        if args.batch_key and batch is None:
            LOG.warning("%s: batch key %r absent; selecting HVG without batch stratification", name, args.batch_key)
        sc.pp.highly_variable_genes(adata, n_top_genes=top, batch_key=batch, flavor="seurat", subset=False)
        LOG.info("%s: %d HVG selected", name, int(adata.var["highly_variable"].sum()))


def _plot_pca(adata, output, color):
    import scanpy as sc
    import matplotlib.pyplot as plt
    fig = sc.pl.pca(adata, color=color if color in adata.obs else None, use_raw=False, show=False, return_fig=True)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)


def _plot_variance(adata, output):
    import scanpy as sc
    import matplotlib.pyplot as plt
    sc.pl.pca_variance_ratio(adata, n_pcs=adata.obsm["X_pca"].shape[1], show=False)
    fig = plt.gcf()
    fig.savefig(output, bbox_inches="tight", dpi=300)
    plt.close(fig)


def _write_pca_scores(adata, output):
    coords = adata.obsm["X_pca"]
    scores = pd.DataFrame(coords, index=adata.obs_names, columns=[f"PC{i + 1}" for i in range(coords.shape[1])])
    metadata = adata.obs.copy(deep=False)
    for column, values in metadata.items():
        if pd.api.types.is_bool_dtype(values):
            metadata[column] = values.map({True: "TRUE", False: "FALSE"})
    table = pd.concat([scores, metadata], axis=1)
    try:
        table.to_csv(output, sep="\t", index_label="samples", na_rep="NA", quoting=csv.QUOTE_NONE, quotechar=None)
    except csv.Error as exc:
        raise ValueError("PCA scores require sample IDs and metadata without tabs or line breaks") from exc


def _path_for_mode(pattern, name, multiple):
    if pattern is None:
        return None
    if "{modality}" in pattern:
        return Path(pattern.replace("{modality}", name))
    if multiple:
        raise ValueError("A modality-specific output must contain '{modality}' when multiple modalities are processed")
    return Path(pattern)


def run_pca(data, args, mods):
    import scanpy as sc

    for name, adata, kind in mods:
        ensure_matrix(adata, name)
        use_hvg = kind == "gex"
        if use_hvg and ("highly_variable" not in adata.var or not adata.var["highly_variable"].any()):
            raise ValueError(f"{name}: highly_variable genes missing; run hvg first")
        n_features = int(adata.var["highly_variable"].sum()) if use_hvg else adata.n_vars
        n_comps = min(args.max_comps, adata.n_obs - 1, n_features - 1)
        if n_comps < 2:
            raise ValueError(f"{name}: insufficient observations/features for PCA")
        # Explicit selection avoids version-specific use_highly_variable behavior.
        sc.pp.pca(adata, n_comps=n_comps, mask_var="highly_variable" if use_hvg else None, svd_solver="arpack", random_state=args.seed)
        LOG.info("%s: PCA computed (%d components)", name, n_comps)
        p = _path_for_mode(args.plot, name, len(mods) > 1)
        if p:
            _plot_pca(adata, p, args.color)
        p = _path_for_mode(args.variance_plot, name, len(mods) > 1)
        if p:
            _plot_variance(adata, p)
        p = _path_for_mode(args.scores, name, len(mods) > 1)
        if p:
            _write_pca_scores(adata, p)
        p = _path_for_mode(args.explained_variance, name, len(mods) > 1)
        if p:
            ratios = np.asarray(adata.uns["pca"]["variance_ratio"], dtype=float)
            variance = pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(len(ratios))], "var": ratios * 100})
            variance.to_csv(p, sep="\t", index=False, float_format="%.1f", na_rep="NA")


def run_neighbors(data, args, mods):
    import scanpy as sc
    import muon as mu
    from sklearn.preprocessing import normalize as l2normalize

    if is_mudata(data) and len(mods) > 1:
        # WNN requires a shared set of cells. We do not subset silently.
        base = data.obs_names
        for name, adata, _ in mods:
            if len(base) != adata.n_obs or not base.equals(adata.obs_names):
                raise ValueError(f"WNN requires identical cells in identical order: {name} differs from MuData.obs_names")
        for name, adata, _ in mods:
            if "X_pca" not in adata.obsm:
                raise ValueError(f"{name}: missing obsm['X_pca']; run pca first")
            if adata.obsm["X_pca"].shape[1] < args.n_pcs:
                raise ValueError(f"{name}: PCA has fewer than --n-pcs={args.n_pcs} components")
            # Keep the original PCA unchanged; L2-normalized copies improve WNN comparability.
            adata.obsm["X_pca_wnn"] = l2normalize(adata.obsm["X_pca"][:, :args.n_pcs], axis=1)
            sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep="X_pca_wnn", random_state=args.seed)
        mu.pp.neighbors(data, n_neighbors=args.n_neighbors, random_state=args.seed)
        LOG.info("Computed multimodal WNN graph from %d modalities", len(mods))
    else:
        name, adata, _ = mods[0]
        if "X_pca" not in adata.obsm:
            raise ValueError(f"{name}: missing obsm['X_pca']; run pca first")
        if adata.obsm["X_pca"].shape[1] < args.n_pcs:
            raise ValueError(f"{name}: PCA has fewer than --n-pcs={args.n_pcs} components")
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs, use_rep="X_pca", random_state=args.seed)
        LOG.info("%s: neighbors graph calculated", name)
        # For single-modality MuData, place the graph at top level for unified UMAP.
        if is_mudata(data):
            import copy
            data.uns["neighbors"] = copy.deepcopy(adata.uns["neighbors"])
            data.obsp["distances"] = adata.obsp["distances"].copy()
            data.obsp["connectivities"] = adata.obsp["connectivities"].copy()


def _umap_owner(data, mods):
    return data if is_mudata(data) else mods[0][1]


def run_umap(data, args, mods):
    import scanpy as sc
    import muon as mu
    import matplotlib.pyplot as plt

    owner = _umap_owner(data, mods)
    if "neighbors" not in owner.uns:
        raise ValueError("Missing neighbors graph; run neighbors first")
    key = owner.uns["neighbors"].get("connectivities_key", "connectivities")
    if key not in owner.obsp:
        raise ValueError(f"Missing obsp[{key!r}] required by neighbors")
    if is_mudata(data):
        mu.tl.umap(data, min_dist=args.min_dist, random_state=args.seed)
        plot_umap = mu.pl.umap
    else:
        sc.tl.umap(owner, min_dist=args.min_dist, random_state=args.seed)
        plot_umap = sc.pl.umap
    LOG.info("UMAP calculated on %s graph", "integrated" if len(mods) > 1 else "single-modality")
    if args.csv:
        coords = pd.DataFrame(owner.obsm["X_umap"], index=owner.obs_names)
        coords.columns = [f"UMAP_{i+1}" for i in range(coords.shape[1])]
        coords.index.name = "cell_barcodes"
        coords.to_csv(args.csv)
    if args.plot:
        color = args.color if args.color in owner.obs else None
        fig = plot_umap(owner, color=color, use_raw=False, show=False, return_fig=True)
        fig.savefig(args.plot, bbox_inches="tight", dpi=300)
        plt.close(fig)
    if args.metadata_plot:
        from matplotlib.backends.backend_pdf import PdfPages
        cols = [c for c in owner.obs if c.startswith("meta_")]
        if not cols:
            LOG.warning("No meta_* columns present for metadata plot")
        else:
            with PdfPages(args.metadata_plot) as pdf:
                for col in cols:
                    fig = plot_umap(owner, color=col, use_raw=False, show=False, return_fig=True)
                    pdf.savefig(fig, bbox_inches="tight", dpi=300)
                    plt.close(fig)


def _gene_list(path):
    return [line.strip() for line in Path(path).read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]


def run_cell_cycle(data, args, mods):
    import scanpy as sc
    import matplotlib.pyplot as plt

    gex = [(name, adata) for name, adata, kind in mods if kind == "gex"]
    if not gex:
        raise ValueError("cell-cycle requires a GEX modality")
    if args.genes:
        all_genes = _gene_list(args.genes)
        s_genes, g2m_genes = all_genes[:args.s_count], all_genes[args.s_count:]
    else:
        s_genes, g2m_genes = _gene_list(args.s_genes), _gene_list(args.g2m_genes)
    if not s_genes or not g2m_genes:
        raise ValueError("Both S and G2M gene lists must be nonempty")
    for name, adata in gex:
        ensure_matrix(adata, name)
        names = adata.var[args.gene_column].astype(str) if args.gene_column in adata.var else adata.var_names.astype(str)
        # Temporary gene symbols for Scanpy; prevent duplicate names from breaking indexing.
        temp = adata.copy()
        temp.var_names = pd.Index(names)
        temp.var_names_make_unique()
        s_present = [x for x in s_genes if x in temp.var_names]
        g_present = [x for x in g2m_genes if x in temp.var_names]
        if not s_present or not g_present:
            raise ValueError(f"{name}: no matching genes in S or G2M list")
        if len(s_present) != len(s_genes) or len(g_present) != len(g2m_genes):
            LOG.warning("%s: matched %d/%d S and %d/%d G2M genes", name, len(s_present), len(s_genes), len(g_present), len(g2m_genes))
        sc.tl.score_genes_cell_cycle(temp, s_genes=s_present, g2m_genes=g_present, use_raw=False)
        for col in ("S_score", "G2M_score", "phase"):
            adata.obs[col] = temp.obs[col].to_numpy()
        LOG.info("%s: cell-cycle scores calculated", name)
        if args.plot or args.sample_plot:
            cycle_genes = list(dict.fromkeys(s_present + g_present))
            subset = temp[:, cycle_genes].copy()
            n_comps = min(50, subset.n_obs - 1, subset.n_vars - 1)
            if n_comps < 2:
                raise ValueError(f"{name}: insufficient cell-cycle genes/cells for PCA plotting")
            sc.pp.pca(subset, n_comps=n_comps, svd_solver="arpack", random_state=args.seed)
            for pattern, color in [(args.plot, "phase"), (args.sample_plot, args.color)]:
                dest = _path_for_mode(pattern, name, len(gex) > 1)
                if dest:
                    _plot_pca(subset, dest, color)


def write_dataset(data, path):
    if not path:
        return
    path = Path(path)
    if is_mudata(data):
        if path.suffix.lower() == ".h5mu":
            # Do not call update() unless necessary: it can alter global metadata.
            data.write_h5mu(path)
        elif path.suffix.lower() == ".h5ad" and len(data.mod) == 1:
            # Preserve integrated/top-level embeddings when extracting a single modality.
            adata = next(iter(data.mod.values())).copy()
            if "X_umap" in data.obsm:
                adata.obsm["X_umap"] = data.obsm["X_umap"].copy()
            if "neighbors" in data.uns:
                import copy
                adata.uns["neighbors"] = copy.deepcopy(data.uns["neighbors"])
                for key in ("distances", "connectivities"):
                    if key in data.obsp:
                        adata.obsp[key] = data.obsp[key].copy()
            adata.write_h5ad(path)
        else:
            raise ValueError("MuData requires .h5mu output; .h5ad is allowed only with exactly one modality")
    elif path.suffix.lower() == ".h5ad":
        data.write_h5ad(path)
    else:
        raise ValueError("AnnData output must be .h5ad")
    LOG.info("Wrote dataset: %s", path)


def positive_int(value):
    x = int(value)
    if x < 1:
        raise argparse.ArgumentTypeError("Must be a positive integer")
    return x


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version="preprocess 1.0.0")
    commands = p.add_subparsers(dest="command", required=True)
    for cmd in ("normalize", "hvg", "pca", "neighbors", "umap", "cell-cycle"):
        sub = commands.add_parser(cmd, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        sub.add_argument("-i", "--input", type=Path, required=True)
        sub.add_argument("-o", "--output", type=Path, help="Optional output .h5ad or .h5mu")
        sub.add_argument("--data-type", choices=("auto", "gex", "protein"), default="auto", help="Override modality detection for standalone .h5ad without Antibody Capture features")
        sub.add_argument("--seed", type=int, default=42)
        if cmd == "normalize":
            sub.add_argument("--protein-raw", type=Path, help="Unfiltered background .h5ad for protein DSB; without it use CLR")
        if cmd == "hvg":
            sub.add_argument("--batch-key", default="sample")
            sub.add_argument("--n-top-genes", type=positive_int, help="Default: max(10, 10%% of GEX features)")
        if cmd == "pca":
            sub.add_argument("--max-comps", type=positive_int, default=50)
            sub.add_argument("--plot", help="PCA PDF; for multimodal input include {modality}")
            sub.add_argument("--variance-plot", help="PCA explained-variance PDF; for multimodal include {modality}")
            sub.add_argument("--scores", help="Tab-separated PCA coordinates and observation metadata; for multimodal include {modality}")
            sub.add_argument("--explained-variance", help="Tab-separated explained variance percentages; for multimodal include {modality}")
            sub.add_argument("--color", default="sample")
        if cmd == "neighbors":
            sub.add_argument("--n-neighbors", type=positive_int, default=20)
            sub.add_argument("--n-pcs", type=positive_int, default=30)
        if cmd == "umap":
            sub.add_argument("--min-dist", type=float, default=0.1)
            sub.add_argument("--plot", help="UMAP PDF")
            sub.add_argument("--csv", help="UMAP coordinates CSV")
            sub.add_argument("--metadata-plot", help="PDF with one UMAP plot for each meta_* observation column")
            sub.add_argument("--color", default="sample")
        if cmd == "cell-cycle":
            genes = sub.add_mutually_exclusive_group(required=True)
            genes.add_argument("--genes", type=Path, help="Combined list: first --s-count genes are S, remaining G2M")
            genes.add_argument("--s-genes", type=Path, help="S-phase genes, paired with --g2m-genes")
            sub.add_argument("--g2m-genes", type=Path)
            sub.add_argument("--s-count", type=positive_int, default=43)
            sub.add_argument("--gene-column", default="gene_symbols")
            sub.add_argument("--plot", help="Cell-cycle PCA colored by phase")
            sub.add_argument("--sample-plot", help="Cell-cycle PCA colored by sample")
            sub.add_argument("--color", default="sample")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
    if args.command == "cell-cycle" and bool(args.s_genes) != bool(args.g2m_genes):
        parser.error("--s-genes and --g2m-genes must be supplied together")
    if args.command == "umap" and not 0 <= args.min_dist <= 1:
        parser.error("--min-dist must be between 0 and 1")
    try:
        data = read_dataset(args.input)
        mods = modalities(data, args.data_type)
        if not mods:
            raise ValueError("No modalities found")
        dispatch = {
            "normalize": normalize,
            "hvg": select_hvg,
            "pca": run_pca,
            "neighbors": run_neighbors,
            "umap": run_umap,
            "cell-cycle": run_cell_cycle,
        }
        dispatch[args.command](data, args, mods)
        write_dataset(data, args.output)
    except (ValueError, KeyError, FileNotFoundError) as exc:
        LOG.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
