#!/usr/bin/env python3
"""Single-cell preprocessing for AnnData (.h5ad) and MuData (.h5mu).

In multimodal MuData, GEX and protein are preprocessed separately, and
weighted nearest neighbors (WNN) integrates modality-specific PCA embeddings.
Only explicitly requested output files are written.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import csv
from html import escape
from inspect import signature
import json
import logging
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import sparse

LOG = logging.getLogger("preprocess")
GEX_NAMES = {"gex", "rna", "gene_expression", "gene-expression"}
PRO_NAMES = {"pro", "prot", "protein", "adt", "citeseq", "cite-seq"}
MODALITY_LABEL = {"gex": "GEX", "protein": "ADT"}
# Metric, report label, and the CLI destinations for its lower and upper bounds.
CELL_CRITERIA = {
    "gex": (
        ("total_counts", "Total counts", "min_counts", "max_counts"),
        ("n_genes_by_counts", "Detected genes", "min_genes", "max_genes"),
        ("pct_counts_mt", "Mitochondrial counts (%)", None, "max_mito_pct"),
    ),
    "protein": (
        ("total_counts", "Total counts", "min_adt_counts", None),
        ("n_genes_by_counts", "Detected antibodies", "min_antibodies", None),
    ),
}


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
    if "feature_types" in data.var:
        is_adt = data.var["feature_types"].eq("Antibody Capture").fillna(False)
        if is_adt.any() and (not is_adt.all() or data_type == "gex"):
            raise ValueError("Use a separate ADT modality for mixed features, or --data-type protein for ADT-only data")
    kind = modality_kind("unknown", data, standalone=True) if data_type == "auto" else data_type
    return [(kind, data, kind)]


def ensure_matrix(adata, label):
    if adata.X is None or adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"{label}: non-empty .X matrix is required")


def validate_counts(adata, label, layer=None):
    ensure_matrix(adata, label)
    if layer is None and ("log1p" in adata.uns or "normalized_pro" in adata.layers):
        raise ValueError(f"{label}: expected raw counts in .X or layers['count'], found normalized data")
    counts = adata.X if layer is None else adata.layers[layer]
    values = counts.data if sparse.issparse(counts) else np.asarray(counts)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError(f"{label}: counts contain negative or non-finite entries")
    if values.dtype.kind not in "biu" and not np.allclose(values, np.round(values), atol=1e-5, rtol=0):
        raise ValueError(f"{label}: {layer or '.X'} contains non-integer values; supply raw counts")
    return counts


def criterion_bounds(args, spec):
    _, _, lower, upper = spec
    return getattr(args, lower) if lower else None, getattr(args, upper) if upper else None


def sample_groups(data, mods, column):
    if column is None:
        return pd.Series("All observations", index=data.obs_names)
    samples = pd.Series(pd.NA, index=data.obs_names, dtype="string")
    frames = [data.obs]
    if is_mudata(data):
        frames.extend(adata.obs for _, adata, _ in mods)
    for obs in frames:
        if column not in obs:
            continue
        values = obs[column].astype("string").reindex(samples.index)
        both = samples.notna() & values.notna()
        if not samples[both].equals(values[both]):
            raise ValueError(f"Conflicting sample identifiers in observation column {column!r}")
        samples = samples.fillna(values)
    if samples.isna().any() or samples.str.strip().eq("").any():
        raise ValueError(f"--sample-column {column!r} requires a usable sample identifier for every cell")
    return samples


def qc_metrics(adata, kind, name):
    import scanpy as sc
    from anndata import AnnData
    from scipy.stats import median_abs_deviation

    layer = "count" if "count" in adata.layers else None
    counts = validate_counts(adata, name, layer)
    # Scanpy removes stored sparse zeros in place. Calculate on an independent matrix.
    work = AnnData(X=counts.copy() if sparse.issparse(counts) else counts)
    groups = []
    if kind == "gex":
        symbols = adata.var.get("gene_symbols", adata.var.get("gene_symbol", adata.var_names.to_series()))
        symbols = symbols.astype("string").str.upper()
        work.var["mt"] = symbols.str.startswith("MT-", na=False).to_numpy(dtype=bool)
        work.var["ribo"] = symbols.str.startswith(("RPS", "RPL"), na=False).to_numpy(dtype=bool)
        work.var["hb"] = symbols.str.match(r"^HB[^(P)]", na=False).to_numpy(dtype=bool)
        groups = ["mt", "ribo", "hb"]
        for group in groups:
            adata.var[group] = work.var[group].to_numpy()
    top = min(20 if kind == "gex" or adata.n_vars >= 20 else 2, adata.n_vars)
    with np.errstate(divide="ignore", invalid="ignore"):
        obs, var = sc.pp.calculate_qc_metrics(work, qc_vars=groups, percent_top=[top])
    for target, metrics in ((adata.obs, obs), (adata.var, var)):
        for metric in metrics:
            target[metric] = metrics[metric].to_numpy()
    soft = np.zeros(adata.n_obs, dtype=bool)
    for metric in ("log1p_total_counts", "log1p_n_genes_by_counts", f"pct_counts_in_top_{top}_genes"):
        values = adata.obs[metric]
        finite = values[np.isfinite(values)]
        flag = np.zeros(adata.n_obs, dtype=bool)
        if not finite.empty:
            flag = ((values - finite.median()).abs() > 5 * median_abs_deviation(finite)).to_numpy()
        adata.obs[f"qc_mad_{metric}"] = flag
        soft |= flag
    return soft


def run_qc(data, args, mods):
    gex = any(kind == "gex" for _, _, kind in mods)
    if args.doublet_column and not gex:
        raise ValueError("--doublet-column requires a GEX modality")
    if not data.obs_names.is_unique or any(not adata.obs_names.is_unique for _, adata, _ in mods):
        raise ValueError("QC requires unique cell identifiers on the observation axes")
    samples = sample_groups(data, mods, args.sample_column)
    aggregate = np.ones(data.n_obs, dtype=bool)
    for name, adata, kind in mods:
        doublet_pass = None
        if kind == "gex" and args.doublet_column:
            labels = adata.obs.get(args.doublet_column)
            if labels is None or labels.isna().any() or labels.astype(str).str.strip().eq("").any():
                raise ValueError(f"{name}: --doublet-column {args.doublet_column!r} requires a label for every cell")
            doublet_pass = labels.astype(str).ne(args.doublet_label)
        # Remove previous criterion flags so QC can be repeated with different bounds.
        old = [f"qc_pass_{spec[0]}" for spec in CELL_CRITERIA[kind]]
        if kind == "gex":
            old.append("qc_pass_doublet")
            adata.var.drop(columns="qc_pass_n_cells_by_counts", errors="ignore", inplace=True)
        old.extend(col for col in adata.obs if col.startswith("qc_mad_pct_counts_in_top_"))
        adata.obs.drop(columns=old, errors="ignore", inplace=True)
        soft = qc_metrics(adata, kind, name)
        passed = np.ones(adata.n_obs, dtype=bool)
        for spec in CELL_CRITERIA[kind]:
            metric = spec[0]
            lower, upper = criterion_bounds(args, spec)
            if lower is None and upper is None:
                continue
            values = adata.obs[metric]
            keep = values.notna()
            if lower is not None:
                keep &= values.ge(lower)
            if upper is not None:
                keep &= values.le(upper)
            adata.obs[f"qc_pass_{metric}"] = keep
            passed &= keep.to_numpy()
            if metric == "pct_counts_mt":
                soft |= ~keep.to_numpy()
        if doublet_pass is not None:
            adata.obs["qc_pass_doublet"] = doublet_pass
            passed &= doublet_pass.to_numpy()
        adata.obs["qc_pass"] = passed
        adata.obs["qc_soft_outlier"] = soft
        adata.var["qc_pass"] = True
        if kind == "gex" and args.min_cells is not None:
            adata.var["qc_pass_n_cells_by_counts"] = adata.var["n_cells_by_counts"].ge(args.min_cells)
            adata.var["qc_pass"] = adata.var["qc_pass_n_cells_by_counts"]
        aggregate &= adata.obs["qc_pass"].reindex(data.obs_names, fill_value=False).to_numpy()
    if is_mudata(data):
        data.obs["qc_pass"] = aggregate
    write_qc_reports(data, args, mods, samples)


def group_positions(samples, obs_names):
    labels = samples.reindex(obs_names)
    return labels.groupby(labels, sort=False).indices


def qc_summaries(data, args, mods, samples):
    summaries = {}
    unit, unit_title = ("ROIs", "ROIs") if args.report_observation_unit == "rois" else ("cells", "Cells")
    report_mods = {kind: [(name, adata) for name, adata, k in mods if k == kind]
                   for kind in CELL_CRITERIA}
    groups = {name: group_positions(samples, adata.obs_names) for name, adata, _ in mods}
    for sample, positions in group_positions(samples, data.obs_names).items():
        admitted = int(data.obs["qc_pass"].iloc[positions].sum())
        total = len(positions)
        row = {f"{unit_title} evaluated": total, "Admitted": admitted, "Excluded": total - admitted}
        bars = []
        for kind, spec_list in CELL_CRITERIA.items():
            # An explicit standalone data type also limits the report to that modality.
            if args.data_type != "auto" and args.data_type != kind:
                continue
            available = report_mods[kind]
            for name, adata in available or [(None, None)]:
                label = MODALITY_LABEL[kind]
                if len(available) > 1:
                    label = f"{label} ({name})"
                obs = None if adata is None else adata.obs.iloc[groups[name].get(sample, [])]
                criteria = [(*spec[:2], *criterion_bounds(args, spec), None) for spec in spec_list]
                if kind == "gex":
                    criteria.append(("doublet", "Doublets", None, None,
                                     args.doublet_label if args.doublet_column else None))
                    row[f"{label} genes: minimum detecting cells"] = (
                        "Unavailable" if obs is None else args.min_cells if args.min_cells is not None else "Inactive"
                    )
                for metric, title, lower, upper, excluded_label in criteria:
                    title = f"{label} {title.lower()}"
                    if lower is not None:
                        row[f"{title}: minimum"] = lower
                    if upper is not None:
                        row[f"{title}: maximum"] = upper
                    if excluded_label is not None:
                        row[f"{title}: excluded label"] = excluded_label
                    flag = f"qc_pass_{metric}"
                    if obs is None or obs.empty or flag not in obs:
                        status = "Unavailable" if obs is None or obs.empty else "Inactive"
                        row[f"{title}: {unit} evaluated"] = status
                        row[f"{title}: {unit} failing"] = status
                        if obs is None or flag not in obs:
                            continue
                    failed = int((~obs[flag]).sum())
                    if not obs.empty:
                        row[f"{title}: {unit} evaluated"] = len(obs)
                        row[f"{title}: {unit} failing"] = failed
                    details = []
                    if lower is not None:
                        details.append(f"Below minimum ({lower}): {int(obs[metric].lt(lower).sum())}")
                    if upper is not None:
                        details.append(f"Above maximum ({upper}): {int(obs[metric].gt(upper).sum())}")
                    if metric in obs and (missing := int(obs[metric].isna().sum())):
                        details.append(f"Undefined metric: {missing}")
                    bars.append((title, failed, "\n".join(details)))
                if obs is not None:
                    missing = total - len(obs)
                    row[f"{unit_title} missing from {label}"] = missing
                    if missing:
                        bars.append((f"{unit_title} missing from {label}", missing, ""))
        summaries[str(sample)] = (row, bars)
    return summaries


def qc_interactive_html(summaries, title):
    parts = [f'<section class="sc-qc-report" aria-label="{title}">', """<style>
.sc-qc-report details { margin: .6em 0; }
.sc-qc-report summary { cursor: pointer; font-weight: bold; }
.sc-qc-report .sc-qc-track { display: flex; height: 1em; background: #dce1e6; max-width: 42em; }
.sc-qc-report .sc-qc-bar { display: block; height: 100%; background: #a43d38; }
.sc-qc-report .sc-qc-admitted { background: #197764; }
</style>"""]
    for i, (sample, (row, bars)) in enumerate(summaries.items()):
        admitted, excluded = row["Admitted"], row["Excluded"]
        total = admitted + excluded
        outcome = f"Admitted: {admitted} / Excluded: {excluded}"
        parts.extend([
            f'<details{" open" if i == 0 else ""}><summary>{escape(sample)}</summary>',
            f'<p>{outcome}</p><div class="sc-qc-track" role="img" aria-label="{outcome}">',
            f'<span class="sc-qc-bar sc-qc-admitted" style="width:{100 * admitted / total:.8f}%" '
            f'title="Admitted: {admitted}"></span>',
            f'<span class="sc-qc-bar" style="width:{100 * excluded / total:.8f}%" '
            f'title="Excluded: {excluded}"></span></div>',
        ])
        if bars:
            parts.append("<p>Criterion failures</p>")
        for title, failed, details in bars:
            caption = f"{title}: {failed}"
            hint = escape(f"{caption}\n{details}".rstrip())
            parts.append(
                f'<div title="{hint}">{escape(caption)}<div class="sc-qc-track" tabindex="0" '
                f'role="img" aria-label="{hint}"><span class="sc-qc-bar" '
                f'style="width:{100 * failed / total:.8f}%"></span></div></div>'
            )
        parts.append("</details>")
    return "\n".join([*parts, "</section>"])


def qc_distribution_pdf(path, args, mods, samples, kind):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.stats import gaussian_kde

    metrics = [(spec[0], spec[1], criterion_bounds(args, spec)) for spec in CELL_CRITERIA[kind]]
    if kind == "gex":
        metrics.append(("pct_counts_ribo", "Ribosomal counts (%)", (None, None)))
    report_title = f"{MODALITY_LABEL[kind]} QC"
    if args.report_observation_unit == "rois":
        report_title = f"ROI pseudobulk {report_title}"
    groups = {name: group_positions(samples, adata.obs_names) for name, adata in mods}
    with PdfPages(path) as pdf:
        for sample in samples.unique():
            fig, axes = plt.subplots(len(mods), len(metrics), figsize=(3.5 * len(metrics), 3 * len(mods)),
                                     squeeze=False, layout="constrained")
            fig.suptitle(f"{report_title} – {sample}")
            for (name, adata), axes_row in zip(mods, axes):
                obs = adata.obs.iloc[groups[name].get(sample, [])]
                for ax, (metric, title, bounds) in zip(axes_row, metrics):
                    values = obs[metric].dropna().to_numpy()
                    ax.set_title(f"{name}: {title}" if len(mods) > 1 else title)
                    if values.size:
                        ax.hist(values, bins=30, density=True, alpha=0.5)
                        if values.min() < values.max():
                            grid = np.linspace(values.min(), values.max(), 128)
                            ax.plot(grid, gaussian_kde(values)(grid))
                    else:
                        ax.text(0.5, 0.5, "No values", transform=ax.transAxes, ha="center")
                    for bound in bounds:
                        if bound is not None:
                            ax.axvline(bound, color="firebrick", linestyle="--")
                    ax.set_ylabel("Density")
            pdf.savefig(fig)
            plt.close(fig)


def qc_violin(path, mods, samples):
    import matplotlib.pyplot as plt

    labels = samples.unique()
    fig, axes = plt.subplots(1, len(mods), figsize=(max(5, len(labels)) * len(mods), 4),
                             squeeze=False, layout="constrained")
    for (name, adata, kind), ax in zip(mods, axes[0]):
        groups = group_positions(samples, adata.obs_names)
        for position, sample in enumerate(labels, 1):
            values = adata.obs["total_counts"].iloc[groups.get(sample, [])].to_numpy()
            if values.size:
                ax.violinplot(values, positions=[position], showmedians=True)
        ax.set_title(f"{MODALITY_LABEL[kind]} ({name})")
        ax.set_ylabel("Raw total counts per cell")
        ax.set_xticks(range(1, len(labels) + 1), labels, rotation=45, ha="right")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_qc_reports(data, args, mods, samples):
    for kind, path in (("gex", args.gex_pdf), ("protein", args.adt_pdf)):
        if path:
            selected = [(name, adata) for name, adata, k in mods if k == kind]
            if not selected:
                raise ValueError(f"Cannot write a {MODALITY_LABEL[kind]} PDF without that modality")
            qc_distribution_pdf(path, args, selected, samples, kind)
    if args.multiqc_violin:
        qc_violin(args.multiqc_violin, mods, samples)
    if not (args.html or args.multiqc_html or args.multiqc_table):
        return
    summaries = qc_summaries(data, args, mods, samples)
    if args.report_observation_unit == "rois":
        report_id, title = "roi_pseudobulk_qc", "ROI pseudobulk QC"
    else:
        report_id, title = "single_cell_qc", "Single-cell QC"
    if args.html or args.multiqc_html:
        view = qc_interactive_html(summaries, title)
        if args.html:
            args.html.write_text(f'<!doctype html><html lang="en"><meta charset="utf-8">'
                                 f'<title>{title}</title><body>{view}</body></html>\n', encoding="utf-8")
        if args.multiqc_html:
            args.multiqc_html.write_text(
                f"<!--\nid: {report_id}_outcomes\nsection_name: {title} outcomes\n-->\n"
                + view + "\n", encoding="utf-8",
            )
    if args.multiqc_table:
        table = {
            "id": f"{report_id}_summary",
            "section_name": f"{title} summary",
            "plot_type": "table",
            "pconfig": {"id": f"{report_id}_summary_table", "title": f"{title} summary"},
            "data": {sample: row for sample, (row, _) in summaries.items()},
        }
        args.multiqc_table.write_text(json.dumps(table, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


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


def axis_mask(frame, key="qc_pass"):
    if key not in frame:
        return np.ones(len(frame), dtype=bool)
    mask = frame[key]
    if not pd.api.types.is_bool_dtype(mask) or mask.isna().any():
        raise ValueError(f"{key!r} must contain a boolean for every element of its axis")
    return mask.to_numpy(dtype=bool)


def qc_rows(data, adata):
    if data is not adata and "qc_pass" in data.obs:
        shared = pd.Series(axis_mask(data.obs), index=data.obs_names)
        return shared.reindex(adata.obs_names, fill_value=False).to_numpy(dtype=bool)
    return axis_mask(adata.obs)


def check_normalization(adata, name, rows, cols):
    key = "normalization_pass"
    if key in adata.obs or key in adata.var:
        current = (key in adata.obs and key in adata.var
                   and np.array_equal(rows, axis_mask(adata.obs, key))
                   and np.array_equal(cols, axis_mask(adata.var, key)))
    else:
        normalized = "log1p" in adata.uns or any(layer in adata.layers for layer in ("normalized_gex", "normalized_pro"))
        current = not normalized or (rows.all() and cols.all())
    if not current:
        raise ValueError(f"{name}: normalization does not match current QC; rerun normalize")


def analysis_subset(data, adata, name):
    from anndata import AnnData

    ensure_matrix(adata, name)
    rows, cols = qc_rows(data, adata), axis_mask(adata.var)
    check_normalization(adata, name, rows, cols)
    if not rows.any() or not cols.any():
        raise ValueError(f"{name}: no cells or features passed QC; adjust the QC criteria")
    if rows.all() and cols.all():
        return adata, rows, cols
    view = adata[rows, cols]
    work = AnnData(X=view.X, obs=view.obs.copy(), var=view.var.copy())
    if "log1p" in adata.uns:
        work.uns["log1p"] = deepcopy(adata.uns["log1p"])
    return work, rows, cols


def restore_matrix(values, rows, cols, shape):
    if np.array_equal(rows, np.arange(shape[0])) and np.array_equal(cols, np.arange(shape[1])):
        return values
    if sparse.issparse(values):
        coo = values.tocoo()
        return sparse.csr_matrix((coo.data, (rows[coo.row], cols[coo.col])), shape=shape)
    result = np.zeros(shape, dtype=values.dtype)
    result[np.ix_(rows, cols)] = values
    return result


def restore_coordinates(values, rows, size):
    if np.array_equal(rows, np.arange(size)):
        return values
    result = np.full((size, values.shape[1]), np.nan, dtype=values.dtype)
    result[rows] = values
    return result


def normalize(data, args, mods):
    import scanpy as sc
    from anndata import AnnData

    raw = extract_raw_protein(args.protein_raw) if args.protein_raw else None
    for name, adata, kind in mods:
        layer = "count" if "count" in adata.layers else None
        counts = validate_counts(adata, name, layer)
        rows, cols = qc_rows(data, adata), axis_mask(adata.var)
        if not rows.any() or not cols.any():
            raise ValueError(f"{name}: normalization requires cells and features that pass QC; adjust the QC criteria")
        selected = counts.copy() if rows.all() and cols.all() else counts[rows, :][:, cols]
        work = AnnData(X=selected, obs=adata.obs.loc[rows, []], var=adata.var.loc[cols, []])
        if kind == "gex":
            sc.pp.normalize_total(work)
            sc.pp.log1p(work)
            adata.uns["log1p"] = work.uns["log1p"]
            normalized_layer = "normalized_gex"
            LOG.info("%s: total-count normalization + log1p", name)
        else:
            import muon as mu

            if raw is not None:
                if (not raw.var_names.is_unique or not work.var_names.is_unique
                        or not work.var_names.isin(raw.var_names).all()):
                    raise ValueError(f"{name}: --protein-raw must contain each admitted antibody with a matching unique identifier")
                # Excluded cells still belong to the cell population, never to the background.
                background = raw[~raw.obs_names.isin(adata.obs_names), work.var_names]
                if background.n_obs < 2 or work.n_vars < 2:
                    raise ValueError(f"{name}: DSB requires at least two background droplets and two admitted antibodies")
                mu.prot.pp.dsb(work, background)
                method = "dsb"
            else:
                # Older Muon releases implement the same Seurat transform without a flavor option.
                work.X = work.X.astype(np.result_type(work.X.dtype, np.float32), copy=False)
                clr = mu.prot.pp.clr
                clr(work, **({"flavor": "seurat"} if "flavor" in signature(clr).parameters else {}))
                method = "clr"
            normalized_layer = "normalized_pro"
            LOG.info("%s: %s protein normalization", name, method)
        if layer is None:
            adata.layers["count"] = counts
        adata.X = restore_matrix(work.X, np.flatnonzero(rows), np.flatnonzero(cols), adata.shape)
        adata.layers[normalized_layer] = adata.X.copy()
        adata.obs["normalization_pass"] = rows
        adata.var["normalization_pass"] = cols


def check_fit_masks(adata, name, rows, cols, command):
    key = f"{command}_pass"
    if key in adata.obs and (not np.array_equal(rows, axis_mask(adata.obs, key))
                             or not np.array_equal(cols, axis_mask(adata.var, key))):
        raise ValueError(f"{name}: {command} was computed with different QC masks; rerun {command}")


def select_hvg(data, args, mods):
    import scanpy as sc

    for name, adata, kind in mods:
        if kind == "protein":
            LOG.info("%s: HVG not applicable to protein; all protein features kept for PCA", name)
            continue
        work, rows, cols = analysis_subset(data, adata, name)
        if "log1p" not in adata.uns:
            raise ValueError(f"{name}: expected log-normalized expression in .X with .uns['log1p']; run normalize first")
        if work.n_vars < 2 or work.n_obs < 2:
            raise ValueError(f"{name}: HVG selection requires at least two admitted cells and genes; adjust QC")
        top = args.n_top_genes or max(int(0.1 * work.n_vars), 10)
        top = min(top, work.n_vars)
        batch = args.batch_key if args.batch_key in work.obs else None
        if args.batch_key and batch is None:
            LOG.warning("%s: batch key %r absent; selecting HVG without batch stratification", name, args.batch_key)
        sc.pp.highly_variable_genes(work, n_top_genes=top, batch_key=batch, flavor="seurat", subset=False)
        if work is not adata:
            fields = ["highly_variable", "means", "dispersions", "dispersions_norm"]
            if batch:
                fields.extend(["highly_variable_nbatches", "highly_variable_intersection"])
            for field in fields:
                values = pd.Series(work.var[field].to_numpy(), index=np.flatnonzero(cols))
                fill = False if pd.api.types.is_bool_dtype(values) else np.nan
                adata.var[field] = values.reindex(range(adata.n_vars), fill_value=fill).to_numpy()
            adata.uns["hvg"] = work.uns["hvg"]
        adata.obs["hvg_pass"] = rows
        adata.var["hvg_pass"] = cols
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
        work, rows, cols = analysis_subset(data, adata, name)
        use_hvg = kind == "gex"
        if use_hvg:
            check_fit_masks(adata, name, rows, cols, "hvg")
            if "highly_variable" not in work.var or not work.var["highly_variable"].any():
                raise ValueError(f"{name}: no admitted highly variable genes; run hvg first or adjust QC")
        n_features = int(work.var["highly_variable"].sum()) if use_hvg else work.n_vars
        n_comps = min(args.max_comps, work.n_obs - 1, n_features - 1)
        if n_comps < 2:
            raise ValueError(f"{name}: PCA requires at least three admitted cells and PCA features, and --max-comps >= 2")
        # Explicit selection avoids version-specific use_highly_variable behavior.
        sc.pp.pca(work, n_comps=n_comps, mask_var="highly_variable" if use_hvg else None, svd_solver="arpack", random_state=args.seed)
        if work is not adata:
            adata.obsm["X_pca"] = restore_coordinates(work.obsm["X_pca"], np.flatnonzero(rows), adata.n_obs)
            adata.varm["PCs"] = restore_matrix(work.varm["PCs"], np.flatnonzero(cols), np.arange(n_comps),
                                              (adata.n_vars, n_comps))
            adata.uns["pca"] = work.uns["pca"]
        adata.obs["pca_pass"] = rows
        adata.var["pca_pass"] = cols
        LOG.info("%s: PCA computed (%d components)", name, n_comps)
        p = _path_for_mode(args.plot, name, len(mods) > 1)
        if p:
            _plot_pca(work, p, args.color)
        p = _path_for_mode(args.variance_plot, name, len(mods) > 1)
        if p:
            _plot_variance(work, p)
        p = _path_for_mode(args.scores, name, len(mods) > 1)
        if p:
            _write_pca_scores(adata, p)
        p = _path_for_mode(args.plot_scores, name, len(mods) > 1)
        if p:
            _write_pca_scores(work, p)
        p = _path_for_mode(args.explained_variance, name, len(mods) > 1)
        if p:
            ratios = np.asarray(adata.uns["pca"]["variance_ratio"], dtype=float)
            variance = pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(len(ratios))], "var": ratios * 100})
            variance.to_csv(p, sep="\t", index=False, float_format="%.1f", na_rep="NA")


def store_graph(target, graph):
    rows = target.obs_names.get_indexer(graph.obs_names)
    neighbors = graph.uns["neighbors"]
    for key in ("connectivities_key", "distances_key"):
        edge_key = neighbors[key]
        target.obsp[edge_key] = restore_matrix(graph.obsp[edge_key], rows, rows, (target.n_obs, target.n_obs))
    target.uns["neighbors"] = deepcopy(neighbors)
    target.obs["neighbors_pass"] = target.obs_names.isin(graph.obs_names)


def run_neighbors(data, args, mods):
    import scanpy as sc
    from anndata import AnnData

    integrated = len(mods) > 1
    has_qc = "qc_pass" in data.obs
    if integrated:
        ids = data.obs_names[axis_mask(data.obs)]
    else:
        adata = mods[0][1]
        ids = adata.obs_names[qc_rows(data, adata)]
    if args.n_neighbors < 2 or len(ids) < 2:
        raise ValueError("Neighbors requires at least two admitted cells and --n-neighbors >= 2")
    if (integrated and len(ids) <= args.n_neighbors) or (has_qc and len(ids) < args.n_neighbors):
        raise ValueError(f"Only {len(ids)} cells are admitted for --n-neighbors={args.n_neighbors}; lower --n-neighbors or adjust QC")
    if integrated:
        import muon as mu
        from mudata import MuData
        from sklearn.preprocessing import normalize as l2normalize

    graphs = {}
    for name, adata, _ in mods:
        if integrated and not has_qc and not data.obs_names.equals(adata.obs_names):
            raise ValueError(f"WNN requires identical cells in identical order without QC: {name} differs from MuData.obs_names")
        rows = adata.obs_names.get_indexer(ids)
        if np.any(rows < 0):
            raise ValueError(f"{name}: some admitted cells are absent; rerun qc")
        admitted, cols = qc_rows(data, adata), axis_mask(adata.var)
        check_normalization(adata, name, admitted, cols)
        check_fit_masks(adata, name, admitted, cols, "pca")
        if "X_pca" not in adata.obsm:
            raise ValueError(f"{name}: missing obsm['X_pca']; run pca first")
        if adata.obsm["X_pca"].shape[1] < args.n_pcs:
            raise ValueError(f"{name}: PCA has fewer than --n-pcs={args.n_pcs} components")
        pcs = adata.obsm["X_pca"][rows, :args.n_pcs]
        if not np.isfinite(pcs).all():
            raise ValueError(f"{name}: missing PCA coordinates for admitted cells; rerun pca")
        graph = AnnData(obs=pd.DataFrame(index=ids))
        if integrated:
            rep = "X_pca_wnn"
            graph.obsm[rep] = l2normalize(pcs, axis=1)
            adata.obsm[rep] = restore_coordinates(graph.obsm[rep], rows, adata.n_obs)
        else:
            rep = "X_pca"
            graph.obsm[rep] = pcs
        sc.pp.neighbors(graph, n_neighbors=args.n_neighbors, n_pcs=None if integrated else args.n_pcs,
                        use_rep=rep, random_state=args.seed)
        store_graph(adata, graph)
        adata.var["neighbors_pass"] = cols
        graphs[name] = graph
    if integrated:
        graph = MuData(graphs)
        mu.pp.neighbors(graph, n_neighbors=args.n_neighbors, random_state=args.seed,
                        n_bandwidth_neighbors=min(20, len(ids) - 1),
                        n_multineighbors=min(200, len(ids) - 1))
        for name, _, _ in mods:
            key = f"{name}:mod_weight"
            if key in graph.obs:
                data.obs[key] = graph.obs[key].reindex(data.obs_names)
        LOG.info("Computed multimodal WNN graph from %d modalities", len(mods))
    if is_mudata(data):
        store_graph(data, graph)


def run_umap(data, args, mods):
    import scanpy as sc
    import matplotlib.pyplot as plt
    from anndata import AnnData

    owner = data
    integrated = len(mods) > 1
    if "neighbors" not in owner.uns:
        raise ValueError("Missing neighbors graph; run neighbors first")
    for name, adata, _ in mods:
        rows, cols = qc_rows(data, adata), axis_mask(adata.var)
        check_normalization(adata, name, rows, cols)
        check_fit_masks(adata, name, rows, cols, "pca")
        check_fit_masks(adata, name, rows, cols, "neighbors")
    graph_rows = axis_mask(owner.obs, "neighbors_pass")
    if not np.array_equal(graph_rows, axis_mask(owner.obs)):
        raise ValueError("Neighbors graph does not match current QC; rerun neighbors")
    ids = owner.obs_names[graph_rows]
    if len(ids) < 4:
        raise ValueError("UMAP requires at least four admitted cells; adjust QC and rerun neighbors")
    neighbors = owner.uns["neighbors"]
    if integrated:
        import muon as mu
        from mudata import MuData

        reps = {}
        for name, adata, _ in mods:
            rep = neighbors["params"]["use_rep"][name]
            reps[name] = AnnData(obs=pd.DataFrame(index=ids), obsm={rep: adata.obsm[rep][adata.obs_names.get_indexer(ids)]})
        work = MuData(reps)
        work.obs = owner.obs.loc[ids].copy()
        umap, plot_umap = mu.tl.umap, mu.pl.umap
    else:
        adata = mods[0][1]
        rep = neighbors["params"].get("use_rep")
        if owner is adata and graph_rows.all():
            work = adata
        else:
            work = AnnData(obs=owner.obs.loc[ids].copy())
            work.obsm[rep] = adata.obsm[rep][adata.obs_names.get_indexer(ids)]
        umap, plot_umap = sc.tl.umap, sc.pl.umap
    if work is not owner:
        work.uns["neighbors"] = deepcopy(neighbors)
        for key in ("distances_key", "connectivities_key"):
            edge_key = neighbors[key]
            work.obsp[edge_key] = owner.obsp[edge_key][graph_rows, :][:, graph_rows]
    umap(work, min_dist=args.min_dist, random_state=args.seed)
    if work is not owner:
        owner.obsm["X_umap"] = restore_coordinates(work.obsm["X_umap"], np.flatnonzero(graph_rows), owner.n_obs)
        owner.uns["umap"] = work.uns["umap"]
    LOG.info("UMAP calculated on %s graph", "integrated" if integrated else "single-modality")
    if args.csv:
        coords = pd.DataFrame(owner.obsm["X_umap"], index=owner.obs_names)
        coords.columns = [f"UMAP_{i+1}" for i in range(coords.shape[1])]
        coords.index.name = "cell_barcodes"
        coords.to_csv(args.csv)
    if args.plot:
        color = args.color if args.color in work.obs else None
        fig = plot_umap(work, color=color, use_raw=False, show=False, return_fig=True)
        fig.savefig(args.plot, bbox_inches="tight", dpi=300)
        plt.close(fig)
    if args.metadata_plot:
        from matplotlib.backends.backend_pdf import PdfPages
        cols = [c for c in work.obs if c.startswith("meta_")]
        if not cols:
            LOG.warning("No meta_* columns present for metadata plot")
        else:
            with PdfPages(args.metadata_plot) as pdf:
                for col in cols:
                    fig = plot_umap(work, color=col, use_raw=False, show=False, return_fig=True)
                    pdf.savefig(fig, bbox_inches="tight", dpi=300)
                    plt.close(fig)


def _gene_list(path):
    return [line.strip() for line in Path(path).read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]


def run_cell_cycle(data, args, mods):
    import scanpy as sc
    import matplotlib.pyplot as plt
    from anndata import AnnData

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
        work, rows, _ = analysis_subset(data, adata, name)
        if work.n_vars < 2:
            raise ValueError(f"{name}: cell-cycle scoring requires at least two admitted genes; adjust QC")
        names = work.var[args.gene_column].astype(str) if args.gene_column in work.var else work.var_names.astype(str)
        # Temporary gene symbols for Scanpy; prevent duplicate names from breaking indexing.
        temp = AnnData(X=work.X, obs=work.obs.copy(), var=work.var.copy())
        if "log1p" in work.uns:
            temp.uns["log1p"] = deepcopy(work.uns["log1p"])
        temp.var_names = pd.Index(names)
        temp.var_names_make_unique()
        s_present = [x for x in s_genes if x in temp.var_names]
        g_present = [x for x in g2m_genes if x in temp.var_names]
        if not s_present or not g_present:
            raise ValueError(f"{name}: no admitted genes match the S or G2M list; check gene symbols and QC")
        if len(s_present) != len(s_genes) or len(g_present) != len(g2m_genes):
            LOG.warning("%s: matched %d/%d S and %d/%d G2M genes", name, len(s_present), len(s_genes), len(g_present), len(g2m_genes))
        try:
            sc.tl.score_genes_cell_cycle(temp, s_genes=s_present, g2m_genes=g_present, use_raw=False)
        except RuntimeError as exc:
            raise ValueError(f"{name}: cell-cycle scoring failed on admitted genes: {exc}") from exc
        for col in ("S_score", "G2M_score"):
            scores = np.full(adata.n_obs, np.nan)
            scores[rows] = temp.obs[col].to_numpy()
            adata.obs[col] = scores
        phase = np.full(adata.n_obs, None, dtype=object)
        phase[rows] = temp.obs["phase"].to_numpy()
        adata.obs["phase"] = pd.Categorical(phase)
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
            rows = data.obs_names.get_indexer(adata.obs_names)
            if "X_umap" in data.obsm:
                adata.obsm["X_umap"] = data.obsm["X_umap"][rows]
            if "neighbors" in data.uns:
                neighbors = data.uns["neighbors"]
                adata.uns["neighbors"] = deepcopy(neighbors)
                for key in ("distances_key", "connectivities_key"):
                    edge_key = neighbors.get(key)
                    if edge_key in data.obsp:
                        adata.obsp[edge_key] = data.obsp[edge_key][rows, :][:, rows]
                if "neighbors_pass" in data.obs:
                    adata.obs["neighbors_pass"] = axis_mask(data.obs, "neighbors_pass")[rows]
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


def nonnegative_int(value):
    x = int(value)
    if x < 0:
        raise argparse.ArgumentTypeError("Must be a non-negative integer")
    return x


def percentage(value):
    x = float(value)
    if not 0 <= x <= 100:
        raise argparse.ArgumentTypeError("Must be between 0 and 100")
    return x


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version="preprocess 1.0.0")
    commands = p.add_subparsers(dest="command", required=True)
    for cmd in ("qc", "normalize", "hvg", "pca", "neighbors", "umap", "cell-cycle"):
        sub = commands.add_parser(cmd, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        sub.add_argument("-i", "--input", type=Path, required=True)
        sub.add_argument("-o", "--output", type=Path, help="Optional output .h5ad or .h5mu")
        sub.add_argument("--data-type", choices=("auto", "gex", "protein"), default="auto", help="Override modality detection for standalone .h5ad")
        sub.add_argument("--seed", type=int, default=42)
        if cmd == "qc":
            sub.add_argument("--min-counts", type=nonnegative_int, help="Minimum total GEX counts per cell")
            sub.add_argument("--max-counts", type=nonnegative_int, help="Maximum total GEX counts per cell")
            sub.add_argument("--min-genes", type=nonnegative_int, help="Minimum detected genes per cell")
            sub.add_argument("--max-genes", type=nonnegative_int, help="Maximum detected genes per cell")
            sub.add_argument("--max-mito-pct", type=percentage, help="Maximum percentage of mitochondrial counts per cell")
            sub.add_argument("--min-cells", type=nonnegative_int, help="Minimum number of cells detecting each GEX gene")
            sub.add_argument("--min-adt-counts", type=nonnegative_int, help="Minimum total ADT counts per cell")
            sub.add_argument("--min-antibodies", type=nonnegative_int, help="Minimum detected antibodies per cell")
            sub.add_argument("--doublet-column", help="Existing GEX observation column with doublet labels")
            sub.add_argument("--doublet-label", default="doublet", help="Label to exclude when --doublet-column is set")
            sub.add_argument("--sample-column", help="Observation column with sample identifiers for reports")
            sub.add_argument("--report-observation-unit", choices=("cells", "rois"), default="cells",
                             help="Observation labels used in QC reports")
            sub.add_argument("--gex-pdf", type=Path, help="Per-sample GEX distributions PDF")
            sub.add_argument("--adt-pdf", type=Path, help="Per-sample ADT distributions PDF")
            sub.add_argument("--html", type=Path, help="Standalone interactive QC report")
            sub.add_argument("--multiqc-html", type=Path, help="Interactive MultiQC fragment (*_mqc.html)")
            sub.add_argument("--multiqc-violin", type=Path, help="MultiQC violin plots (*_mqc.png)")
            sub.add_argument("--multiqc-table", type=Path, help="MultiQC summary table (*_mqc.json)")
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
            sub.add_argument("--plot-scores", help="Tab-separated PCA scores for computed observations; for multimodal include {modality}")
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
    if args.command == "qc":
        for metric in ("counts", "genes"):
            lower, upper = getattr(args, f"min_{metric}"), getattr(args, f"max_{metric}")
            if lower is not None and upper is not None and lower > upper:
                parser.error(f"--min-{metric} must be <= --max-{metric}")
    try:
        data = read_dataset(args.input)
        mods = modalities(data, args.data_type)
        if not mods:
            raise ValueError("No modalities found")
        dispatch = {
            "qc": run_qc,
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
