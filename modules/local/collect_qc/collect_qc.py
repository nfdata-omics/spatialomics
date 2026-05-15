
import importlib
import importlib.metadata
import argparse
import sys
import yaml

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
pio.templates.default = "plotly_white"

def collect_qc_metrics_files(metrics_files, output_file):
    """
    Collect QC metrics from metrics summary files.

    Parameters
    ----------
    metrics_files : list of str
        List of file paths to metrics summary CSV files.

    output_file : str
        Path to output combined QC metrics CSV file.

    Returns
    -------
    None
        Saves the combined annotated metrics to a unique CSV file.
    """

    # Read all the metrics csv file with Pandas and store them in a list of DataFrames
    dfs = []
    for file in metrics_files:
        df = pd.read_csv(file)
        dfs.append(df)

    # Concatenate all DataFrames
    combined = pd.concat(dfs, ignore_index=True)

    # Save to a single CSV file
    combined.to_csv(output_file, index=False)


def plot_qc_distributions(annotated_obs_files, output_file):
    """
    Plot QC distributions from annotated obs summary files.

    Parameters
    ----------
    annotated_obs_files : list of str
        List of file paths to annotated obs summary CSV files.

    output_file : str
        Path to output HTML file containing QC distribution plots.

    Returns
    -------
    None
        Saves the QC distribution plots to an HTML file.
    """

    # Read all the annotated obs csv file with Pandas and store them in a list of DataFrames
    dfs = []
    for file in annotated_obs_files:
        df = pd.read_csv(file)
        dfs.append(df)

    # Concatenate all DataFrames
    combined = pd.concat(dfs, ignore_index=True)

    metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("<html><head><meta charset='utf-8'><title>QC plots</title></head><body>\n")

        for i, metric in enumerate(metrics):

            fig = qc_histogram_figure(combined, metric)
            if i == 0:
                f.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
            else:
                f.write(fig.to_html(full_html=False, include_plotlyjs=False))

        f.write("</body></html>\n")


def nice_bin_width(raw_width, minimum_width=1.0):
    """Round a raw bin width up to a readable 1/2/5 x 10^n width."""
    raw_width = max(float(raw_width), minimum_width)
    exponent = np.floor(np.log10(raw_width))
    scale = 10**exponent
    fraction = raw_width / scale

    if fraction <= 1:
        nice_fraction = 1
    elif fraction <= 2:
        nice_fraction = 2
    elif fraction <= 5:
        nice_fraction = 5
    else:
        nice_fraction = 10

    return nice_fraction * scale


def histogram_edges(values, max_bins=100, minimum_bin_width=1.0):
    """Build shared histogram edges with a capped bin count and minimum bin width."""
    values = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    if values.size == 0:
        return None

    min_value = float(values.min())
    max_value = float(values.max())
    if min_value == max_value:
        half_width = minimum_bin_width / 2
        return np.array([min_value - half_width, min_value + half_width], dtype=float)

    target_bins = min(max_bins, max(10, int(np.sqrt(values.size))))
    bin_width = nice_bin_width((max_value - min_value) / target_bins, minimum_bin_width)

    start = np.floor(min_value / bin_width) * bin_width
    stop = np.ceil(max_value / bin_width) * bin_width
    edges = np.arange(start, stop + bin_width, bin_width, dtype=float)

    if edges.size > max_bins + 1:
        bin_width = nice_bin_width((max_value - min_value) / max_bins, minimum_bin_width)
        start = np.floor(min_value / bin_width) * bin_width
        stop = np.ceil(max_value / bin_width) * bin_width
        edges = np.arange(start, stop + bin_width, bin_width, dtype=float)

    if edges.size < 2:
        edges = np.array([min_value - minimum_bin_width / 2, min_value + minimum_bin_width / 2])

    return edges


def qc_histogram_figure(combined, metric):
    """Create an interactive normalized histogram without embedding raw observations."""
    fig = go.Figure()
    palette = px.colors.qualitative.Plotly
    edges = histogram_edges(combined[metric])
    if edges is None:
        fig.update_layout(title=metric)
        return fig

    centers = (edges[:-1] + edges[1:]) / 2
    widths = edges[1:] - edges[:-1]

    for index, (sample, sample_df) in enumerate(combined.groupby("sample", sort=False)):
        values = pd.to_numeric(sample_df[metric], errors="coerce").dropna().to_numpy(dtype=float)
        if values.size == 0:
            continue

        counts, _ = np.histogram(values, bins=edges)
        fractions = counts / values.size
        customdata = np.column_stack((edges[:-1], edges[1:], counts, fractions))
        color = palette[index % len(palette)]
        fig.add_trace(
            go.Bar(
                x=centers,
                y=fractions,
                width=widths,
                name=str(sample),
                marker={"color": color, "line": {"width": 0}},
                opacity=0.4,
                customdata=customdata,
                hovertemplate=(
                    "sample=%{fullData.name}<br>"
                    f"{metric}=%{{customdata[0]:.3g}}-%{{customdata[1]:.3g}}<br>"
                    "count=%{customdata[2]:,}<br>"
                    "fraction=%{customdata[3]:.3f}<extra></extra>"
                ),
            )
        )

    fig.update_layout(
        title=metric,
        xaxis_title=metric,
        yaxis_title="Fraction of observations",
        barmode="overlay",
        bargap=0,
        hovermode="x unified",
        legend_title_text="Sample",
    )

    return fig


def versions_yaml(process_name, list_of_libs=None):
    """
    Generate YAML formatted string with versions of relevant libraries.

    Parameters
    ----------
    process_name : str
        Process name to use as key in the versions dictionary.
    list_of_libs : list of str, optional
        List of specific library names to include in the versions dictionary.

    Returns
    -------
    str
        YAML formatted string containing library versions and Python version.
    """

    versions = {}
    versions[process_name] = {}

    versions[process_name]['python'] = f"{sys.version_info.major}" \
        f".{sys.version_info.minor}.{sys.version_info.micro}"

    for lib in list_of_libs:
        try:
            version = importlib.metadata.version(lib)
        except importlib.metadata.PackageNotFoundError:
            try:
                module = importlib.import_module(lib)
                version = getattr(module, '__version__', 'unknown')
            except (ImportError, AttributeError):
                version = None
        if version is not None:
            versions[process_name][lib] = version

    return yaml.dump(versions)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Spatial quality control")
    parser.add_argument("--qc-metrics", nargs="+",
                        help="List of QC metrics summary CSV files to collect")
    parser.add_argument("--annotated-obs", nargs="+",
                        help="List of annotated obs summary CSV files to collect")
    parser.add_argument("--versions-dict", type=str,
                        help="Return dictionary of versions used by the module and exit")

    args = parser.parse_args()

    if args.versions_dict:

        libs = [
            "numpy",
            "pandas",
            "plotly",
        ]
        print(versions_yaml(args.versions_dict, libs))

    else:

        qc_metrics_csv_files = args.qc_metrics
        collect_qc_metrics_files(qc_metrics_csv_files, "qc_metrics.csv")

        annotated_obs_csv_files = args.annotated_obs
        plot_qc_distributions(annotated_obs_csv_files, "qc_distributions.html")
