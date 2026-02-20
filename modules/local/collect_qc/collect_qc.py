
import importlib
import importlib.metadata
import argparse
import sys
import yaml

import pandas as pd
import plotly.express as px
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

            fig = px.violin(combined, y=metric, color="sample", box=True, title=metric)
            if i == 0:
                f.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
            else:
                f.write(fig.to_html(full_html=False, include_plotlyjs=False))


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
            "pandas",
        ]
        print(versions_yaml(args.versions_dict, libs))

    else:

        qc_metrics_csv_files = args.qc_metrics
        collect_qc_metrics_files(qc_metrics_csv_files, "qc_metrics.csv")

        annotated_obs_csv_files = args.annotated_obs
        plot_qc_distributions(annotated_obs_csv_files, "qc_distributions.html")
