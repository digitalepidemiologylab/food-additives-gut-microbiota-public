import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import os

from additives_microbiota.plotting import process_string

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

path_asv_taxonomy = "~/WORK/qiime2_2024/qiime_processed_files/taxonomy_rarefied-table_2024_09/taxonomy.tsv"

otuID_taxonomy = pd.read_csv(path_asv_taxonomy, sep="\t", index_col=0)
otuID_taxonomy.head()
taxonomy_dict = dict(zip(otuID_taxonomy.index, otuID_taxonomy["Taxon"]))

prevalence = pd.read_csv(
    "~/WORK/qiime2_2024/analysis/outputs/prevalence_microbes.csv",
).set_index("#OTU ID")
prevalence = pd.Series(prevalence["0"].values, index=prevalence.index)

dict_ASV_names = dict(
    zip(
        prevalence.sort_values(ascending=False).index,
        ["ASV_" + str(i) for i in range(len(prevalence))],
    )
)

# ------------------------------------------------------------------------------
# Define functions
# ------------------------------------------------------------------------------


def transform_counts_data(counts_data, method="clr", pseudocount=None, log_base="e"):
    """
    Transform a feature table.

    Parameters
    ----------
    counts_data : pd.DataFrame
        Samples x Features table (rows are samples after your .T).
    method : str
        'clr' (for read counts compositional analysis) OR
        'log1p' (useful for PICRUSt continuous predicted abundances) OR
        'relative_abundance'.
    pseudocount : float or None
        Only used for 'clr'. If None, uses a small data-driven value.
    log_base : str
        'e' for natural log, '10' for log10. Only used for 'log1p'.

    Returns
    -------
    pd.DataFrame
    """
    if method == "clr":
        from skbio.stats.composition import clr

        df = counts_data.copy()

        if pseudocount is None:
            # small, data-driven pseudocount: half of the smallest positive value
            positive = df.to_numpy()
            positive = positive[positive > 0]
            if positive.size == 0:
                raise ValueError("CLR transform requires at least one positive value.")
            pseudocount = 0.5 * positive.min()

        df = df + pseudocount

        df_mb = pd.DataFrame(
            clr(df), index=counts_data.index, columns=counts_data.columns
        )

    elif method == "log1p":
        if log_base == "e":
            df_mb = np.log1p(counts_data)
        elif log_base == "10":
            df_mb = np.log10(counts_data + 1.0)
        else:
            raise ValueError("Invalid log_base. Choose 'e' or '10'.")

    elif method == "relative_abundance":
        df_mb = counts_data.div(counts_data.sum(axis=1), axis=0)

    else:
        raise ValueError(
            "Invalid method. Choose 'clr', 'log1p', or 'relative_abundance'."
        )

    return df_mb


# def preprocess_picrust_data(df, data_type, method="log1p", min_prev=0.10, min_var=1e-8):
#     """
#     Preprocess pathway data by:
#     - transforming counts data,
#     - filtering pathways by minimum prevalence,
#     - filtering by minimum variance.

#     Args:
#         pathway_data (pd.DataFrame): Pathway abundance matrix (samples x pathways).
#         method (str): Transformation method for counts data (default 'log1p').
#     """
#     data_transformed = transform_counts_data(df, method=method)

#     id_dict = {
#         p: f"{data_type}_" + str(e) for e, p in enumerate(data_transformed.columns)
#     }
#     rev_id_dict = {v: k for k, v in id_dict.items()}

#     data_transformed.columns = data_transformed.columns.map(id_dict)

#     is_nonzero = data_transformed.ne(0)
#     keep_features = is_nonzero.mean(axis=0) >= min_prev
#     data_transformed_filt = data_transformed.loc[:, keep_features].copy()

#     # Drop near-zero variance pathways after filtering
#     data_transformed_filt = data_transformed_filt.loc[
#         :, data_transformed_filt.var(axis=0) > min_var
#     ]

#     return data_transformed_filt, id_dict, rev_id_dict


def preprocess_picrust_data(df, data_type, method="log1p", min_prev=0.10, min_var=1e-8):
    """
    Preprocess pathway data by:
    - transforming counts data
    - filtering pathways by minimum prevalence
    - filtering by minimum variance

    Returns prevalence indexed by ORIGINAL pathway names.
    """
    data_transformed = transform_counts_data(df, method=method)

    id_dict = {
        p: f"{data_type}_" + str(e) for e, p in enumerate(data_transformed.columns)
    }
    rev_id_dict = {v: k for k, v in id_dict.items()}

    data_transformed.columns = data_transformed.columns.map(id_dict)

    is_nonzero = data_transformed.ne(0)
    prevalence = is_nonzero.mean(axis=0)

    keep_features = prevalence >= min_prev
    data_transformed_filt = data_transformed.loc[:, keep_features].copy()
    prevalence_filt = prevalence.loc[keep_features].copy()

    data_transformed_filt = data_transformed_filt.loc[
        :, data_transformed_filt.var(axis=0) > min_var
    ]
    prevalence_filt = prevalence_filt.loc[data_transformed_filt.columns]

    prevalence_filt.index = prevalence_filt.index.map(rev_id_dict)

    return data_transformed_filt, id_dict, rev_id_dict, prevalence_filt


def taxonID_shortName(featureID, level="genus"):
    """
    Shorten the name of a taxa to the genus or species level if possible.
    """
    taxa_name = otuID_taxonomy.loc[featureID].values[0].replace("; ", ";")
    asv_id = dict_ASV_names[featureID]

    regular = shorten_taxa_name(extract_taxa_name(taxa_name, level))
    # regular = regular + " (" + asv_id + ")"
    regular = asv_id + " (" + regular + ")"

    if (
        regular.startswith("g__COE1")
        or regular.startswith("s__UBA")
        or regular.startswith("g__G11")
        or regular.startswith("g__C-53")
    ):
        return shorten_taxa_name(extract_taxa_name(taxa_name, "family")) + ";" + regular
    elif (
        regular.startswith("s__CAG")
        or regular.startswith("s__PeH")
        or regular.startswith("g__SFM")
        or regular.startswith("s__QALR01")
    ):
        return shorten_taxa_name(extract_taxa_name(taxa_name, "order")) + ";" + regular

    return regular


def shorten_taxa_name(taxa_name):
    """
    Shorten the name of a taxa to the genus or species level if possible.
    """
    if taxa_name.startswith("species__"):
        return taxa_name.split("__")[1]
    elif taxa_name.startswith("genus__"):
        return taxa_name.split("__")[1]
    elif taxa_name.startswith("family__"):
        return taxa_name.split("__")[1]
    else:
        return taxa_name


def extract_taxa_name(taxonomy_string, level="species"):
    """
    Extracts the name at a specified taxa level from a taxonomy string.
    If the name at that level is absent, it searches higher levels until a name is found.
    """
    # Define the order of levels to facilitate searching upwards
    levels_order = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    levels_abbr = {
        "domain": "d",
        "phylum": "p",
        "class": "c",
        "order": "o",
        "family": "f",
        "genus": "g",
        "species": "s",
    }

    # Split the taxonomy string into components
    taxa_parts = taxonomy_string.split(";")

    # Dictionary to hold each part of the taxonomy with its level as key
    taxa_dict = {}
    for part in taxa_parts:
        taxa_level, _, taxa_name = part.partition("__")
        if taxa_level:  # Check that taxa_level is not empty
            taxa_dict[taxa_level[0]] = taxa_name

    # Search for the desired level and upwards if necessary
    # level_abbr = levels_abbr.get(level, "")
    for i in range(levels_order.index(level), -1, -1):
        current_level = levels_order[i]
        if (
            levels_abbr[current_level] in taxa_dict
            and taxa_dict[levels_abbr[current_level]]
        ):
            return (
                f"{levels_abbr[current_level]}__{taxa_dict[levels_abbr[current_level]]}"
            )

    return "Name not found"


def _unpack_hdi_and_filter(df, col):
    df = df.copy()
    df[["lower", "upper"]] = df[col].str.split(",", expand=True)
    # remove ( from lower and ) from upper and convert to float
    df.lower = df.lower.str[1:].astype("float")
    df.upper = df.upper.str[:-1].astype("float")

    df["credible"] = np.where((df.lower > 0) | (df.upper < 0), "yes", "no")

    df.upper = df.upper - df[col.replace("hdi", "mean")]
    df.lower = df[col.replace("hdi", "mean")] - df.lower

    return df


def _get_top_n_feats(df, n):
    if df.shape[0] < 2 * n:
        df_for_display = df
    else:
        bottomn = df[:n]
        topn = df[-1 * n :]  # noqa: E203
        df_for_display = pd.concat([bottomn, topn])

    return df_for_display


def _display_top_n_feats(
    df,
    n,
    yvar,
    xvar,
    xlab,
    ylab,
    title,
    ax=None,
    figsize=(6, 10),
    annotate_plot_args=None,
):
    if df.shape[0] < 2 * n:
        df_for_display = df
    else:
        bottomn = df[:n]
        topn = df[-1 * n :]  # noqa: E203
        df_for_display = pd.concat([bottomn, topn])

    df_for_display["point_color"] = df_for_display[xvar].apply(
        lambda x: "darkorange" if x < 0 else "purple"
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot points with colors
    sns.stripplot(
        data=df_for_display,
        y=yvar,
        x=xvar,
        hue="point_color",
        palette={"darkorange": "darkorange", "purple": "purple"},
        legend=False,
        ax=ax,
    )

    # Plot error bars with colors for each point
    for idx, row in df_for_display.iterrows():
        color = "darkorange" if row[xvar] < 0 else "purple"
        ax.errorbar(
            x=row[xvar],
            y=row[yvar],
            xerr=[[row["lower"]], [row["upper"]]],
            ls="none",
            color=color,
            alpha=0.5,
        )

    ax.grid(axis="y", color="grey", linestyle="--", linewidth=0.3)
    ax.axvline(x=0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_ylabel("")
    ax.set_xlabel(xlab, fontsize=12)
    ax.set_title(title)

    if annotate_plot_args is not None:
        ax.text(
            annotate_plot_args["textpos"][0],
            annotate_plot_args["textpos"][1],
            annotate_plot_args["text"],
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
            ha="left",
        )


def birdman_plot_single_var(
    df_inf,
    var,
    action="plot",
    n=50,
    plot_path=None,
    annotate_plot_args={"textpos": (0.3, 0.95), "text": ""},
    ax=None,
):
    sub_df = _unpack_hdi_and_filter(df_inf, var + "_hdi")
    sub_df.rename_axis(index="Feature", inplace=True)
    sub_df = sub_df.sort_values(by=var + "_mean")
    # print("Filtered Shape: " + str(sub_df.loc[sub_df["credible"] == "yes"].shape))

    xlab = "Mean Effect Size for " + process_string(var)  # e.g. var = "host_age[T.34]_"
    ylab = "Features"
    df_for_display = sub_df.reset_index()
    df_for_display = df_for_display.loc[df_for_display.credible == "yes"]

    if action == "plot":
        df_for_display["Feature"] = df_for_display["Feature"].map(taxonID_shortName)
        _display_top_n_feats(
            df_for_display,
            n,
            "Feature",
            var + "_mean",
            xlab,
            ylab,
            process_string(var) + f" (Top {n} and Bottom {n} Features)",
            ax=ax,
            annotate_plot_args=annotate_plot_args,
        )

        if plot_path is not None:
            plt.savefig(plot_path, dpi=300, bbox_inches="tight")
            plt.close()
    elif action == "get_top_n":
        return _get_top_n_feats(df_for_display, n)
    else:
        raise ValueError("Invalid action")


def compute_log_ratios_using_significant_features(
    counts_data, var, significant_feats, meta, var_list
):
    """
    Compute log ratios using significant features.

    Parameters:
    counts_data (pd.DataFrame): The counts data.
    significant_feats (pd.DataFrame): The significant features.
    var (str): The variable to compute log ratios for.
    meta (pd.DataFrame): The metadata.

    Returns:
    log_ratios (pd.DataFrame): A DataFrame with log ratios for each feature and sample.
    corr (float): The pearson correlation between the log ratios and the variable.
    corr_dict (dict): A dictionary with the pearson correlation and p-value between the log ratios and each variable in var_list.
    """

    numerator_asvs = significant_feats[significant_feats[f"{var}_mean"] > 0][
        "Feature"
    ].values
    denominator_asvs = significant_feats[significant_feats[f"{var}_mean"] < 0][
        "Feature"
    ].values
    numerdata = counts_data.T.loc[numerator_asvs].sum(axis=0)
    denomdata = counts_data.T.loc[denominator_asvs].sum(axis=0)

    # Removing samples with 0 from numerdata and denomdata
    numerdata = numerdata[numerdata > 0]
    denomdata = denomdata[denomdata > 0]

    common_samples = numerdata.index.intersection(denomdata.index).intersection(
        meta["sample-id"].values
    )

    if len(common_samples) == 0:
        return None, None, None, None, None

    numerdata = numerdata.loc[common_samples]
    denomdata = denomdata.loc[common_samples]

    # Compute Log Ratios
    log_ratios = np.log10(numerdata / denomdata)
    log_ratios = log_ratios.to_frame().reset_index()
    log_ratios.columns = ["sample-id", "Log Ratio"]
    log_ratios = log_ratios.merge(meta, on="sample-id")

    # find pearson correlation between log_ratios and Var
    corr, pval = pearsonr(log_ratios["Log Ratio"], log_ratios[var])
    print(
        f"Common Samples: {common_samples.shape}",
        f"\t{var}",
        f"\t Pearson Correlation {var}:",
        round(corr, 2),
        f"\t P-value: {pval:.2e}",
    )

    corr_dict = {}
    for v in var_list:
        corr_stat, p_val = pearsonr(log_ratios["Log Ratio"], log_ratios[v])
        corr_dict[v] = {"correlation": corr_stat, "pvalue": p_val}

    return log_ratios, corr, corr_dict, numerator_asvs, denominator_asvs


def fetch_folders(directory):
    """
    Fetches all folders within the specified directory.

    Parameters:
    directory (str): The path to the directory to search for folders.

    Returns:
    List[str]: A list of folder names within the specified directory.
    """
    try:
        # Get the list of all items in the specified directory
        all_items = os.listdir(directory)
        # Filter and return only the directories
        folders = [
            item for item in all_items if os.path.isdir(os.path.join(directory, item))
        ]
        return folders
    except FileNotFoundError:
        print(f"The directory {directory} does not exist.")
        return []
    except PermissionError:
        print(f"Permission denied to access {directory}.")
        return []
