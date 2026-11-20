from scipy.stats import mannwhitneyu, ttest_ind, kruskal, f_oneway
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from itertools import combinations
import re
from statsmodels.stats.outliers_influence import variance_inflation_factor

from additives_microbiota.plotting import process_string_pretty_vars

# from additives_microbiota.plotting import process_string


def significance_test(data, group_col, value_col, alpha=0.05, non_parametric=False):
    groups = data[group_col].unique()
    if len(groups) == 2:
        group1 = data[data[group_col] == groups[0]][value_col]
        group2 = data[data[group_col] == groups[1]][value_col]
        if non_parametric:
            stat, p_val = mannwhitneyu(group1, group2)
            test_name = "Mann-Whitney U"
        else:
            stat, p_val = ttest_ind(group1, group2)
            test_name = "t-test"
        return p_val < alpha, p_val, test_name, stat
    elif len(groups) > 2:
        group_values = [data[data[group_col] == grp][value_col] for grp in groups]
        if non_parametric:
            stat, p_val = kruskal(*group_values)
            test_name = "Kruskal-Wallis"
        else:
            stat, p_val = f_oneway(*group_values)
            test_name = "ANOVA"
        return p_val < alpha, p_val, test_name, stat
    return False, 1.0, "No test", 0


def create_quartiles(df, column_name, ncut=4):
    """Add quartile categories (1-4) for a specified column in a DataFrame."""
    df_copy = df.copy()
    df_copy[column_name + "_quartile"] = pd.qcut(
        df_copy[column_name], q=ncut, labels=[f"Q{i}" for i in range(1, ncut + 1)]
    )
    return df_copy


def add_stat_annotation(ax, data, x, y, pairs, test="Mann-Whitney", text_format="star"):
    """Add statistical annotation bars to the plot."""
    annotator = Annotator(ax, pairs, data=data, x=x, y=y)
    annotator.configure(
        test=test,
        text_format=text_format,
        line_width=0.5,  # Thinner lines
        line_height=0.02,
    )  # Slightly shorter bars
    annotator.apply_and_annotate()


def _natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(s))]


# def test_quartile_differences(
#     df,
#     column_name,
#     value_col,
#     extremes=True,
#     alpha=0.05,
#     non_parametric=True,
#     ax=None,
#     palette="Spectral",
#     make_quartile=True,
#     ncut=4,
#     showfliers=True,
# ):
#     """
#     Test differences between quartiles of a specified column and create a boxplot with significance bars.
#     Ensures x-order, pairs, and labels are consistent.
#     """
#     df = df.copy()
#     if make_quartile:
#         df = create_quartiles(df, column_name, ncut=ncut)
#         quartile_col = column_name + "_quartile"
#     else:
#         quartile_col = column_name

#     if extremes:
#         df = df[df[quartile_col].isin(["Q1", "Q4"])]
#         if isinstance(df[quartile_col].dtype, pd.CategoricalDtype):
#             df[quartile_col] = df[quartile_col].cat.remove_unused_categories()

#     if ax is None:
#         fig, ax = plt.subplots(1, 1, figsize=(3, 3))

#     # Decide canonical order
#     if isinstance(df[quartile_col].dtype, pd.CategoricalDtype):
#         groups_order = list(df[quartile_col].cat.categories)
#     else:
#         groups_order = sorted(df[quartile_col].unique(), key=_natural_key)

#     # Enforce that order on the data
#     df[quartile_col] = pd.Categorical(
#         df[quartile_col], categories=groups_order, ordered=True
#     )

#     # Boxplot with fixed order and deprecation-safe hue usage
#     sns.boxplot(
#         data=df,
#         x=quartile_col,
#         y=value_col,
#         hue=quartile_col,
#         order=groups_order,
#         hue_order=groups_order,
#         dodge=False,
#         palette=palette,
#         ax=ax,
#         linewidth=0.5,
#         fliersize=1,
#         width=0.5,
#         legend=False,
#         showfliers=showfliers,
#     )

#     # Pairs from the same order
#     if extremes:
#         pairs = [("Q1", "Q4")]
#     else:
#         pairs = list(combinations(groups_order, 2))

#     ax.set_xlabel(process_string(column_name, newLineSep=4))
#     ax.set_ylabel(process_string(value_col, newLineSep=4))

#     test = "Mann-Whitney" if non_parametric else "t-test"
#     add_stat_annotation(ax, df, x=quartile_col, y=value_col, pairs=pairs, test=test)

#     is_sig, p_val, test_name, stat = significance_test(
#         df, quartile_col, value_col, alpha=alpha, non_parametric=non_parametric
#     )

#     if is_sig:
#         print(
#             f"{column_name} -> {value_col}\tSignificant differences exist.",
#             f"{test_name} statistic: {stat:.3f} p-value: {p_val:.3f}",
#         )
#     else:
#         print(f"{column_name} -> {value_col}\tNo significant differences.")

#     ylims = ax.get_ylim()
#     yrange = ylims[1] - ylims[0]
#     text_y = ylims[1] + yrange * 0.05

#     test_abbrev = {
#         "Kruskal-Wallis": "KW",
#         "Mann-Whitney U": "MW",
#         "ANOVA": "ANOVA",
#         "t-test": "t-test",
#     }
#     abbreviated_name = test_abbrev.get(test_name, test_name)

#     if test_name == "Kruskal-Wallis":
#         test_info = f"${abbreviated_name}:$ $H={stat:.2f}$, $p={p_val:.3f}$"
#     elif test_name == "Mann-Whitney U":
#         test_info = f"${abbreviated_name}:$ $U={stat:.2f}$, $p={p_val:.3f}$"
#     else:
#         test_info = f"${abbreviated_name}:$ $t={stat:.2f}$, $p={p_val:.3f}$"

#     ax.text(
#         0,
#         text_y,
#         test_info,
#         horizontalalignment="left",
#         verticalalignment="bottom",
#         fontsize=8,
#     )

#     ax.set_ylim(ylims[0], ylims[1] + yrange * 0.15)
#     return ax


def test_group_differences(
    df,
    group_col,
    value_col,
    specific_groups=None,
    alpha=0.05,
    non_parametric=True,
    ax=None,
    palette="Spectral",
    compare_all_pairs=True,
    showfliers=True,
    add_counts=True,
    xlabel=None,
    ylabel=None,
    verbose=True,
    box_width=0.5,
    sharey=False,
    log_transform=False,
):
    df = df.copy()

    if specific_groups is not None:
        df = df[df[group_col].isin(specific_groups)]
        if isinstance(df[group_col].dtype, pd.CategoricalDtype):
            df[group_col] = df[group_col].cat.remove_unused_categories()

    # Decide which value column to use (original vs log1p transform)
    value_col_used = value_col
    used_log_transform = False
    if log_transform:
        value_col_used = value_col + "_log"
        df[value_col_used] = np.log1p(df[value_col])
        used_log_transform = True

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(3, 3), sharey=sharey)

    if isinstance(df[group_col].dtype, pd.CategoricalDtype):
        groups_order = list(df[group_col].cat.categories)
    else:
        groups_order = sorted(df[group_col].unique(), key=_natural_key)

    df[group_col] = pd.Categorical(df[group_col], categories=groups_order, ordered=True)

    sns.boxplot(
        data=df,
        x=group_col,
        y=value_col_used,
        hue=group_col,
        order=groups_order,
        hue_order=groups_order,
        dodge=False,
        palette=palette,
        ax=ax,
        linewidth=0.5,
        showfliers=showfliers,
        fliersize=1,
        width=box_width,
        legend=False,
    )

    if compare_all_pairs:
        pairs = list(combinations(groups_order, 2))
    else:
        pairs = [
            (groups_order[i], groups_order[i + 1]) for i in range(len(groups_order) - 1)
        ]

    if add_counts:
        cnts = df[group_col].value_counts()
        labels = [f"{grp}\n(n={cnts.get(grp, 0)})" for grp in groups_order]
        ax.set_xticks(range(len(groups_order)))
        ax.set_xticklabels(labels)

    if xlabel:
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel(process_string_pretty_vars(group_col, newLineSep=4))

    if ylabel:
        if used_log_transform:
            ax.set_ylabel(f"{ylabel}\n(log1p transformed)")
        else:
            ax.set_ylabel(ylabel)
    else:
        base_label = process_string_pretty_vars(value_col, newLineSep=4)
        if used_log_transform:
            ax.set_ylabel(f"{base_label}\n(log1p transformed)")
        else:
            ax.set_ylabel(base_label)

    test = "Mann-Whitney" if non_parametric else "t-test"

    if len(groups_order) > 1:
        add_stat_annotation(
            ax, df, x=group_col, y=value_col_used, pairs=pairs, test=test
        )

    is_sig, p_val, test_name, stat = significance_test(
        df, group_col, value_col_used, alpha=alpha, non_parametric=non_parametric
    )

    if verbose:
        if is_sig:
            print(
                f"{group_col} -> {value_col_used}\tSignificant differences exist.",
                f"{test_name} statistic: {stat:.3f} p-value: {p_val:.3f}",
            )
        else:
            print(f"{group_col} -> {value_col_used}\tNo significant differences.")

    ylims = ax.get_ylim()
    yrange = ylims[1] - ylims[0]
    text_y = ylims[1] + yrange * 0.05

    test_abbrev = {
        "Kruskal-Wallis": "KW",
        "Mann-Whitney U": "MW",
        "ANOVA": "ANOVA",
        "t-test": "t-test",
    }
    abbreviated_name = test_abbrev.get(test_name, test_name)

    if test_name == "Kruskal-Wallis":
        test_info = f"${abbreviated_name}:$ $H={stat:.2f}$, $p={p_val:.3f}$"
    elif test_name == "Mann-Whitney U":
        test_info = f"${abbreviated_name}:$ $U={stat:.2f}$, $p={p_val:.3f}$"
    else:
        test_info = f"${abbreviated_name}:$ $t={stat:.2f}$, $p={p_val:.3f}$"

    ax.text(
        0,
        text_y,
        test_info,
        horizontalalignment="left",
        verticalalignment="bottom",
        fontsize=8,
    )

    ax.set_ylim(ylims[0], ylims[1] + yrange * 0.15)
    return ax


def test_quartile_differences(
    df,
    column_name,
    value_col,
    extremes=True,
    alpha=0.05,
    non_parametric=True,
    ax=None,
    palette="Spectral",
    make_quartile=True,
    ncut=4,
    showfliers=True,
    add_counts=True,
    xlabel=None,
    ylabel=None,
    sharey=False,
    log_transform=False,
):
    """
    Test differences between quartiles of a specified column and create a boxplot with significance bars.
    Wrapper around test_group_differences for backward compatibility.
    """
    df = df.copy()
    if make_quartile:
        df = create_quartiles(df, column_name, ncut=ncut)
        quartile_col = column_name + "_quartile"
    else:
        quartile_col = column_name

    specific_groups = ["Q1", "Q4"] if extremes and make_quartile else None

    return test_group_differences(
        df=df,
        group_col=quartile_col,
        value_col=value_col,
        specific_groups=specific_groups,
        alpha=alpha,
        non_parametric=non_parametric,
        ax=ax,
        palette=palette,
        compare_all_pairs=not extremes,
        showfliers=showfliers,
        add_counts=add_counts,
        xlabel=xlabel,
        ylabel=ylabel,
        sharey=sharey,
        log_transform=log_transform,
    )


def calculate_vif(model):
    X_exog = model.model.exog
    result = pd.DataFrame(
        {
            "variable": model.model.exog_names,
            "VIF": [
                variance_inflation_factor(X_exog, i) for i in range(X_exog.shape[1])
            ],
        }
    )

    res_no_intercept = result[result["variable"] != "Intercept"]
    if res_no_intercept["VIF"].max() > 5:
        print(res_no_intercept[res_no_intercept["VIF"] > 5])
    return res_no_intercept
