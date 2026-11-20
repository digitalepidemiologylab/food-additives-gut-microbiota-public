import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor
from patsy import dmatrices
from statsmodels.stats.diagnostic import het_breuschpagan

from additives_microbiota.plotting import process_string


def fit_regression_model(
    data,
    independent_vars,
    dependent_var,
    demographic_vars=None,
    robust_regression=False,
    do_robust_se=False,
):
    """
    Fits a model using independent variables and optional demographic variables.
    If robust_regression is True, uses robust regression (RLM) instead of OLS.
    """
    # Create formula for model
    formula = f"{dependent_var} ~ " + " + ".join(independent_vars)
    if demographic_vars:
        formula += " + " + " + ".join(demographic_vars)

    if robust_regression:
        # Use Robust Linear Model with Huber's T norm
        model = smf.rlm(formula=formula, data=data, M=sm.robust.norms.HuberT()).fit()
    else:
        # Use standard OLS
        if do_robust_se:
            model = smf.ols(formula=formula, data=data).fit(cov_type="HC3")
        else:
            model = smf.ols(formula=formula, data=data).fit()

    return model


def calculate_vif(model, data):
    """Calculate VIF for each predictor in a fitted OLS model, handling categorical variables."""

    # Recreate the design matrix used in the model
    formula = model.model.formula
    y, X = dmatrices(formula, data=data, return_type="dataframe")

    # Calculate VIF (excluding intercept)
    vif_data = pd.DataFrame()
    vif_data["Variable"] = X.columns[1:]  # Skip intercept
    vif_data["VIF"] = [
        variance_inflation_factor(X.values, i) for i in range(1, X.shape[1])
    ]

    return vif_data.sort_values("VIF", ascending=False)


def modeling_stat_tests(
    data, independent_vars, dependent_var, demographic_vars=None, bootstrap_var=None
):
    """
    Perform statistical testing for an additive effect on a dependent variable.
    Returns a dictionary of statistical results:
    1. Non-robust model 2. Robust model 3. Heteroskedasticity test 4. VIF results 5. Bootstrap CI
    """
    stat_res = {}
    model = fit_regression_model(
        data,
        independent_vars,
        dependent_var,
        demographic_vars=demographic_vars,
        robust_regression=False,
        do_robust_se=False,
    )
    stat_res["non_robust_model"] = model.summary()

    # robust_model = fit_regression_model(
    #     data,
    #     independent_vars,
    #     dependent_var,
    #     demographic_vars=demographic_vars,
    #     robust_regression=True,
    # )
    robust_model = fit_regression_model(
        data,
        independent_vars,
        dependent_var,
        demographic_vars=demographic_vars,
        robust_regression=False,
        do_robust_se=True,
    )
    stat_res["robust_model"] = robust_model.summary()

    # bp = het_breuschpagan(robust_model.resid, robust_model.model.exog)
    bp = het_breuschpagan(model.resid, model.model.exog)
    bp_labels = ["Lagrange multiplier stat", "LM p-value", "F-stat", "F p-value"]
    stat_res["heteroskedasticity_test"] = dict(zip(bp_labels, bp))

    # vif_results = calculate_vif(robust_model, data)
    vif_results = calculate_vif(model, data)
    stat_res["vif_results"] = vif_results

    #################################################
    # Bootstrap
    if bootstrap_var:
        n_boot = 1000
        boot_coefs = []
        for _ in range(n_boot):
            bs = data.sample(frac=1, replace=True)
            m = fit_regression_model(
                bs, [bootstrap_var], dependent_var, demographic_vars=demographic_vars
            )
            boot_coefs.append(m.params[bootstrap_var])
        lower, upper = np.percentile(boot_coefs, [2.5, 97.5])
        stat_res["bootstrap_ci"] = (lower, upper)

    #################################################
    # Remove Influential Points
    stat_res["removed_influential_rows_robust_model"], clean_model = (
        remove_influential_points(
            model,
            data,
            independent_vars,
            demographic_vars,
            stat_res,
            dependent_var=dependent_var,
        )
    )

    return stat_res, model, robust_model, clean_model


def remove_influential_points(
    model, meta, mixtures, demographic_vars, stat_res, dependent_var="shannon_entropy"
):
    """
    Identifies and removes influential points using Cook's distance,
    then refits the model on the cleaned data.
    """
    influence = model.get_influence()
    cooks = influence.cooks_distance[0]
    print("len(cooks):", len(cooks))

    idx = model.model.data.row_labels
    cooks_series = pd.Series(cooks, index=idx)

    # print("Top 5 Cook's distances:\n", cooks_series.nlargest(5))

    # 1️⃣ Compute threshold and identify all high‑influence rows
    threshold = 4 / len(cooks)
    high_infl = cooks_series[cooks_series > threshold].index
    print("len(high_infl):", len(high_infl))

    # 2️⃣ Drop them from your DataFrame before fitting
    clean_meta = meta.drop(index=high_infl)

    # 3️⃣ Refit the model on the cleaned data
    clean_model = fit_regression_model(
        clean_meta,
        mixtures,
        dependent_var=dependent_var,
        demographic_vars=demographic_vars,
        robust_regression=False,
        do_robust_se=True,
    )
    return clean_model.summary(), clean_model


def forest_plot(
    results,
    features=None,
    title="Forest plot of cleaned model coefficients",
    figsize=(4, 3),
    r2_annot=True,
    r2_annot_loc=(0.98, 0.98),
    pretty_var_names=None,
    highlight_vars=None,
    highlight_color="crimson",
    base_color="purple",
):
    ci = results.conf_int(alpha=0.05)
    df = (
        pd.DataFrame({"coef": results.params, "lower": ci[0], "upper": ci[1]})
        .drop("Intercept", errors="ignore")
        .loc[lambda x: ~x.index.str.contains("batch")]
    )

    if features is not None:
        df = df.loc[[f for f in features if f in df.index]]

    if pretty_var_names is None:
        pretty_var_names = {}
    if highlight_vars is None:
        highlight_vars = []

    err_lower = df["coef"] - df["lower"]
    err_upper = df["upper"] - df["coef"]

    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(len(df))

    # split into highlighted vs not
    is_hi = df.index.isin(highlight_vars)

    # non-highlighted
    ax.errorbar(
        df.loc[~is_hi, "coef"],
        y[~is_hi],
        xerr=[err_lower[~is_hi], err_upper[~is_hi]],
        fmt="o",
        capsize=4,
        markersize=5,
        color=base_color,
    )

    # highlighted
    ax.errorbar(
        df.loc[is_hi, "coef"],
        y[is_hi],
        xerr=[err_lower[is_hi], err_upper[is_hi]],
        fmt="o",
        capsize=4,
        markersize=6,
        color=highlight_color,
    )

    ax.axvline(0, linestyle="--", color="grey", alpha=0.5)
    ax.set_yticks(y)

    processed_labels = [
        pretty_var_names.get(idx, process_string(idx)) for idx in df.index
    ]
    ax.set_yticklabels(processed_labels)

    if r2_annot:
        r2 = results.rsquared
        ax.text(
            r2_annot_loc[0],
            r2_annot_loc[1],
            f"$R^2 = {r2:.2f}$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=10,
        )

    ax.set_xlabel("Coefficient estimate (Shannon diversity)", fontsize=10.5)
    ax.set_title(title)
    fig.tight_layout()
    return fig, ax
