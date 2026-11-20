import pandas as pd
import numpy as np

# from os import listdir
import pickle
import itertools
from typing import List, Any

# from config.variables import (
#     ALL_VARS,
#     get_cv_vars,
# )

# from utils.data_processing import combine_feature_importances, process_string

from utils.model_utils import (
    model_feature_classifier_performance,
    model_feature_regressor_performance,
)
from tqdm import tqdm

TOP_N_FEATS = 40
NUM_ITER_CLR = 50  # Modify this to reduce the number of iterations!
NUM_ITER_REG = 50  # Modify this to reduce the number of iterations!

# --------------- Variables -----------------------------
slurm_dir = "./slurm_outputs"

additive_vars = [
    "preservatives_antioxidants",
    "thickeners_stabilizers_emulsifiers_gelling_agents",
    "colorants",
    "ph_regulators_anti_caking_leavening_agents",
    "other_additives",
    "sweeteners_glazing_agents",
    "flavor_enhancers",
    "total_additive_intake",
]

mixtures = ["Mixture1", "Mixture2", "Mixture3"]

# --------------- Metadata -----------------------------
meta = pd.read_csv("data/meta_zscores.csv", index_col=0)

# --------------- Counts data -----------------------------
counts_data_species = pd.read_csv(
    "data/counts_data_with_taxonomy.csv",
    index_col=0,
)

counts_data_asv = pd.read_csv(
    "data/table_rarefied.tsv",
    sep="\t",
    index_col=0,
    skiprows=1,
).T

# --------------- Birdman ASVs -----------------------------


def pick_birdman_asv(counts_data, vars):
    """
    Picks the differentially abundant ASVs from the birdman analysis from the vars list
    """
    significant_feats = set()
    for var in vars:
        diff_abun_path = f"{slurm_dir}/{var}/results/differential_taxa.tsv"
        diff_abun = pd.read_csv(diff_abun_path, sep="\t", index_col=0)
        significant_feats = significant_feats.union(set(diff_abun["Feature"].values))
    return counts_data[list(significant_feats)]


# --------------- Assign low/high for Classifier -----------------------------


def assign_low_high(meta, var, low_desc="25%", high_desc="75%"):
    """
    Assigns 'low', 'high', or 'mid' labels to rows based on threshold values.
    """
    meta = meta.copy()

    # Get descriptive statistics
    desc = meta[var].describe()
    low_threshold = desc[low_desc]
    high_threshold = desc[high_desc]

    # Initialize the variable column if it doesn't exist
    if "variable" not in meta.columns:
        meta["variable"] = "mid"

    # Assign labels based on thresholds
    meta.loc[meta[var] <= low_threshold, "variable"] = "low"
    meta.loc[meta[var] >= high_threshold, "variable"] = "high"

    return meta


# --------------- Run classification analysis -----------------------------
clr_res_low_high_auroc = {}
clr_res_low_high_auprc = {}
clr_res_quintiles_auroc = {}
clr_res_quintiles_auprc = {}
clr_res_quintiles_featImp = {}
clr_res_low_high_featImp = {}

var_list_clr = (
    [i for i in meta.columns if "HEI_minus" in i]
    + [i for i in meta.columns if "fruits_minus" in i]
    + [i for i in meta.columns if "vegetables_minus" in i]
    + ["HEI", "fruits_fg_eaten", "vegetables_fg_eaten", "fast_food_fg_eaten"]
)

# Classification with low/high
for var in var_list_clr:
    print(var)
    m = assign_low_high(meta, var)
    # c = pick_birdman_asv(counts_data_asv, var_list_clr)
    auroc_iters, auprc_iters, featImp_iters, fpr_iters, tpr_iters = (
        model_feature_classifier_performance(
            meta=m,
            mb_data=counts_data_species,
            target_col="variable",
            levels=["low", "high"],
            num_iter=NUM_ITER_CLR,
            make_quartiles=False,
        )
    )
    top_features = (
        pd.concat(featImp_iters)
        .groupby("Feature")
        .mean()
        .sort_values(by="Importance", ascending=False)[:TOP_N_FEATS]
    )
    clr_res_low_high_auroc[var] = auroc_iters
    clr_res_low_high_auprc[var] = auprc_iters
    clr_res_low_high_featImp[var] = featImp_iters
    clr_res_low_high_featImp[var] = top_features


pd.DataFrame(clr_res_low_high_auroc).to_csv("outputs/clr_res_low_high_auroc.csv")
pd.DataFrame(clr_res_low_high_auprc).to_csv("outputs/clr_res_low_high_auprc.csv")

with open("outputs/clr_res_low_high_featImp.pkl", "wb") as f:
    pickle.dump(clr_res_low_high_featImp, f)

# Classification with quintiles
for var in var_list_clr:
    # c = pick_birdman_asv(counts_data_asv, var_list_clr)
    try:
        auroc_iters, auprc_iters, featImp_iters, fpr_iters, tpr_iters = (
            model_feature_classifier_performance(
                meta=meta,
                mb_data=counts_data_species,
                target_col=var,
                levels=["Q1", "Q5"],
                num_iter=NUM_ITER_CLR,
                make_quartiles=True,
                quartile_num=5,
            )
        )

        top_features = (
            pd.concat(featImp_iters)
            .groupby("Feature")
            .mean()
            .sort_values(by="Importance", ascending=False)[:TOP_N_FEATS]
        )

        clr_res_quintiles_auroc[var] = auroc_iters
        clr_res_quintiles_auprc[var] = auprc_iters
        clr_res_quintiles_featImp[var] = featImp_iters
        clr_res_quintiles_featImp[var] = top_features

    except Exception as e:
        print(f"Classification with quintiles: Error for {var}: {e}")

# Saving results
pd.DataFrame(clr_res_quintiles_auroc).to_csv("outputs/clr_res_quintiles_auroc.csv")
pd.DataFrame(clr_res_quintiles_auprc).to_csv("outputs/clr_res_quintiles_auprc.csv")

with open("outputs/clr_res_quintiles_featImp.pkl", "wb") as f:
    pickle.dump(clr_res_quintiles_featImp, f)


# --------------- Run regression analysis -----------------------------
def run_regression_analysis(main_var, meta, counts_data, var_list):
    """
    Run regression analysis for a main variable and its interactions with additives.

    Parameters:
    - main_var: String, the main variable to analyze
    - meta: DataFrame, metadata
    - counts_data: DataFrame, count data
    - var_list: List, variables to analyze
    - NUM_ITER_REG: Integer, number of iterations for regression

    Returns:
    - Dictionary containing correlation results
    """
    spr_corrs_results = {}
    if main_var == "HEI":
        var = "HEI"
    else:
        var = f"{main_var}_fg_eaten"
    spr_corrs = model_feature_regressor_performance(
        meta,
        counts_data,
        var,
        sample_id_col="sample-id",
        objective="reg:squarederror",
        num_iter=NUM_ITER_REG,
        verbose=False,
        scale_target=False,
    )
    spr_corrs_results[var] = spr_corrs

    print(f"{main_var}:", round(np.mean(spr_corrs), 2), round(np.std(spr_corrs), 2))
    #############

    for additive in tqdm(var_list):
        if additive == "fast_food":
            spr_corrs_fast_food = model_feature_regressor_performance(
                meta,
                counts_data,
                "fast_food_fg_eaten",
                sample_id_col="sample-id",
                objective="reg:squarederror",
                num_iter=NUM_ITER_REG,
                verbose=False,
                scale_target=False,
            )
            spr_corrs_results["fast_food_fg_eaten"] = spr_corrs_fast_food

        spr_corrs = model_feature_regressor_performance(
            meta,
            counts_data,
            f"{main_var}_minus_{additive}",
            sample_id_col="sample-id",
            objective="reg:squarederror",
            num_iter=NUM_ITER_REG,
            verbose=False,
            scale_target=False,
        )
        spr_corrs_results[f"{main_var}_minus_{additive}"] = spr_corrs

        print(
            f"({main_var}) - ({additive}):",
            round(np.mean(spr_corrs), 2),
            round(np.std(spr_corrs), 2),
        )

    return spr_corrs_results


# Run analysis for fruits, HEI, and vegetables
main_vars = ["HEI", "fruits", "vegetables"]
var_list_reg = additive_vars + mixtures + ["fast_food"]

for main_var in main_vars:
    print(f"\nRunning analysis for {main_var}...")

    # Pick the Birdman ASVs that are differentially abundant for the main variable
    # counts_data_specific = pick_birdman_asv(counts_data_asv, var_list_clr)

    spr_corrs_results = run_regression_analysis(
        main_var, meta, counts_data_species, var_list_reg
    )

    spr_corrs_df = pd.DataFrame(
        {key: pd.Series(values) for key, values in spr_corrs_results.items()}
    )
    spr_corrs_df.to_csv(f"outputs/spr_corrs_main_var_{main_var}.csv")
    print(f"Results saved to outputs/spr_corrs_main_var_{main_var}.csv")

# --------------- Run Classification with Additives Combinations -----------------------------


def get_all_combinations(
    items: List[Any], min_size: int = 2, max_size: int = None
) -> List[tuple]:
    """
    Generate all combinations of items from min_size to max_size elements.

    Parameters:
    - items: List of items to generate combinations from
    - min_size: Minimum number of items in each combination (default: 2)
    - max_size: Maximum number of items in each combination (default: length of items)

    Returns:
    - List of tuples, where each tuple is a combination
    """
    if max_size is None:
        max_size = len(items)

    # Ensure min_size and max_size are valid
    min_size = max(1, min_size)
    max_size = min(len(items), max_size)

    # Generate all combinations from min_size to max_size
    all_combinations = []
    for r in range(min_size, max_size + 1):
        combinations_of_size_r = list(itertools.combinations(items, r))
        all_combinations.extend(combinations_of_size_r)

    return all_combinations


few_vars = [
    "thickeners_stabilizers_emulsifiers_gelling_agents",
    "preservatives_antioxidants",
    "ph_regulators_anti_caking_leavening_agents",
    "sweeteners_glazing_agents",
    "colorants",
]

combinations = get_all_combinations(few_vars)
print(f"Total number of combinations: {len(combinations)}")

ml_res_clr_low_high_combo = {}
ml_res_clr_quintiles_combo = {}

for var_combo in tqdm(combinations):
    m_copy = meta.copy()
    m_copy["combo_additives"] = m_copy[list(var_combo)].sum(axis=1)
    print(var_combo)
    m_copy = assign_low_high(m_copy, "combo_additives")
    try:
        auroc_iters, auprc_iters, featImp_iters, fpr_iters, tpr_iters = (
            model_feature_classifier_performance(
                m_copy,
                counts_data_species,
                "variable",
                levels=["low", "high"],
                num_iter=NUM_ITER_CLR,
                make_quartiles=False,
            )
        )
        ml_res_clr_low_high_combo[var_combo] = np.mean(auroc_iters)
    except Exception as e:
        print(e)
        continue

pd.DataFrame([ml_res_clr_low_high_combo], index=["mean_auroc"]).to_csv(
    "outputs/ml_res_clr_low_high_combo.csv"
)

for var_combo in combinations:
    print(var_combo)
    m_copy = meta.copy()
    m_copy["combo_additives"] = m_copy[list(var_combo)].sum(axis=1)
    m_copy = assign_low_high(m_copy, "combo_additives")
    try:
        auroc_iters, auprc_iters, featImp_iters, fpr_iters, tpr_iters = (
            model_feature_classifier_performance(
                m_copy,
                counts_data_species,
                "variable",
                levels=["Q1", "Q5"],
                num_iter=NUM_ITER_CLR,
                make_quartiles=True,
                quartile_num=5,
            )
        )
        ml_res_clr_quintiles_combo[var_combo] = np.mean(auroc_iters)
    except Exception as e:
        print(e)
        continue

pd.DataFrame([ml_res_clr_quintiles_combo], index=["mean_auroc"]).to_csv(
    "outputs/ml_res_clr_quintiles_combo.csv"
)

# --------------- Run Classification with Mixtures Combinations -----------------------------

combinations_mixtures = get_all_combinations(mixtures, min_size=1)

ml_res_clr_mixtures_low_high_combo = {}
ml_res_clr_mixtures_quintiles_combo = {}

for combo in combinations_mixtures:
    print(combo)
    m_copy = meta.copy()
    m_copy["combo_mixtures"] = m_copy[list(combo)].sum(axis=1)
    m_copy = assign_low_high(m_copy, "combo_mixtures")
    try:
        auroc_iters, auprc_iters, featImp_iters, fpr_iters, tpr_iters = (
            model_feature_classifier_performance(
                m_copy,
                counts_data_species,
                "variable",
                levels=["low", "high"],
                num_iter=NUM_ITER_CLR,
                make_quartiles=False,
            )
        )
        ml_res_clr_mixtures_low_high_combo[combo] = np.mean(auroc_iters)
    except Exception as e:
        print(e)
        continue

pd.DataFrame([ml_res_clr_mixtures_low_high_combo], index=["mean_auroc"]).to_csv(
    "outputs/ml_res_clr_mixtures_low_high_combo.csv"
)

for combo in combinations_mixtures:
    print(combo)
    m_copy = meta.copy()
    m_copy["combo_mixtures"] = m_copy[list(combo)].sum(axis=1)
    m_copy = assign_low_high(m_copy, "combo_mixtures")
    try:
        auroc_iters, auprc_iters, featImp_iters, fpr_iters, tpr_iters = (
            model_feature_classifier_performance(
                m_copy,
                counts_data_species,
                "variable",
                levels=["Q1", "Q5"],
                num_iter=NUM_ITER_CLR,
                make_quartiles=True,
                quartile_num=5,
            )
        )
        ml_res_clr_mixtures_quintiles_combo[combo] = np.mean(auroc_iters)
    except Exception as e:
        print(e)
        continue

pd.DataFrame([ml_res_clr_mixtures_quintiles_combo], index=["mean_auroc"]).to_csv(
    "outputs/ml_res_clr_mixtures_quintiles_combo.csv"
)
