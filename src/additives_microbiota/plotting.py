import pandas as pd
from pathlib import Path

pretty_var_names = {
    "age": "Age (years)",
    "bmi": "BMI (kg/m²)",
    "C(smoking, Treatment(reference='non-smoker'))[T.former]": "Smoking: former",
    "C(smoking, Treatment(reference='non-smoker'))[T.smoker]": "Smoking: current",
    # "language[T.latin]": "Language region (French/italian)",
    "C(language, Treatment(reference='latin'))[T.german]": "Language region (German)",
    # "urbanity[T.semi_urban]": "Semi-urban",
    # "urbanity[T.urban]": "Urban",
    "C(urbanity, Treatment(reference='rural'))[T.semi_urban]": "Semi-urban",
    "C(urbanity, Treatment(reference='rural'))[T.urban]": "Urban",
    "defecate_quantity_per_day": "Defecation frequency",
    "C(gender)[T.male]": "Gender: male",
    "C(antibiotics_treatment)[T.True]": "Antibiotics History",
    "eaten_quantity_in_gram": "Total eaten quantity (g/day)",
    "energy_kcal_eaten": "Energy (kcal/day)",
    "alcohol_eaten": "Alcohol (g/day)",
}

pretty_var_names_simple = {
    "carb_eaten": "Carbohydrate",
    "defecate_quantity_per_day": "Defecation frequency",
    "bmi": "BMI",
    "age_group_2": "Age group",
    "thickeners_stabilizers_emulsifiers_gelling_agents": "TSEGs",
    "sweeteners_glazing_agents": "SG",
    "preservatives_antioxidants": "PA",
    "ph_regulators_anti_caking_leavening_agents": "pH/AC/L",
    "other_additives": "Other additives",
    "flavor_enhancers": "Flavor enhancers",
    "colorants": "Colorants",
    "total_additive_intake": "Total additives",
    #
    "vegetables_minus_colorants": "Vegetables - Colorants",
    "vegetables_minus_preservatives_antioxidants": "Vegetables - PA",
    "vegetables_minus_thickeners_stabilizers_emulsifiers_gelling_agents": "Vegetables - TSEGs",
    "vegetables_minus_ph_regulators_anti_caking_leavening_agents": "Vegetables - pH/AC/L",
    "vegetables_minus_flavor_enhancers": "Vegetables - Flavor enhancers",
    "vegetables_minus_sweeteners_glazing_agents": "Vegetables - SG",
    "vegetables_minus_other_additives": "Vegetables - Other additives",
    "vegetables_minus_total_additive_intake": "Vegetables - Total additives",
    "vegetables_minus_fast_food": "Vegetables - Fast food",
    "vegetables_minus_Mixture1": "Vegetables - Mixture1",
    "vegetables_minus_Mixture2": "Vegetables - Mixture2",
    "vegetables_minus_Mixture3": "Vegetables - Mixture3",
    #
    "fruits_minus_colorants": "Fruits - Colorants",
    "fruits_minus_preservatives_antioxidants": "Fruits - PA",
    "fruits_minus_thickeners_stabilizers_emulsifiers_gelling_agents": "Fruits - TSEGs",
    "fruits_minus_ph_regulators_anti_caking_leavening_agents": "Fruits - pH/AC/L",
    "fruits_minus_flavor_enhancers": "Fruits - Flavor enhancers",
    "fruits_minus_sweeteners_glazing_agents": "Fruits - SG",
    "fruits_minus_other_additives": "Fruits - Other additives",
    "fruits_minus_total_additive_intake": "Fruits - Total additives",
    "fruits_minus_fast_food": "Fruits - Fast food",
    "fruits_minus_Mixture1": "Fruits - Mixture1",
    "fruits_minus_Mixture2": "Fruits - Mixture2",
    "fruits_minus_Mixture3": "Fruits - Mixture3",
    #
    "HEI_minus_colorants": "HEI - Colorants",
    "HEI_minus_preservatives_antioxidants": "HEI - PA",
    "HEI_minus_thickeners_stabilizers_emulsifiers_gelling_agents": "HEI - TSEGs",
    "HEI_minus_ph_regulators_anti_caking_leavening_agents": "HEI - pH/AC/L",
    "HEI_minus_flavor_enhancers": "HEI - Flavor enhancers",
    "HEI_minus_sweeteners_glazing_agents": "HEI - SG",
    "HEI_minus_other_additives": "HEI - Other additives",
    "HEI_minus_total_additive_intake": "HEI - Total additives",
    "HEI_minus_fast_food": "HEI - Fast food",
    "HEI_minus_Mixture1": "HEI - Mixture1",
    "HEI_minus_Mixture2": "HEI - Mixture2",
    "HEI_minus_Mixture3": "HEI - Mixture3",
}


def process_string_pretty_vars(
    s,
    ignore_vars=[
        "PDI_Quintile",
        "hPDI_Quintile",
        "HEI",
        "BMI",
    ],
    newLineSep=7,
):
    if s in ignore_vars:
        return s
    if s in pretty_var_names_simple:
        return pretty_var_names_simple[s]
    if s.startswith("LR_"):
        return "Log Ratio (" + process_string_pretty_vars(s.replace("LR_", "")) + ")"
    s = s.replace("_fg", "")
    # Step 1 & 2: Split the string by "_" and replace "_" with a space, then capitalize the first letter
    parts = s.replace("_eaten", "").split("_")
    # capitalized_parts = [part.capitalize() for part in parts]
    capitalized_parts = [part[0].upper() + part[1:] if part else "" for part in parts]

    # Step 3: If there are more than 3 items, insert "\n" in the middle
    if newLineSep and len(capitalized_parts) > newLineSep and "CV" not in s:
        mid_index = len(capitalized_parts) // 2
        processed_string = (
            " ".join(capitalized_parts[:mid_index])
            + "\n"
            + " ".join(capitalized_parts[mid_index:])
        )
    elif "CV" in s:
        processed_string = (
            capitalized_parts[0] + " (" + " ".join(capitalized_parts[1:]) + ")"
        )
    else:
        processed_string = " ".join(capitalized_parts)

    if "Ph " in processed_string:
        processed_string = processed_string.replace("Ph ", "pH ")

    return processed_string


def process_string(
    s,
    ignore_vars=[
        "PDI_Quintile",
        "hPDI_Quintile",
        "HEI",
        "BMI",
    ],
    newLineSep=7,
):
    if s in ignore_vars:
        return s
    if s == "bmi":
        return "BMI"
    if s == "defecate_quantity_per_day":
        return "Defecation frequency"
    s = s.replace("_fg", "")
    # Step 1 & 2: Split the string by "_" and replace "_" with a space, then capitalize the first letter
    parts = s.replace("_eaten", "").split("_")
    # capitalized_parts = [part.capitalize() for part in parts]
    capitalized_parts = [part[0].upper() + part[1:] if part else "" for part in parts]

    # Step 3: If there are more than 3 items, insert "\n" in the middle
    if newLineSep and len(capitalized_parts) > newLineSep and "CV" not in s:
        mid_index = len(capitalized_parts) // 2
        processed_string = (
            " ".join(capitalized_parts[:mid_index])
            + "\n"
            + " ".join(capitalized_parts[mid_index:])
        )
    elif "CV" in s:
        processed_string = (
            capitalized_parts[0] + " (" + " ".join(capitalized_parts[1:]) + ")"
        )
    else:
        processed_string = " ".join(capitalized_parts)

    if "Ph " in processed_string:
        processed_string = processed_string.replace("Ph ", "pH ")

    return processed_string


ecode_category_colors = {
    "colorants": "gold",
    "preservatives_antioxidants": "#66b3ff",
    "thickeners_stabilizers_emulsifiers_gelling_agents": "#ffcc99",
    "ph_regulators_anti_caking_leavening_agents": "#b3de69",
    "flavor_enhancers": "#ffb3e6",
    "sweeteners_glazing_agents": "#c2c2f0",
    "other_additives": "#c2f0c2",
    "unknown": "#cccccc",
}


def get_ecode_category(code):
    base_path = Path(__file__).resolve().parent
    mapping_path = base_path / "../../data/additive_name_mapping_v3.csv"

    additive_name_mapping = pd.read_csv(mapping_path)

    if code in additive_name_mapping["ecode"].values:
        category = additive_name_mapping.loc[
            additive_name_mapping["ecode"] == code, "category"
        ].values[0]

        if category in ["preservatives", "antioxidants"]:
            return "preservatives_antioxidants"
        elif category in ["sweeteners", "glazing_agents"]:
            return "sweeteners_glazing_agents"
        else:
            return category
    else:
        return "unknown"
