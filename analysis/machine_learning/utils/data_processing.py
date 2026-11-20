import pandas as pd


def create_quartiles(df, column_name, q=4):
    """
    Create quartiles for a numerical column in a pandas DataFrame.

    Parameters:
    - df: pandas.DataFrame containing the data.
    - column_name: String, the name of the column to be binned into quartiles.
    - q: Integer, the number of quantiles. Default is 4 (quartiles).

    Returns:
    - A pandas Series with the quartile assignment for each entry in the column.
    """
    quartiles = pd.qcut(df[column_name], q, labels=[f"Q{i}" for i in range(1, q + 1)])
    return quartiles.astype(str)


def combine_feature_importances(dfs):
    """
    Combine multiple DataFrames containing feature importances.

    Parameters:
    - dfs (list of pd.DataFrame): List of DataFrames to combine.

    Returns:
    - pd.DataFrame: Combined DataFrame with unique 'Feature' column and
                   separate 'importance' columns for each input DataFrame.
    """
    combined_df = None

    for i, df in enumerate(dfs):
        # Rename 'importance' column uniquely
        df_renamed = df.rename(columns={"Importance": f"importance_{i}"})
        df_renamed = df_renamed[:50]

        if combined_df is None:
            combined_df = df_renamed
        else:
            # Merge with the combined DataFrame on 'feature', using outer join
            combined_df = pd.merge(combined_df, df_renamed, on="Feature", how="outer")

    return combined_df.fillna(0)
