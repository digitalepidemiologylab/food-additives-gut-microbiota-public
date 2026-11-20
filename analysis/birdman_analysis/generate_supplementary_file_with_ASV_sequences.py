import pandas as pd
import os

slurm_outputs_dir = "./slurm_outputs/adjustedBatch"


def read_fasta_to_dataframe(fasta_file):
    """
    Reads a FASTA file and returns a DataFrame with hashes (sequence IDs) and sequences.

    Parameters:
    fasta_file (str): Path to the FASTA file.

    Returns:
    pd.DataFrame: DataFrame with 'Hash' and 'Sequence' columns.
    """
    # Initialize lists to store the hashes and sequences
    hashes = []
    sequences = []

    # Read the FASTA file manually
    with open(fasta_file, "r") as file:
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                    sequence = ""
                hashes.append(line[1:])  # Remove '>' and store the hash
            else:
                sequence += line
        sequences.append(sequence)  # Add the last sequence

    # Create and return a DataFrame
    return pd.DataFrame({"Hash": hashes, "Sequence": sequences})


fasta_file_path = (
    "../../../dietary-microbiota-paper/qiime/exported_rep_seqs/dna-sequences.fasta"
)
rdf = read_fasta_to_dataframe(fasta_file_path)
rdict = rdf.set_index("Hash")["Sequence"].to_dict()
del rdf

########################################################
# Differential Abundance Taxa Combining With ASV Sequences
########################################################

vars_in_birdman_analysis = os.listdir(slurm_outputs_dir)

for e, var in enumerate(vars_in_birdman_analysis):
    file_path = f"{slurm_outputs_dir}/{var}/results/differential_taxa.tsv"

    if not os.path.exists(file_path):
        print(var, "differential_taxa.tsv does not exist")
        continue

    diff_abun_df = pd.read_csv(
        file_path,
        sep="\t",
        index_col=0,
    )
    diff_abun_df.rename(
        columns={
            f"{var}_mean": "Variable_mean",
            f"{var}_std": "Variable_std",
            f"{var}_hdi": "Variable_hdi",
        },
        inplace=True,
    )
    diff_abun_df["type"] = diff_abun_df["Variable_mean"].map(
        lambda i: "numerator" if i > 0 else "denominator"
    )
    diff_abun_df["Variable_name"] = var
    diff_abun_df["Sequence"] = diff_abun_df["Feature"].map(rdict)

    if len(diff_abun_df) < 30:
        print(var, "has fewer ASV")
        continue

    if e == 0:
        diff_abun_df_for_comparison = diff_abun_df
    else:
        diff_abun_df_for_comparison = pd.concat(
            [diff_abun_df_for_comparison, diff_abun_df]
        )

    print(var, "done")


diff_abun_df_for_comparison["Variable_name"] = diff_abun_df_for_comparison[
    "Variable_name"
].map(lambda i: i.replace("fg_eaten", "eaten"))
diff_abun_df_for_comparison.to_csv(
    "../outputs/supplementary_files/Supplementary_Data_6.tsv", sep="\t", index=False
)
