__author__ = 'Eric T Dawson'
import argparse

try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError:
    print "You'll need to install PANDAS, NumPy and Matplotlib to utilize these scripts"
    print "This can be done by: "
    print "sudo easy_install pip"
    print "pip install numpy"
    print "pip install matplotlib"
    print "pip install pandas"
from FirehosePreprocess import rename_samples
from FirehosePreprocess import get_only_tumor_samples
from FirehosePreprocess import process_file
from FirehosePreprocess import df_to_filename
from FirehosePreprocess import print_stats
from FirehosePreprocess import filename_to_df


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="infile", default=None, type=str, required=True,
                        help="An Illumina methylation27 or 450 microarray file for preprocessing.")
    return parser.parse_args()


def process_methylation(df):
    """Takes in a pandas DataFrame and performs a variety of transformations
    to make the data easier to work with:
    1. Removes all genomic location/chromosome information
    2. Reindexes the beta values by gene symbol
    3. Removes any unknown genes or those probes matching multiple genes
    4. Drops any genes of samples with no data
    5. Removes any non-primary tumour samples
    6. Renames tumour samples to the first 12 characters of their TCGA identifier
    7. Removes any extraneous header lines
    8. Prints some basic stats about the DataFrame."""

    print df
    exit(0)
    ## Keep only columns referring to Beta Values
    col_to_bool = {x[0]: x[1] for x in zip(df.columns, np.array(df.ix["Composite Element REF"] == "Beta_value"))}
    betas = df[[x for x in col_to_bool if col_to_bool[x]]]

    ## Grab the gene symbol column, fill any missing gene names with "NA" and make the gene symbol column
    ## the index of the new beta value dataframe
    col_to_bool = {x[0]: x[1] for x in zip(df.columns, np.array(df.ix["Composite Element REF"] == "Gene_Symbol"))}
    genes = df[[x for x in col_to_bool if col_to_bool[x]][1]]
    genes = genes.fillna("NA")

    betas.index = genes

    ## Change the name of the beta value index to "Gene" rather than "Gene Symbol"
    ## Drop any samples or genes with no data.
    betas.index.name = "Gene"
    betas = betas.dropna(how="all", axis=(0, 1))
    betas = betas.drop([x for x in betas.index if x == "NA" or ";" in x or "?" in x])

    df = betas
    df.index.name = "Gene"
    df = get_only_tumor_samples(df)
    df = rename_samples(df)
    df = process_file(df)

    print_stats(df, "Methylation")
    return df

if __name__ == "__main__":
    args = parse_args()
    dat = filename_to_df(args.infile)
    dat = process_methylation(dat)
    out_name = ".".join([".".join((args.infile).split(".")[:-1]), "processed.txt"])
    df_to_filename(dat, out_name)

