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
                        help="A CNV array file for processing.")
    return parser.parse_args()


def process_cnv_file(df):
    """Processes a CNV file:
    1. Removes any non primary tumour samples.
    2. Renames samples to their 12-character TCGA identifiers
    3. Drops unknown genes from analysis
    4. Drops any extraneous rows known to exist in Firehose data
    5. Prints some basic information about the data."""
    ## Only keep tumor samples (no metastases/normals)
    df = get_only_tumor_samples(df)
    df = rename_samples(df)
    # Remove unknown genes
    df.drop([x for x in df.index if "?" in x], axis=0, inplace=True)

    df = process_file(df)
    print_stats(df, "CNV")
    return df

if __name__ == "__main__":
    args = parse_args()
    df = filename_to_df(args.infile)
    process_cnv_file(df)
    out_name = ".".join([".".join((args.infile).split(".")[:-1]), "processed.txt"])
    df_to_filename(df, out_name)