__author__ = 'Eric T Dawson'
import os
import sys
import re
import argparse
import multiprocessing as mp
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
from FirehosePreprocess import z_score_transform


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="infile", default=None, type=str, required=True,
                        help="An mRNA microarray file from the Broad Firehose for preprocessing "
                             + ", preferably from an Affymetrix 8x15k or similar array.")
    return parser.parse_args()


def process_rna_array(df):
    """Transforms an Affymetrix mRNA microarray from the Broad Firehose
    for easier analysis. """
    df = get_only_tumor_samples(df)
    df = df.drop([x for x in df.index if "?" in x], axis=0)
    df.rename(index={x: x.split("|")[0] if not "SLC35E2|728661" in x else "SLC35E2B" for x in df.index}, inplace=True)
    df = rename_samples(df)
    df = df.dropna(axis=[0, 1], how="all")

    df = process_file(df)
    print_stats(df, "messenger RNA Microarray")
    return df

if __name__ == "__main__":
    args = parse_args()
    dat = filename_to_df(args.infile)
    dat = process_rna_array(dat)
    out_name = ".".join([".".join((args.infile).split(".")[:-1]), "processed.txt"])
    df_to_filename(dat, out_name)

    # zed = z_score_transform(dat)
    # out_name = ".".join([".".join((args.infile).split(".")[:-1]), "zscore.processed.txt"])
    # df_to_filename(zed, out_name)