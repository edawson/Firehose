__author__ = 'Eric T Dawson'

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
from Firehose.FirehosePreprocess import rename_samples
from Firehose.FirehosePreprocess import get_only_tumor_samples
from Firehose.FirehosePreprocess import process_file
from Firehose.FirehosePreprocess import df_to_filename
from Firehose.FirehosePreprocess import print_stats
from Firehose.FirehosePreprocess import filename_to_df
from FirehosePreprocess import z_score_transform


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="infile", default=None, type=str, required=True,
                        help="An Illumina methylation27 or 450 microarray file for preprocessing.")
    return parser.parse_args()


def process_rna_seq(df):
    ## Only keep tumor samples (no metastases/normals)
    df = get_only_tumor_samples(df)
    ## Remove ? genes
    df.drop([x for x in df.index if "?" in x], axis=0, inplace=True)

    ## Change labels and
    ## Fix SLC misnomer
    # SLC35E2|728661 to SLC35E2B|728661
    df = df.rename(index={x: x.split("|")[0] if not "SLC35E2|728661" == x else "SLC35E2B" for x in df.index})
    df = rename_samples(df)
    df = process_file(df)
    df = df.dropna(axis=[0, 1], how="all")
    print_stats(df, "messenger RNA Seq")
    return df

if __name__ == "__main__":
    args = parse_args()
    dat = filename_to_df(args.infile)
    dat = process_rna_seq(dat)
    out_name = ".".join([".".join((args.infile).split(".")[:-1]), "processed.txt"])
    df_to_filename(dat, out_name)

    zed = z_score_transform(dat)
    out_name = ".".join([".".join((args.infile).split(".")[:-1]), "zscore.processed.txt"])
    df_to_filename(zed, out_name)

