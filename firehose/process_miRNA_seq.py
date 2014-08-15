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
from FirehosePreprocess import z_score_transform


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="infile", default=None, type=str, required=True,
                        help="A microRNA seq file from the Broad Firehose for preprocessing.")
    return parser.parse_args()


def process_mirnaseq(df):
    """Takes in a data file containing miRNAseq data
    and processes it as follows:
    1. Remove any non-human miRNAs
    2. Remove any samples not from primary tumors
    3. Remove the MIMAT identifier from miRNA names
    4. Rename tumour samples to their 12-character TCGA sample identifier
    5. Remove any extraneous header lines known to exist in TCGA data
    6. Print some basic stats and write to outfile
    7. Z-score transform the data and write out as a separate file"""

    ## Only keep tumor samples (no metastases/normals)
    df = get_only_tumor_samples(df)

    ## Remove any non human samples
    df = df.drop([x for x in df.index if "hsa" not in x], axis=0)

    ## Change labels to their shorter identifier
    df = df.rename(index={x: x.split("|")[0] for x in df.index})

    df = rename_samples(df)
    df = process_file(df)
    print_stats(df, "microRNA-Seq")
    return df

if __name__ == "__main__":
    args = parse_args()
    dat = filename_to_df(args.infile)
    dat = process_mirnaseq(dat)
    out_name = ".".join([".".join((args.infile).split(".")[:-1]), "processed.txt"])
    df_to_filename(dat, out_name)

    zed = z_score_transform(dat)
    out_name = ".".join([".".join((args.infile).split(".")[:-1]), "zscore.processed.txt"])
    df_to_filename(zed, out_name)