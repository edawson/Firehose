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
                        help="An NCI MAF file for processing.")
    return parser.parse_args()