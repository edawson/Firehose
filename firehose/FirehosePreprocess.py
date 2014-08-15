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


tar_dir = "/media/obelisk/tcga/"
cancer_type = "BRCA"
cancer_dir = "/home/eric/data/tcga/BRCA_preprocessing/"

cancers = ["BLCA",
           "BRCA",
           "CESC",
           "COAD",
           "DLBC",
           "HNSC",
           "KICH",
           "KIRC",
           "LAML",
           "PAAD",
           "PRAD",
           "COADREAD",
        ]

dir_patterns = ["miRseq_Mature_Preprocess.Level_3",
                "Mutation_Assessor.Level_4",
                "mRNAseq_Preprocess.Level_3",
                "CopyNumber_Gistic2.Level_4",
                "Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3"
                ]

file_patterns = [".miRseq_mature_RPM_log2.txt",
                 ".uncv2.mRNAseq_RSEM_normalized_log2.txt",
                 ".medianexp.txt",
                 "all_data_by_genes.txt",
                 ".methylation__humanmethylation27"
                ]

tmpdir = cancer_dir + "tmp"


def run_analyses():
    os.chdir(cancer_dir)

    mirna_seq = filename_to_df(cancer_dir + cancer_type + ".miRseq_mature_RPM_log2.txt")
    mirna_seq = process_mirnaseq(mirna_seq)
    zscore_transformed_mirna = z_score_transform(mirna_seq)
    df_to_filename(mirna_seq, cancer_type + ".miRNA_seq_processed.txt")
    df_to_filename(zscore_transformed_mirna, cancer_type + ".miRNA_seq_zscore_processed.txt")

    rna_seq = filename_to_df(cancer_dir + cancer_type + ".uncv2.mRNAseq_RSEM_normalized_log2.txt")
    rna_seq = process_rna_seq(rna_seq)
    z_transformed_rna = z_score_transform(rna_seq)
    df_to_filename(rna_seq, cancer_type + ".mRNA_seq_processed.txt")
    df_to_filename(z_transformed_rna, cancer_type + ".mRNA_seq_zscore_processed.txt")

    try:
        rna_arr = filename_to_df(cancer_dir + cancer_type + ".medianexp.txt")
        rna_arr = process_rna_array(rna_arr)
        df_to_filename(rna_arr, cancer_type + ".mRNA_array_process.txt")
    except IOError:
        print "mRNA Array data not found; continuing..."

    cnv = filename_to_df("all_thresholded.by_genes.txt")
    cnv = process_cnv_file(cnv)
    df_to_filename(cnv, cancer_type + ".cnv_processed.txt")

    # TODO METHYLATION DATA MAY BE 450 or 27...NOT SURE WHICH TO USE
    met = filename_to_df(cancer_dir + cancer_type + ".methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
    met = process_methylation(met)
    df_to_filename(met, cancer_type + ".methylation_processed.txt")
    return rna_seq


def find_firehose_files():
    return


def z_score_transform(df):
    """Transforms an entire dataframe from raw values
    to z-scores based on the mean of the entire dataframe"""
    dat = df.values
    no_nans = dat[~np.isnan(dat)]
    mean = np.mean(no_nans)
    sd = np.std(no_nans)
    ret = df.applymap(lambda x: (x - mean) / sd)

    return ret  # TODO


def print_stats(df, data_type):
    """Prints basic statistics about the passed dataframe
    such as sample number, element number and the number of miRNAs present"""
    print "Data type: ", data_type
    print "Number of primary tumor samples from individual patients: ", len(df.columns)
    print "Number of miRNAs: ", len([x for x in df.index if "hsa-mir" in x.lower() or "hsa-let" in x.lower()])
    print "Number of Genes: ", len(df.index)


def filename_to_df(infile, dframes=None):
    """Takes in a filename for a tab-delimited file and returns
    a pandas DataFrame with the first column of the file as the index.
    An optional list-like object may be passed and the returned datframe
    will be appended to it."""
    f = open(infile, "rU")
    df = pd.read_csv(f, delimiter="\t", index_col=0)
    f.close()
    if dframes is not None:
        dframes.append(df)
    return df


def df_to_filename(df, filename):
    """Write to the passed dataframe to <filename>
    as a tab-separated csv-style file. Floats are limited
    to five decimals of precision for sanity. Missing values are
    written as 'NA'."""
    with open(filename, "w") as f:
        df.to_csv(f, na_rep="NA", float_format="%.5f", sep="\t")


def is_num_or_zero(x):
    try:
        return float(x) != 0
    except ValueError:
        return False


## TODO DOCS
## Processes a Copy Number Variation files
## according to Guanming's specs:
## 1. Remove any non-primary tumour samples
## 2. Rename samples to their 12-character TCGA identifiers
## 3. Drop unknown genes in the sample
## 4. Drop any extraneous rows known to exist in TCGA data
## 5. Print some basic stats for CNV data
def process_cnv_file(df):
    ## Only keep tumor samples (no metastases/normals)
    df = get_only_tumor_samples(df)
    df = rename_samples(df)
    # Remove unknown genes
    df.drop([x for x in df.index if "?" in x], axis=0, inplace=True)

    df = process_file(df)
    print_stats(df, "CNV")
    return df


def rename_samples_for_parity(dframes):
    ret = map(rename_samples, dframes)
    return ret  # TODO


## TODO
## Takes the union of a list of dataframes
## containing genomic data. Dataframes should be in
## element x sample form (e.g. rows=genes and cols=patients)
def union_of_data(dframes):
    dframes = rename_samples_for_parity(dframes)
    ret = pd.concat(dframes, join="outer", axis=1)
    print_stats(ret, "Union")
    return ret


## TODO
## Takes the intersection of a list of dataframes
## containing genomic data. Dataframes should be in
## element x sample form (e.g. rows=genes and cols=patients)
def intersection_of_data(dframes):
    # max_frame = filter(lambda a, b: a if (len(a.columns) > len(b.columns)) else b, dframes)
    # ret_cols = max_frame.columns
    # for d in dframes:
    #     ret_cols = [x for x in ret_cols if x in d.columns]
    # ret = [d[ret_cols] for d in dframes]
    dframes = rename_samples_for_parity(dframes)
    ret = pd.concat(dframes, join="inner", axis=1)
    print_stats(ret, "Intersection")
    return ret


## Takes a single dataframe and returns
## a sliced dataframe containing only columns with
## primary tumor data as indicated by TCGA convention
## (*-01[A-Z]*)
def get_only_tumor_samples(df):
    return df[[x for x in df.columns if (re.match("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01[A-Z]?", x) is not None)]]


## Takes a dataframe and returns a matching
## datframe with the columns renamed to only
## their TCGA center and patient identifier
## (i.e. TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})
def rename_samples(df):
    df = df.rename(columns=lambda x: x[0:12])
    return df


## Takes in a data file containing miRNASeq Data
## and processes it according to Guanming's specs:
## 1. Remove all samples that are not from primary tumors
## 2. Remove any miRNAs that are not human
## 3. Remove the MIMAT identifier from the miRNA names, keeping
## only the standard miRNA names.
## (NOTA BENE: no miRBase renaming is performed)
## 4. Rename tumour samples to their 12-character center+sample
## TCGA identifiers.
## 5. Remove and additional header lines that are known to exist in
## Firehose files.
## 6. Print some basic stats on the file.
def process_mirnaseq(df):
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


# TODO need to get max vals for genes
## Takes a single dataframe extracted from a TCGA methylation file
## and performs transformations according to Guanming's specs:
## 1. Remove all genomic location/chromosome information, maintaing
## only beta values
## 2. Reindex beta values by gene symbol
## 3. Drop any values which map to multiple genes or which map to an
## unknown gene
## 4. Drop any sample or genes which have no data
## 5. Remove any non-primary tumour samples
## 6. Rename tumour samples to their 12-character TCGA identifier
## 7. Remove any extraneous header lines known to exist in Firehose files
## 8. Print some basic stats
def process_methylation(df):
    # Keep only columns referring to Beta Values
    col_to_bool = {x[0]: x[1] for x in zip(df.columns, np.array(df.ix["Composite Element REF"] == "Beta_value"))}
    betas = df[[x for x in col_to_bool if col_to_bool[x]]]

    col_to_bool = {x[0]: x[1] for x in zip(df.columns, np.array(df.ix["Composite Element REF"] == "Gene_Symbol"))}
    genes = df[[x for x in col_to_bool if col_to_bool[x]][1]]
    genes = genes.fillna("NA")
    betas.index = genes

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


# TODO
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


#TODO
def process_rna_array(df):
    df = get_only_tumor_samples(df)
    df = df.drop([x for x in df.index if "?" in x], axis=0)
    df.rename(index={x: x.split("|")[0] if not "SLC35E2|728661" in x else "SLC35E2B" for x in df.index}, inplace=True)
    df = rename_samples(df)
    df = df.dropna(axis=[0, 1], how="all")
    df = process_file(df)
    print_stats(df, "messenger RNA Microarray")
    return df


def process_rppa(df):
    df = get_only_tumor_samples(df)
    df = df.drop([x for x in df.index if "?" in x], axis=0)
    df = df.rename(index={x: x.split("|")[0] if not "SLC35E2|728661" in x else "SLC35E2B" for x in df.index})
    df = rename_samples(df)
    df = df.dropna(axis=[0, 1], how="all")
    df = process_file(df)
    print_stats(df, "Reverse-Phase Protein Array")
    return df


def merge_tumor_dfs(dframes):
    return pd.concat(dframes)


def process_file(df):
    drops = ["Composite Element REF"]
    df = df[[col for col in df.columns if "TCGA" in col]]
    try:
        df = df.drop(drops)
    except ValueError:
        pass
    df = df.dropna(axis=[0, 1], how="all")
    return df


if __name__ == "__main__":

    ## Check CNV, miRNA, mRNA, methylation files
    # os.chdir(cancer_dir)
    #
    # mirna_seq = filename_to_df(cancer_dir + cancer_type + ".miRseq_mature_RPM_log2.txt")
    # mirna_seq = process_mirnaseq(mirna_seq)
    # zscore_transformed_mirna = z_score_transform(mirna_seq)
    # df_to_filename(mirna_seq, cancer_type + ".miRNA_seq_processed.txt")
    # df_to_filename(zscore_transformed_mirna, cancer_type + ".miRNA_seq_zscore_processed.txt")
    #
    rna_seq = filename_to_df(cancer_dir + cancer_type + ".uncv2.mRNAseq_RSEM_normalized_log2.txt")
    rna_seq = process_rna_seq(rna_seq)
    df_to_filename(rna_seq, "out.txt")
    print len(rna_seq.index), " UNIQUE: ", len(set(rna_seq.index))
    # z_transformed_rna = z_score_transform(rna_seq)
    # df_to_filename(rna_seq, cancer_type + ".mRNA_seq_processed.txt")
    # df_to_filename(z_transformed_rna, cancer_type + ".mRNA_seq_zscore_processed.txt")
    #
    # rna_arr = filename_to_df(cancer_dir + cancer_type + ".medianexp.txt")
    # rna_arr = process_rna_array(rna_arr)
    # df_to_filename(rna_arr, cancer_type + ".mRNA_array_process.txt")
    #
    # cnv = filename_to_df("all_data_by_genes.txt")
    # cnv = process_cnv_file(cnv)
    # df_to_filename(cnv, cancer_type + ".cnv_processed.txt")
    #
    # met = filename_to_df("cancer_dir + cancer_type + ".methylation__humanmethylation27__jhu" +
    # "_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
    # met = process_methylation(met)
    # df_to_filename(met, cancer_type + ".methylation_processed.txt")
