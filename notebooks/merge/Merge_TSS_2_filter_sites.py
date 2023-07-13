import sys
import numpy as np
import pandas as pd
import argparse
from collections import defaultdict
import time
from time import gmtime, strftime
# input: list of samples
# output: merged table
import warnings
warnings.filterwarnings('ignore')
passed_suffix = ".called.ben.passed.bed"
filtered_suffix = ".called.ben.filtered_out.bed"

def Bootstrpping(x):
    print(x)
    print(x.shape)
    print("==========")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    #Required
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i","--input",dest="input",required=True,help="Input table")
    group_required.add_argument("-s","--samplesheet",dest="samplesheet",required=True,help="Input samplesheet")
    group_required.add_argument("-o","--output",dest="output",required=False, default="out.csv",help="Output")  # 1. Matrix of counts; matrix of tpm; matrix of details
    group_required.add_argument("--no-ambiguous",dest="ambiguous",required=False, action="store_true", default=False,help="If use ambiguous sites?")  # 1. Matrix of counts; matrix of tpm; matrix of details
    options = parser.parse_args()

    print("Date: {}".format(strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    all_samples = defaultdict(list)
    data_out = {}
    dict_all_filtered_sites = {}
    with open(options.samplesheet, "r") as input:
        for line in input.readlines():
            line = line.strip().split("\t") 
            sample = line[0]
            prefix = line[1]
            # replicate = line[2]
            all_samples[sample].append(prefix)

    df_init = pd.read_csv(options.input, index_col=[0,1,2,3,4,5,6], sep=",", header=[0,1], low_memory=False)
    if options.ambiguous == True:
        indexes = []
        for idx in df_init.index:
            if idx[3].startswith("*") == False:
                indexes.append(idx)
        df_init = df_init.loc[indexes]
    
    used_indexes = []
    N = 0
    for sample, replicates in all_samples.items():
        if len(replicates)==1:
            subdf = df_init[[i for i in df_init.columns if i[0] in replicates]]
            xs_passed = subdf.xs("is_passed", axis=1, level=1).fillna(False)
            indexes = xs_passed[xs_passed.any(axis=1)==True].index
            used_indexes.extend(indexes)
        else:
            subdf = df_init[[i for i in df_init.columns if i[0] in replicates]]
            xs_passed = subdf.xs("is_passed", axis=1, level=1).fillna(False)
            xs_passed["Gene"] = xs_passed.index.get_level_values(3)
            # if xs_passed.all(axis=1).all() == True:
            #     used_indexes.extend(used_indexes)
            # else:
                # indexes = xs_passed[xs_passed.any(axis=1)==True].index
                # subdf = subdf.loc[indexes]
            xs_counts = subdf.xs("Count", axis=1, level=1).fillna(0)
            xs_counts["Gene"] = xs_counts.index.get_level_values(3)
            for gene in xs_counts["Gene"].unique():
                xs_counts_gene = xs_counts[xs_counts["Gene"] == gene][replicates]
                xs_passed_gene = xs_passed[xs_passed["Gene"] == gene][replicates]
                if xs_passed_gene.all(axis=1).all() == True:
                     used_indexes.extend(xs_passed_gene.index)
                else:
                    xs_pseudocounts = xs_counts_gene.copy()
                    coverage_sums = xs_counts_gene.sum(axis=1)
                    freq = coverage_sums/coverage_sums.sum()
                    sample_sums = xs_counts_gene.sum(axis=0)
                    for sample in replicates:
                        xs_pseudocounts[sample] = freq * sample_sums[sample]
                    xs_pseudocounts_lower = np.round(xs_pseudocounts * 0.25, 0)
                    xs_pseudocounts_upper = np.round(xs_pseudocounts * 4, 0)
                    # print(xs_pseudocounts)
                    xs_pseudocounts_out = (xs_counts_gene>=xs_pseudocounts_lower) & (xs_counts_gene<=xs_pseudocounts_upper)
                    xs_pseudocounts_out_filtered = xs_pseudocounts_out[xs_pseudocounts_out.all(axis=1) == True]

                    # cusum
                    xs_counts_gene["Sum"] = xs_counts_gene.sum(axis=1)
                    xs_counts_gene = xs_counts_gene.sort_values(by="Sum", ascending=False)
                    xs_counts_gene["Cusum"] = xs_counts_gene["Sum"].cumsum()/xs_counts_gene["Sum"].sum()
                    xs_pseudocounts_out_filtered["Cusum"] = xs_counts_gene["Cusum"]
                    # if xs_pseudocounts_out_filtered[xs_pseudocounts_out_filtered["Cusum"]>0.95].shape[0] > 1:
                    #     used_indexes.extend(xs_pseudocounts_out_filtered.loc[xs_pseudocounts_out_filtered["Cusum"]<=0.95].index)
                    # else:
                    used_indexes.extend(xs_pseudocounts_out_filtered.index)
                N += 1
                sys.stderr.write("{}\r".format(N))
    df_out = df_init.loc[df_init.index.intersection(used_indexes)]
    df_out.to_csv(options.output)
    print("Finished at: {}".format(strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

