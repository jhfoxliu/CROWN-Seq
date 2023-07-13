import sys
import numpy as np
import pandas as pd
import argparse
from collections import defaultdict
import time
from time import gmtime, strftime
# input: list of samples
# output: merged table

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
    group_required.add_argument("-o","--output",dest="output",required=False, default="out.csv",help="Output prefix")  # 1. Matrix of counts; matrix of tpm; matrix of details
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
    df_count = df_init.xs("Count", axis=1 , level=1)

    for sample, replicates in all_samples.items():
        data_out[sample] = {}
        data_out[sample]["Count"] = {}
        data_out[sample]["TPM"] = {}
        data_out[sample]["Fraction"] = {}
        data_out[sample]["Passed"] = {}
        subdf = df_count[replicates]
        for idx, row in subdf.iterrows():
            data_out[sample]["Count"][idx] = int(row.sum())
        
    df_out = pd.DataFrame.from_dict({(i,j): data_out[i][j] 
                                    for i in data_out.keys() 
                                    for j in data_out[i].keys()},
                                    orient='columns')
    df_out[("INFO", "Gene")] = df_out.index.get_level_values(3)
    df_out[("INFO", "Base")] = df_out.index.get_level_values(4)
    df_out[("INFO", "Biotype")] = df_out.index.get_level_values(5)
    count_sample = df_out.xs("Count", axis=1, level=1)
    count_sample["Gene"] = count_sample.index.get_level_values(3)
    for sample in all_samples.keys():
        # pre-filter
        df_out[(sample, "Fraction")] = count_sample.groupby(by="Gene")[sample].transform(lambda x: x/x.sum())
        # df_out[(sample, "Fraction")] = df_out[(sample, "Fraction")].fillna(None)
        df_out[(sample, "TPM")] = df_out[(sample, "Count")] / (df_out[(sample, "Count")].sum() + 0.) * 1000000
        df_out[(sample, "Passed")] = False
        df_out.loc[df_out[(sample, "TPM")]>= 0.1,  (sample, "Passed")]  = True
    df_out.to_csv(options.output + ".merged.csv")
    df_out.xs("Fraction", axis=1, level=1).to_csv(options.output + ".xs_frac.csv")
    df_out.xs("TPM", axis=1, level=1).to_csv(options.output + ".xs_tpm.csv")
    df_out.xs("Count", axis=1, level=1).to_csv(options.output + ".xs_count.csv")

    xs_df_out = df_out.xs("Count", axis=1, level=1)
    xs_df_out["Base"] = xs_df_out.index.get_level_values(4)
    xs_df_out.groupby(by="Base")[list(all_samples.keys())].sum().to_csv(options.output+".stats.csv")