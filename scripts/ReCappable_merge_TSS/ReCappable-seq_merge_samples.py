import pandas as pd
import argparse
from collections import defaultdict
import time
from time import gmtime, strftime
# input: list of samples
# output: merged table

passed_suffix = ".called.ben.passed.bed"
filtered_suffix = ".called.ben.filtered_out.bed"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    #Required
    group_required = parser.add_argument_group("Required")
    # group_required.add_argument("-i","--input",dest="input",required=True,help="Input table")
    group_required.add_argument("-s","--samplesheet",dest="samplesheet",required=True,help="Input samplesheet")
    group_required.add_argument("-o","--output",dest="output",required=False, default="out.csv",help="Output prefix")  # 1. Matrix of counts; matrix of tpm; matrix of details
    group_required.add_argument("--no-ambiguous",dest="ambiguous",required=False, action="store_true", default=False,help="If use ambiguous sites?")  # 1. Matrix of counts; matrix of tpm; matrix of details
    options = parser.parse_args()

    print("Date: {}".format(strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    all_samples = defaultdict(list)
    data_out = {}
    dict_all_sites = {}
    with open(options.samplesheet, "r") as input:
        for line in input.readlines():
            line = line.strip().split("\t") 
            sample = line[0]
            prefix = line[1]
            # replicate = line[2]
            all_samples[sample].append(prefix)
            # 1. get all sites

            data_out[prefix] = {}
            data_out[prefix]["Count"] = {}
            data_out[prefix]["Zscore"] = {}
            data_out[prefix]["all_reads"] = {}
            data_out[prefix]["is_passed"] = {}

            df_temp = pd.read_csv(prefix+passed_suffix, sep="\t", index_col=[0,1,2,3,4,5,13], header=0, low_memory=False)
            for idx, row in df_temp.iterrows():
                if options.ambiguous == True:
                    if row.index[3].startswith("*") == False:
                        dict_all_sites[idx] = 1
                else:
                    dict_all_sites[idx] = 1

    print("Number of sites analyzing: {}".format(len(dict_all_sites)))
    print("Mode: ambiguous gene annotations not allowed? {}".format(options.ambiguous))
    # 2. fetch everything
    with open(options.samplesheet, "r") as input:
        for line in input.readlines():
            line = line.strip().split("\t") 
            sample = line[0]
            prefix = line[1]
            print("Running: {}".format(prefix))
            df_passed = pd.read_csv(prefix+passed_suffix, index_col=[0,1,2,3,4,5,13], sep="\t", header=0, low_memory=False)
            indexes = [i for i in df_passed.index if i in dict_all_sites]
            df_passed = df_passed.loc[indexes]
            print("Sites passed filters: {}".format(df_passed.shape[0]))
            for idx, row in df_passed.iterrows():
                data_out[prefix]["Count"][idx] = row["count"]
                data_out[prefix]["Zscore"][idx] = row["zscore"]
                data_out[prefix]["all_reads"][idx] = row["all_reads"]
                data_out[prefix]["is_passed"][idx] = True

            df_failed = pd.read_csv(prefix+filtered_suffix, index_col=[0,1,2,3,4,5,13], sep="\t", header=0, low_memory=False)
            indexes = [i for i in df_failed.index if i in dict_all_sites]
            df_failed = df_failed.loc[indexes]
            print("Sites failed to pass filters w/ signals: {}".format(df_failed.shape[0]))
            for idx, row in df_failed.iterrows():
                data_out[prefix]["Count"][idx] = row["count"]
                data_out[prefix]["Zscore"][idx] = row["zscore"]
                data_out[prefix]["all_reads"][idx] = row["all_reads"]
                data_out[prefix]["is_passed"][idx] = False
            print("="*20)
    df_out = pd.DataFrame.from_dict({(i,j): data_out[i][j] 
                                    for i in data_out.keys() 
                                    for j in data_out[i].keys()},
                                    orient='columns')
    df_out.to_csv(options.output)