import os, sys, time, argparse
import numpy as np
import pandas as pd
from collections import defaultdict
import scipy.stats
import statsmodels.stats.multitest
import tqdm 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="BS_hisat2",fromfile_prefix_chars='@',description=__doc__,formatter_class=argparse.RawTextHelpFormatter)
    #Required
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i","--input",dest="input",required=True,help="Input gtf file (raw signals)")
    group_required.add_argument("-o","--output",dest="output",required=True,help="Output prefix")
    group_required.add_argument("--n_iters",dest="n_iters", type=int, default=1000, required=False,help="Iteration times,default=10000")
    group_required.add_argument("--p-thresh",dest="p_thresh", type=float, default=0.05, required=False,help="P-value threshold, default=0.05")
    group_required.add_argument("--q-thresh",dest="q_thresh", type=float, default=0.1, required=False,help="Q-value threshold (FDR-BH), default=0.1")
    group_required.add_argument("--min-count",dest="min_count", type=int, default=10, required=False,help="Minimal read counts, default=10")
    # group_required.add_argument("--min-prop",dest="min_prop", type=float, default=0.025, required=False,help="Minimal proportion within a gene, default=0.025")
    group_required.add_argument("--min-rel-prop",dest="min_rel_prop", type=float, default=0.025, required=False,help="Minimal proportion within a gene, default=0.025")
    options = parser.parse_args()

    dict_records = defaultdict(dict)
    # 1. Group all signals by gene
    with open(options.input, "r") as input:
        for line in input.readlines():  # The gtf is under 1 Gb
            # 19      KO_ctrl_1.gtf   tssgtf  58347629        58347629        0.607198311697  -       5       A1BG    A
            line = line.strip().split("\t")
            chr = line[0]
            pos_1 = int(line[2])
            tpm = float(line[7])
            strand = line[5]
            count = int(line[4])
            gene, gene_biotype, trans_biotype = line[3].split("|")
            base = line[6]
            dict_records[(gene, gene_biotype)][(chr, pos_1, strand)] = (count, tpm, base)

    # 2. Within a gene, compute bootstrapping. If only <= 3 signals were detected ...?
    with open(options.output+".ben.passed.bed", "w") as output_passed, \
         open(options.output+".ben.filtered_out.bed", "w") as output_filtered_out:
        # Header here:
        print("chr\tpos_0\tpos_1\tgene\tbase\tstrand\tcount\tTPM\tprop\tsignals\tzscore\tall_reads\tall_tpms\tgene_biotype", file=output_passed) # bed out
        print("chr\tpos_0\tpos_1\tgene\tbase\tstrand\tcount\tTPM\tprop\tsignals\tzscore\tall_reads\tall_tpms\tgene_biotype", file=output_filtered_out) # bed out
        
        genes_to_analyze = dict_records.keys()
        pbar = tqdm.tqdm(genes_to_analyze)

        for gene in pbar:
            keys = []
            read_counts = []
            tpms = []
            bases = []
            info = []

            for key, (count, tpm, base) in dict_records[gene].items():
                keys.append(key)
                read_counts.append(count)
                tpms.append(tpm)
                bases.append(base)
                info.append(key)

            n_items = len(tpms)
            tpms = np.array(tpms)
            tpms_max = tpms.max()
            rel_prop_max = tpms/tpms.max()
            rel_prop = tpms/tpms.sum()
            tpm_total = np.sum(tpms)
            counts_total = np.sum(read_counts)

            if len(tpms) >= 3:
                zscores = (tpms - np.mean(tpms))/np.std(tpms)
                for k, z, tpm, rel1, rel2, count, base in zip(*[info, zscores, tpms, rel_prop_max, rel_prop, read_counts, bases]):
                    if z > 1 and tpm > 0.5 and rel1 >= 0.05 and tpms_max > 1:
                        print("{chr}\t{pos_0}\t{pos_1}\t{gene}\t{base}\t{strand}\t{count}\t{tpm}\t{prop}\t{signal_count}\t{zscore}\t{counts_total}\t{tpm_total}\t{gene_bio}".format(
                            chr=k[0], pos_0=k[1]-1, pos_1=k[1], strand=k[2], gene=gene[0], gene_bio=gene[1],
                            count=count, tpm=tpm, base=base, prop=rel2,
                            signal_count=n_items, zscore=z, tpm_total=tpm_total, counts_total=counts_total),
                            file=output_passed) # bed out
                    else:
                        print("{chr}\t{pos_0}\t{pos_1}\t{gene}\t{base}\t{strand}\t{count}\t{tpm}\t{prop}\t{signal_count}\t{zscore}\t{counts_total}\t{tpm_total}\t{gene_bio}".format(
                            chr=k[0], pos_0=k[1]-1, pos_1=k[1], strand=k[2], gene=gene[0], gene_bio=gene[1],
                            count=count, tpm=tpm, base=base, prop=rel2,
                            signal_count=n_items, zscore=z, tpm_total=tpm_total, counts_total=counts_total),
                            file=output_filtered_out) # bed out
            else:
                for k, tpm, rel1, rel2, count, base in zip(*[info, tpms, rel_prop_max, rel_prop, read_counts, bases]):
                    z = 999
                    if tpm > 0.5 and rel1 >= 0.05 and tpms_max > 1:
                        print("{chr}\t{pos_0}\t{pos_1}\t{gene}\t{base}\t{strand}\t{count}\t{tpm}\t{prop}\t{signal_count}\t{zscore}\t{counts_total}\t{tpm_total}\t{gene_bio}".format(
                            chr=k[0], pos_0=k[1]-1, pos_1=k[1], strand=k[2], gene=gene[0], gene_bio=gene[1],
                            count=count, tpm=tpm, base=base, prop=rel2,
                            signal_count=n_items, zscore=z, tpm_total=tpm_total, counts_total=counts_total),
                            file=output_passed) # bed out
                    else:
                        print("{chr}\t{pos_0}\t{pos_1}\t{gene}\t{base}\t{strand}\t{count}\t{tpm}\t{prop}\t{signal_count}\t{zscore}\t{counts_total}\t{tpm_total}\t{gene_bio}".format(
                            chr=k[0], pos_0=k[1]-1, pos_1=k[1], strand=k[2], gene=gene[0], gene_bio=gene[1],
                            count=count, tpm=tpm, base=base, prop=rel2,
                            signal_count=n_items, zscore=z, tpm_total=tpm_total, counts_total=counts_total),
                            file=output_filtered_out) # bed out