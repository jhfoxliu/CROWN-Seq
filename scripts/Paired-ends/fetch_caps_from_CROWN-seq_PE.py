from Bio import SeqIO
import pandas as pd
import argparse
import pysam
import pandas as pd
from Bio.Seq import reverse_complement

if __name__ == "__main__":
    description = """pileup_tRNA_bases"""
    parser = argparse.ArgumentParser(prog="pileup_tRNA_bases",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-r","--ref",dest="reference",required=True,help="reference fasta")
    group_required.add_argument("-l","--list",dest="site_list",required=True,help="site list")
    group_required.add_argument("-b","--bam",dest="bam",required=True,help="input bam, sorted")
    group_required.add_argument("-o","--output",dest="output",required=True,help="output")
    # group_optional = parser.add_argument_group("Optional")
    options = parser.parse_args()

    reference = options.reference
    BAM = options.bam
    output = {"Ref":{}, "A": {}, "T": {}, "C":{}, "G": {}}
    reference_genome = {}
    for seq in SeqIO.parse(reference,"fasta"):
        reference_genome[seq.id] = str(seq.seq).upper()

    with open(options.site_list, "r") as input_sites, pysam.AlignmentFile(BAM, "rb") as input_BAM:
        for line in input_sites.readlines():
            line = line.strip().split("\t")
            chr, pos_1, strand = line
            if chr == "M":
                chr = "MT"
            pos_1 = int(pos_1)
            KEY = (chr, pos_1, strand)
            output["A"][KEY] = 0
            output["T"][KEY] = 0
            output["C"][KEY] = 0
            output["G"][KEY] = 0
            if strand == "+":
                output["Ref"][KEY] = reference_genome[chr][pos_1 -1]
            else:
                output["Ref"][KEY] = reverse_complement(reference_genome[chr][pos_1 -1])

            for read in input_BAM.fetch(chr, pos_1-1, pos_1, multiple_iterators=True):  # , flag_filter=1536
                alignment_pairs = read.get_aligned_pairs(with_seq=True)
                if read.is_read1 == True and strand == "+" and read.is_reverse == False:
                    used = alignment_pairs[0]
                    query_pos_0 = used[0]
                    ref_pos_0 = used[1]
                    query_base = read.query_sequence[0]
                    if ref_pos_0 == pos_1 - 1:
                        if query_base == "A":
                            output["A"][KEY] += 1
                        elif query_base == "T":
                            output["T"][KEY] += 1
                        elif query_base == "C":
                            output["C"][KEY] += 1
                        elif query_base == "G":
                            output["G"][KEY] += 1
                elif read.is_read1 == True and strand == "-" and read.is_reverse == True:
                    used = alignment_pairs[-1]
                    query_pos_0 = used[0]
                    ref_pos_0 = used[1]
                    query_base = read.query_sequence[-1]
                    query_base = reverse_complement(query_base)
                    if ref_pos_0 == pos_1 - 1:
                        if query_base == "A":
                            output["A"][KEY] += 1
                        elif query_base == "T":
                            output["T"][KEY] += 1
                        elif query_base == "C":
                            output["C"][KEY] += 1
                        elif query_base == "G":
                            output["G"][KEY] += 1

df = pd.DataFrame(output)
df.to_csv(options.output)
