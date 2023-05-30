import sys
from Bio import SeqIO
from Bio.Seq import reverse_complement

# cmd
# python <script> input_table ref_genome output_control_list output_site_control_pairs

reference_genome = {}
for seq in SeqIO.parse(sys.argv[2],"fasta"):
    reference_genome[seq.id] = str(seq.seq).upper()

black_list = {}
with open(sys.argv[1], "r") as input:
    for line in input.readlines():
        chr, pos, strand = line.strip().split("\t")
        pos = int(pos)
        pos_0 = pos - 1
        key = (chr, pos_0, strand)
        black_list[key] = 1

with open(sys.argv[1], "r") as input, open(sys.argv[3], "w") as output_ctrl, open(sys.argv[4], "w") as output:
    for line in input.readlines():
        chr, pos, strand = line.strip().split("\t")
        pos = int(pos)
        pos_0 = pos - 1
        pos_test = pos_0 

        ctrl_pos = None
        if strand == "+":
            while ctrl_pos is None:
                pos_test = pos_test + 1
                key = (chr, pos_test, strand)
                if key not in black_list:
                    
                    base = reference_genome[chr][pos_test]
                    if base == "A":
                        ctrl_pos = pos_test
        else:
            while ctrl_pos is None:
                pos_test = pos_test - 1
                key = (chr, pos_test, strand)
                if key not in black_list:
                    base = reference_genome[chr][pos_test]
                    if base == "T":
                        ctrl_pos = pos_test
        print("{}\t{}\t{}".format(chr, ctrl_pos + 1, strand), file=output_ctrl)
        print("{}\t{}\t{}\t{}\t{}\t{}".format(chr, pos, strand, chr, ctrl_pos + 1, strand), file=output)