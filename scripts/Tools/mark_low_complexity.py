import sys
from Bio import SeqIO
from Bio.Seq import reverse_complement
import regex
reference = {}
for seq in SeqIO.parse(sys.argv[1], "fasta"):
    reference[seq.id] = str(seq.seq).upper()

with open(sys.argv[2], "r") as input:
    for line in input.readlines():
        line = line.strip().split("\t")
        chr, pos, strand = line
        reject = False
        pos = int(pos)
        if reference.get(chr) is None:
            print("{}\t{}\t{}\tUnknown".format(chr, pos, strand))
        if strand == "+":
            flanking = reference[chr][pos-11: pos+10].replace("A", "G")
            res = regex.finditer("(GGGGGGGG){e<=1}", flanking, overlapped=True)
            for i in res:
                if i.start() <= 10 and i.end() >= 10:
                    reject = True
                    break
            print("{}\t{}\t{}\t{}".format(chr, pos, strand, reject))
        else:
            flanking = reverse_complement(reference[chr][pos-11: pos+10]).replace("A", "G")
            res = regex.finditer("(GGGGGGGG){e<=1}", flanking, overlapped=True)            
            for i in res:
                if i.start() <= 10 and i.end() >= 10:
                    reject = True
                    break
            print("{}\t{}\t{}\t{}".format(chr, pos, strand, reject))