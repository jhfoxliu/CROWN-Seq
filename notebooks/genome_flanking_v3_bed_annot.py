#!/share/public/apps/bin/python

#Jianheng Liu @ Zhanglab, SYSU
#Aug, 2019
#Email: liujh26@mail2.sysu.edu.cn
#Usage: Get flanking sequence for site list

import argparse
from Bio import SeqIO
from Bio.Seq import reverse_complement
import pysam
import os, sys

def read_genome(fn):
	RefSeq = {}
	for seq in SeqIO.parse(fn,"fasta"):
		RefSeq[seq.id] = seq.seq
	return RefSeq

def get_flanking_region_pysam(ref_genome,chr,start,end,strand,flanking=None,up=0,down=0):
	if flanking is not None:
		up = flanking
		down = flanking
	# start 0-based
	# end 1-based

	try:
		if strand == "+":
			seq = ref_genome.fetch(reference=chr, start=start-up, end=end+down) # [start, end) 0-based
			if seq and len(seq) == end - start + up + down:
				return seq.upper()
		elif strand == "-":
			seq = ref_genome.fetch(reference=chr, start=start-down, end=end+up)
			seq = reverse_complement(seq)
			if seq and len(seq) == end - start + up + down:
				return seq.upper()
	except ValueError:
		pass
	except IndexError:
		pass


def handle_region(reference, method, idx_chr=0, idx_start=1, idx_end=2, idx_strand=5):
	with open(options.input,'r') as input:
		line = input.readline()
		while (line):
			row = line.strip().split(options.join_spacer)
			chr = row[idx_chr]
			start = int(row[idx_start]) # 0-based
			end = int(row[idx_end]) # 1-based
			strand = row[idx_strand]
			seq_out = method(reference,chr,start,end,strand,flanking=options.flanking,up=options.up,down=options.down)
			if seq_out:
				if options.add_chr:
					chr = "chr" + chr
				print("{}\t{}".format("\t".join(row), seq_out))
			line = input.readline()


__doc__ = """
Get flanking sequence.

Input:
-I site (default)    chr[\\t]pos_1[\\t]strand
-I site_0            chr[\\t]pos_0[\\t]strand
-I bed               chr[\\t]pos_0[\\t]pos_1[\\t]strand
-I bed6              chr[\\t]pos_1[\\t]pos_1[\\t].[\\t].[\\t]strand
-I region            chr[\\t]start_1[\\t]end_1[\\t]strand
-I slice_bed         chr[\\t]start_0[\\t]start_1[\\t]strand
-I slice_bed6        chr[\\t]start_0[\\t]start_1[\\t].[\\t].[\\t]strand

Output:
fasta (to stdout)
"""
if __name__ == "__main__":
	description = __doc__
	parser = argparse.ArgumentParser(prog="genome_flanking",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	# input
	group_required = parser.add_argument_group("Inputs")
	group_required.add_argument("-i","--input",dest="input", required=False,help="Site list")
	# group_required.add_argument("-I",dest="input_mode", required=False,default = "site",  choices=["site", "bed", "bed6", "region", "slice_bed", "slice_bed6"],help="Input mode, default=site, single site: site, bed, bed6; a region: region, slice_bed, slice_bed6")
	group_required.add_argument("-f","--fasta",dest="fasta",default=None,help="Reference fasta")
	group_required.add_argument("-s",dest="shortcut",default=None,help="shortcuts, prior to -f, use --show-shortcut to view shortcut lists")
	group_required.add_argument("-F",dest="flanking",type=int,help="Flanking bases, prior to -U and -D")
	group_required.add_argument("-U",dest="up",type=int,default=0,help="Upstream bases")
	group_required.add_argument("-D",dest="down",type=int,default=0,help="Downstream bases")
	group_required.add_argument("-J",dest="join_spacer",default="\t",help="spacer for input, default=\"\\t\"")
	# Mode
	group_mode = parser.add_argument_group("Mode")
	group_mode.add_argument("-M",dest="memory",action="store_true",default=False,help="Load fasta to memory to reduce IO")
	group_mode.add_argument("-S",dest="spacer",default="|",help="spacer for output, default=\"|\"")
	group_mode.add_argument("--no-na", dest="nona", default=False, action="store_true", help="Ignore sites without chromosome")
	# Other
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--add-chr",dest="add_chr",default=False, action="store_true",help="Add chr to reference, default=False")
	group_other.add_argument("--ucsc",dest="ucsc",default=False, action="store_true",help="Set ids as UCSC genome browser keys")
	group_other.add_argument("--samtools",dest="samtools",default="samtools",help="samtools path, default is environment")
	group_other.add_argument("--show-shortcut",dest="showShortCut",default=False,action="store_true",help="View shortcuts")
	group_other.add_argument("--version",action="version",version="%(prog)s 3.0")
	options = parser.parse_args()
	
	
	shortcuts = {\
	# "GRCh37": "/share/public1/data/liujh/database/human/ensembl_release75/genome/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa",
	# "GRCh37.75": "/share/public1/data/liujh/database/human/ensembl_release75/genome/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa",
	# "human":"/share/public1/data/liujh/database/human/ensembl_release75/genome/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa",
    "GRCh38":"/home/fox/Database/Genome/Ensembl/GRCh38_r104/Homo_sapiens.GRCh38.dna_sm.primary_assembly.format.fa",
	#"hg19":"/mnt/c/Users/sysul/wsl_documents/Database/hg19/hg19.fa",
	#"hg38":"/mnt/c/Users/sysul/wsl_documents/Database/hg38/hg38.fa",
	"GRCm39": "/home/fox/Database/Genome/Ensembl/GRCm39_r108/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa",
	"mm10": "/home/fox/Database/Genome/Ensembl/GRCm38_r102/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa",
	# "BDGP5.78":"/share/public1/data/liujh/database/fly/BDGP5.78/genome/Drosophila_melanogaster.BDGP5.dna_sm.toplevel.format.fa",
	# "fly":"/share/public1/data/liujh/database/fly/BDGP5.78/genome/Drosophila_melanogaster.BDGP5.dna_sm.toplevel.format.fa",
	# "Zv9.78":"/share/public1/data/liujh/database/zebrafish/zv9_ensembl78/genome/Danio_rerio.Zv9.dna_sm.toplevel.format.fa",
	# "Zv9":"/share/public1/data/liujh/database/zebrafish/zv9_ensembl78/genome/Danio_rerio.Zv9.dna_sm.toplevel.format.fa",
	# "zebrafish":"/share/public1/data/liujh/database/zebrafish/zv9_ensembl78/genome/Danio_rerio.Zv9.dna_sm.toplevel.format.fa",
	# "XT9.1":"/share/public1/data/liujh/database/xenopus_tropicalis/xenTro9.1/genome/XT9_1.fa",
	# "XT":"/share/public1/data/liujh/database/xenopus_tropicalis/xenTro9.1/genome/XT9_1.fa",
	# "XL9.2":"/share/public1/data/liujh/database/xenopus_laevis/xenLae_9.2/genome/XL9_2.fa",
	# "XL":"/share/public1/data/liujh/database/xenopus_laevis/xenLae_9.2/genome/XL9_2.fa",
	}
	
	if options.showShortCut == True:
		for key, values in shortcuts.items():
			print("\t".join([key, values]))
		sys.exit()
	elif options.fasta is None and options.shortcut is None:
		sys.stderr.write("Exit with any reference.\n")
		sys.exit()
	
	if options.shortcut is not None:
		if options.shortcut in shortcuts:
			fasta_file = shortcuts.get(options.shortcut)
		else:
			raise KeyError("shortcut not in list.")
	else:
		fasta_file = options.fasta

	if os.path.isfile(fasta_file+".fai") == False:
		sys.stderr.write("Missing faidx. Try to build it for reference.\n")
		os.system("{samtools} faidx {fasta}".format(samtools=options.samtools, fasta=fasta_file))
	with pysam.FastaFile(fasta_file) as ref_genome:
		handle_region(ref_genome, get_flanking_region_pysam)

