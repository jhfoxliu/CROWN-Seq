import pysam
import sys
from collections import defaultdict

# Single end
data = defaultdict(int)
N = 0
with pysam.AlignmentFile(sys.argv[1], "rb") as SAM:
    for read in SAM:
        if read.get_tag("NH") == 1 and read.is_secondary == False and read.is_supplementary == False and read.is_unmapped == False:
            if read.is_read1 and read.is_reverse == False:
                align_pairs = read.get_aligned_pairs(matches_only=False)
                useful_base = align_pairs[0]
                if useful_base[1] is not None:
                    key = (read.reference_name, useful_base[1], "+") # 0-based 
                    data[key] += 1
                    N += 1
                    sys.stderr.write("{}\r".format(N))
            elif read.is_read1 and read.is_reverse == True:
                align_pairs = read.get_aligned_pairs(matches_only=False)
                useful_base = align_pairs[-1]
                if useful_base[1] is not None:
                    key = (read.reference_name, useful_base[1], "-") # 0-based 
                    data[key] += 1
                    N += 1
                    # sys.stderr.write("{}\r".format(N))
for key, values in data.items():
    print("{}\t{}\t{}\t{}\t{}\t{}".format(key[0], key[1], key[1] + 1, ".", values, key[2]))
