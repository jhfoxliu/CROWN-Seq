import sys

total_count = 0.0
with open(sys.argv[1], "r") as input:
    for line in input.readlines():
        line = line.strip().split("\t")
        count = int(line[4])
        total_count += count

with open(sys.argv[1], "r") as input:
    for line in input.readlines():
        line = line.strip().split("\t")
        count = int(line[4])
        tpm = count/total_count * 1000000
        print("{}\t{}".format("\t".join(line), tpm ))