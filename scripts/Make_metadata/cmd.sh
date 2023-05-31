# convert gtf file to annotation file
python gtf2anno_plus.py -i {name}.gtf > {name}.anno

# extract tss information from anno file
python anno_to_tss -i {name}.anno > {name}.tss

# convert gtf file to genelist file
python gtf2genelist.py -i {name}.gtf -f transcripts.fa > {name}.genelist

# extract transcription information from genelist

awk 'OFS="\t"{print $6, $7-1, $8, $2, $3, $9}' {name}.genelist | bedtools sort -i | lss > {name}.tx.bed