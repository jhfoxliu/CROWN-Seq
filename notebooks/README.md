## Description

### Workflow: 

(1) *Denovo* TSS call from CROWN-seq data
(2) Merge TSS list from different runs
(3) Use the *denovo* called site list to fetch all 5' end reads at TSS

(1) and (2) are optional if you can provide a site list)


### Notebook orders:

(1) TSS call:

Inputs: 

* CROWN-seq reads

Notebooks:

* For HEK293T (WT): `CROWN-Seq_find_TSS_HEK_WT.ipynb`
* For A549 (WT): `CROWN-Seq_find_TSS_A549_WT.ipynb`

Outputs:

* TSS sites passed filters, BED file (e.g., `HEK_KO.called.ben.passed.bed`)
* TSS sites failed to pass the filters, BED file (e.g., `HEK_WT.called.ben.filtered_out.bed`)

(2) Merge TSS list:

Input: 

* Bed files from (1)

Notebook:

* See notebook in ./merge/ folder

Outputs:

* TSS site list (e.g. `CROWN_sites.txt`)

(3) Run m6Am call:

* For HEK293T (WT): `CROWN-Seq_measure_m6Am_HEK293T_WT.ipynb`
* For A549 (WT): `CROWN-Seq_measure_m6Am_A549_WT.ipynb` 

Outputs:

* A pileup table for A/U/C/G counts (e.g., `HEK293T_WT.CROWN.csv`)

(4) Downstream analysis is based on the pileups.

