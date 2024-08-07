# LIFE4137-Individual-Project
Anchored GWAS Identifies Pathways Linked to Leaf and Seed Ionome in 1,164 Arabidopsis Accessions

## Workflow
### 1. Run GWAS pipeline on Ada
e.g. sbatch GWAS_leaf_Ni60.sh
### 2. Retrieve the top 10 significant SNPs and files
### 3. Annotate the likely region of candidate genes
+ #### Create 500 bp windows around the peaks <br>
get_region.R
+ #### Intersect the regions with gff annotation file on Linux or search on genome browser
```
# convert to bed files for running bedrolls
awk 'OFS="\t"{print $1,$2,$3,".",".","."}' leaf_As75_peaks.txt>leaf_As75_peaks.bed

# intersect the peak regions with annotation file
bedtools intersect -wb -a leaf_As75_peaks.bed -b TAIR10_GFF3_genes.gff > leaf_As75_annotated_peaks.gff
```
#### Genome Browser <br>
https://jbrowse2.arabidopsis.org/index.html?session=local-xlYTiOwW-6_kJOC6JSWTx
### 4. Get the gene descriptions from TAIR 
https://www.arabidopsis.org/
### 5. Visualization of global and anchored GWAS results
Global_Manhatton.R <br>
accession_Manhatton.R
