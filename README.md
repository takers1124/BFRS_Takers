# BFRS_Takers
Data analysis protocol for Takers's MS research at BFRS (DNA amplicon sequencing of fungal ITS1 barcode)

Part 1 - bioinformatics with Python using QIIME2 environment

(1) 
Install QIIME2.
I used version 2023.7 for initial analysis.

(2)
Import raw DNA sequences (which I will send in a zip file via email). 
These sequences had primer & barcodes removed (by Novogene).
Many ways to upload, I used a manifest file (included in zip file). 

(3) 
View sequence QC, denoise with Dada2, construct feature tables.

FeatureTable[Frequency] = count of sequence reads (frequency) for each sample.

FeatureTable[Sequence] = matches DNA sequence with each unique "feature". 

Terms ("feature" = ASV = amplicon sequence variant = OTU = operational taxonomic unit).

(4)
Incorporate metadata (sites, treatments) and generate preliminary visualizations. 

(5) 
Create initial diversity metrics (alpha rarefaction and phylogenetic tree), export FeatureTable & tree for use in R, and run ANCOM-BC analysis.

(6)
Assign taxonomy to ASVs. 
References sequence data to UNITE database with machine learning. 
I had to connect to a remote server (AWS) to do this.  

(7)
Assign functional groups. 
References taxonomic groups (assigned with UNITE) to FunGuild database. 

Part 2 - statistical analysis with R 

(8) 
Export feature table and taxonomy table (with functional groups) as CSV files and import into R.

(9) 
Install required packages

(10)


