# BFRS_Takers
Data analysis protocol for Takers's MS research at BFRS (DNA amplicon sequencing of fungal ITS1 barcode)

Part 1 - bioinformatics with Python using QIIME2 environment

(1) 
Install QIIME2.
https://docs.qiime2.org/2024.5/install/native/
I used version 2023.7 for initial analysis.

(2)
Import raw DNA sequences (which I will send in a zip file via email). 
These sequences had primer & barcodes removed (by Novogene).
Many ways to upload, we will use a manifest file (included in zip file). 

(3) 
View sequence QC, denoise with Dada2, construct feature tables.
FeatureTable[Frequency] = count of sequence reads (frequency) for each sample.
FeatureTable[Sequence] = matches DNA sequence with each unique "feature". 
Terms ("feature" = ASV = amplicon sequence variant = OTU = operational taxonomic unit).

(4)
Incorporate metadata (sites, treatments) and generate preliminary visualizations. 

(5)
Assign taxonomy to ASVs. 
References my dataset to UNITE database with machine learning. 
I had to connect to a remote server (AWS) to do this. I might just send you the finished csv file and bypass this step. 

(6)
Assign functional groups. 
References my taxonomic groups to FunGuild database. 

Part 2 - statistical analysis with R 

(1) 
Export feature table and taxonomy table (with functional groups) as CSV files and import into R.

(2) 
Install required packages

(3)


