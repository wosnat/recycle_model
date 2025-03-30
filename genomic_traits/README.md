# Genomic traits

What are the genomic traits of the strains used in the experiment with emphasis on those correlating with the mechanisms tested in the models


The traits are taken from the paper:

 Zoccarato, L., Sher, D., Miki, T., Segrè, D., & Grossart, H. P. (2022). 
 A comparative whole-genome approach identifies bacterial traits for marine microbial interactions. 
 Communications Biology, 5(1), 276.

The data in this folder was taken from the supp data and from the accompanying github repo: 
https://github.com/lucaz88/genome_comparison_code



# reanalysis by luca (located in C:\Users\Osnat\OneDrive - University of Haifa\Documents\results\10cc\Luca reannotation)

 I finished the analyses on 13 genomes (all the yellow one except Thiomicrospira sp). I uploaded the results to this link:

https://bokubox.boku.ac.at/#8dc1f076613b67837b952f542c86bd81

Inside the archive you can find the following data:
• Folder ‘4_gnm_analysis’ contains several reports of MAGs:
        ◦ checkM2 for completeness and quality
        ◦ gtdbtk for the taxonomic classification using Genome Database Taxonomy (GTDB)
• Folder ‘5_gnm_annotation’ contains:
        ◦ A Prokka subfolder containing the gene predictions for each MAG in both nucleotide and translated peptide formats, as well as the GFF and GBK files.
        ◦ Outputs of multiple functional annotations: Prokka, KEGG Orthologies (KO), KEGG Modules (KM), KEGG manual (manual annotation of DHPS) and taurine metabolisms), dbCAN for the CAZy enzymes, antiSMASH for secondary metabolites, BioV for membrane transporters and some other specific databases targeting pathways involving DMSP and vibrioferrin metabolisms, phytohormone production.
        ◦ The ‘MASTER_table.tsv’ file which consolidates the results from all annotation tools into a single table.
• Folder ‘6_plots’ contains html iterative plots that may help exploring some of the provided dataset.

Don't mind the foldernames starting with 4_, 5_,.. as the workflow is included in another pipeline.


# list of EC enzymes 
downloded from Brenda https://www.brenda-enzymes.org/download.php
on 30/3/2025

# list of ROS enzymes

Johnson, Lisa A., and Laura A. Hug. "Distribution of reactive oxygen species defense mechanisms across domain bacteria." Free Radical Biology and Medicine 140 (2019): 93-102.‏
