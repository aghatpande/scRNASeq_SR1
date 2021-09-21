# scRNASeq_SR1

scRNASeq project with chronic lymphocytic leukemia (CLL) data from Dr YL @ FCCC

all data is in the "data"" folder. The old script of SR's seurat analysis from gRED is in the "prev_work" folder.

The data supplied is already "filtered_feature_bc_matrix". Are already filtered by cellranger and dont need to be filtered?

There were multiple Rmd files generated. The final analysis is in the files: "CLL_scRNASeq_alt_filter.Rmd" and "DE_between_CLL_ctrl_CLL_treat.Rmd".

The other rmd files are intermediate analyses that helped produce the final analysis. 
The results are stored on the local computer (macbookair) in the results folder of the project

The SessionInfo is stored in the sessioninfo_scrnaseq_07282021.txt file ending in latest date.

The slide deck is in the cll_seurat_analysis_slides_ver2.zip file and is a Keynote file and the full list of 1685 differentially expressed, annotated genes is in "results/CLL_DE_genes.csv" which can be imported into Excel for viewing.

21st Sept 2021: the latest analysis is in "CLL_separate_analysis_30thaug21.Rmd". The latest slide deck is "...ver3" slides