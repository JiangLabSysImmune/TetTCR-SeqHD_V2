# Assign Putative Tetramer Specificities

Use these scripts to assign putative epitope specificities to each cell using BD Rhapsody output (mRNA counts, TCR, SampleTag calls, etc.), tetramer_dbec.csv, and other meta data (sample_sheet.txt, peptide.csv).

The seurat.R script is an example of a single Rhapsody experiment that can be used as a template for analyzing new data in RStudio. This script will output the necessary csv files to run with the totalvi scripts to merge, background correct, batch correct, integrate (RNA + protein), and normalize expression data.

# Dependecies:
Installation of the following R packages is required:

Seurat (V3)
ggplot2
sctransform
cowplot
reshape2
RColorBrewer
dplyr
tidyr
cluster
Matrix
fitdistrplus
stats
cgwtools
tiblle
gt
ggridges
ggrepel
gridExtra
pheatmap
ggpubr

#END
