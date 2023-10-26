# Assign Putative Tetramer Specificities

Use these scripts to assign putative epitope specificities to each cell using BD Rhapsody output (mRNA counts, TCR, SampleTag calls, etc.), tetramer_dbec.csv, and other meta data (sample_sheet.txt, peptide.csv).

The `seurat.R` script is an example of a single Rhapsody experiment that can be used as a template for analyzing new data in RStudio. This script will output the necessary csv files to run with the totalvi scripts to merge, background correct, batch correct, integrate (RNA + protein), and normalize expression data.

# Dependecies:
Installation of the following R packages is required:<br>

Seurat V3 (https://satijalab.org/seurat/)<br>
sctransform (https://satijalab.org/seurat/articles/sctransform_vignette.html)<br>
ggplot2<br>
cowplot<br>
reshape2<br>
RColorBrewer<br>
dplyr<br>
tidyr<br>
cluster<br>
Matrix<br>
fitdistrplus<br>
stats<br>
cgwtools<br>
tiblle<br>
gt<br>
ggridges<br>
ggrepel<br>
gridExtra<br>
pheatmap<br>
ggpubr<br>

#END
