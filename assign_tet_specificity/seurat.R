# example code to pre-process individual experiments before batch correction
# main function is to assign putative specificities to cells/clonotypes (tetDemux)
# each experiment is slightly modified based on quality thresholds, etc.
# this script is for 22-11-09 experiment
# modify paths, etc. before use

library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(reshape)
library(RColorBrewer)
library(dplyr)
library(tidyr)

setwd('C:/Users/mjmalone/Desktop/run.196/mjmalone_20221109_SVAR/rhapsody/')
source("C:/Users/mjmalone/Desktop/scripts/seurat_functions.R")
source("C:/Users/mjmalone/Desktop/scripts/tetDemux.R")

experiment_id <- "S-COV_22-11-09"

#specify files
tcr_file <- 'JiLL2850_VDJ_perCell.csv'
rhapsody_file <- 'Combined_JiLL2850_DBEC_MolsPerCell.csv'
st_file <- "JiLL2850_Sample_Tag_Calls.csv"
tetramer_file <- "tetramer_dbec.csv"
sample_sheet_file <- "sample_sheet.txt"
peptide_file <- "peptide.csv"

#load files
rhapsody <- read.csv(rhapsody_file,comment.char = "#",row.names = 1,check.names = FALSE)
sampletag <- read.csv(st_file,comment.char = "#",row.names = 1)
tcr <- read.csv(tcr_file,comment.char = "#",row.names = 1)
tet <- read.csv(tetramer_file,row.names = 1,check.names = FALSE)
sample_sheet <- read.csv(sample_sheet_file,sep = '\t')
peptide <- read.csv(peptide_file)

#meta_data
meta_data <- sampletag
meta_data <- gather.meta(meta_data,sample_sheet,experiment_id)

#remove multiplets and undetermined sample tag calls
rhapsody <- rhapsody[rownames(sampletag[!grepl('Multiplet|Undetermined',sampletag$Sample_Tag),]),]
##### remove low umi count cells (need to adjust manually)
hist(log10(rowSums(rhapsody)),breaks=50)
#define cutoff
thresh <- 2.6
abline(v=thresh)
cutoff <- 10 ^ thresh
rhapsody <- rhapsody[rowSums(rhapsody)>=cutoff,]
meta_data <- meta_data[rownames(rhapsody),]
#tetramer (match gene expr rownames)
tet <- tet[rownames(rhapsody),]
rownames(tet) <- rownames(rhapsody)
tet[is.na(tet)] <- 0
#add missing any missing peptides from ref
if (length(setdiff(peptide$name,colnames(tet))) > 0){
  tet[,setdiff(peptide$name,colnames(tet))] <- 0
}
tet <- tet[,peptide$name]
meta_data <- cbind.data.frame(meta_data,tet)
#calculate clone frequencies
meta_data <- tcr.clones.rhapsody(meta_data,tcr[rownames(meta_data),])

#rename/combine genes+abs
expr_matrix <- annotate.genes(rhapsody)

#create seurat object and mappings
seurat.all <- create.seurat(expr_matrix,meta_data,experiment_id,
                            resolution = .6)

####################################

#### assign tetramer specificity

tetcounts <- as.sparse(t(tet))
tetcounts <- tetcounts+1
tet_assay <- CreateAssayObject(counts = tetcounts)
seurat.all$combined[["tet"]] <- tet_assay
seurat.all$combined <- NormalizeData(seurat.all$combined, assay = "tet", normalization.method = "LogNormalize")
#Demultiplexing tetramer
seurat.all$combined <- tetDemux(seurat.all$combined,assay="tet",positive.quantile = 0.99, nsamples = 300, dist_cutoff = 0.1, SNR_cutoff = 3)
#refined tet calls
seurat.all$combined <- tet.call.refined(seurat.obj = seurat.all$combined)
  
#antigen boxplots ~ good QC plots
tet.boxplots(seurat.all$combined,
             unique(seurat.all$combined@meta.data$tet.call.refined))

#save RDS
saveRDS(object = seurat.all,file = "seurat_all.rds")

#export mrna, abseq, metadata for totalvi
rna_cols <- rownames(expr_matrix)[!grepl("AbSeq$",rownames(expr_matrix))]
abseq_cols <- rownames(expr_matrix)[grepl("AbSeq$",rownames(expr_matrix))]
rna <- t(expr_matrix[rna_cols,])
abseq <- t(expr_matrix[abseq_cols,])
write.csv(rna,file = "mrna.csv",quote = FALSE)
write.csv(abseq,file = "abseq.csv",quote = FALSE)
write.csv(data.frame(seurat.all$combined@meta.data),file = "meta_data.csv",quote = FALSE)

Idents(seurat.all$combined) <- "seurat_clusters"

#########################################
### Some Prelim Analysis Plots #####

#clusters
DimPlot(seurat.all$combined,group.by = "seurat_clusters",pt.size = 1)
#sample id
DimPlot(seurat.all$combined, label = FALSE,group.by = 'donor.id',
        cols = colorRampPalette(brewer.pal(n=8,name="Set1"))(11))
DimPlot(seurat.all$combined, label = FALSE,group.by = 'time.point',
        cols = c("magenta","cyan","grey50"))
#combined
#FeaturePlot(seurat.all$combined,features = c(colnames(abseq)),
#            cols = c('blue','cyan','green','yellow','orange','red'))
#plot all abseq
FeaturePlot(seurat.all$combined,features = rownames(expr_matrix)[grepl("AbSeq",rownames(expr_matrix))],
            cols = c('blue','cyan','green','yellow','orange','red'))

DimPlot(seurat.all$combined,group.by = "clonality",cols = c("red","blue","grey80"))
FeaturePlot(seurat.all$combined,features = "f")
DimPlot(seurat.all$combined,group.by = "donor.id",shape.by = "time.point",split.by = "origin",ncol = 3)
DimPlot(seurat.all$combined,group.by = "time.point")

#END
