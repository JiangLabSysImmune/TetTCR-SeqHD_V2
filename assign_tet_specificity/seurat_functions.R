library(ggplot2)
library(reshape2)
library(Seurat)
library(sctransform)
library(cowplot)
library(ggrepel)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)

gather.meta <- function(md,sample,experiment_id){
  cells <- rownames(md)
  md$cell <- rownames(md)
  md <- merge(md,sample,by = 'Sample_Tag')
  rownames(md) <- md$cell
  md <- md[cells,]
  md$experiment <- experiment_id
  
  return(md)
}

origin <- function(seurat.obj){
  df <- seurat.obj@meta.data
  df$origin <- ifelse(grepl("_",df$tet_call),"CR",
                      ifelse(grepl("hCOV-S",df$tet_call),"hCOV",
                             ifelse(grepl("CMV",df$tet_call),"CMV",
                                    ifelse(grepl("EBV",df$tet_call),"EBV",
                                           ifelse(grepl("MART1",df$tet_call),"Human",
                                                  ifelse(grepl("IAV",df$tet_call),"IAV",
                                                         ifelse(grepl("Negative",df$tet_call),"Negative",
                                                                ifelse(grepl("Sticky",df$tet_call),"Sticky",
                                                                       ifelse(grepl("^COV-",df$tet_call),"SARS-2",
                                                                              ifelse(grepl("HUMAN",df$tet_call),"Human",
                                                                                     ifelse(grepl("HCV",df$tet_call),"HCV","Other")))))))))))
  
  return(df)
}

annotate.genes <- function(df){
  gene_annotation <- as.data.frame(colnames(df))
  gene_annotation <- as.data.frame(gene_annotation)
  colnames(gene_annotation) <- 'target'
  gene_annotation$target <- as.character(gene_annotation$target)
  gene_annotation$gene_short_name <- 'gene'
  for (i in seq(1,length(gene_annotation$target))){
    if (grepl('pAbO',gene_annotation$target[i])){
      gene <- gene_annotation$target[i]
      name <- paste(paste(unlist(strsplit(gene,split = '\\|'))[1],unlist(strsplit(gene,split = '\\|'))[2],sep='-'),"AbSeq",sep='-')
    }else{
      gene <- gene_annotation$target[i]
      name <- unlist(strsplit(gene,split = '\\|'))[1]
    }
    gene_annotation[i,'gene_short_name'] <- name
  }
  
  colnames(df) <- gene_annotation$gene_short_name
  df <- rowsum(t(df), group = colnames(df), na.rm = T)
  
  return(df)
}

tcr.clones <- function(tcra,tcrb,meta.data){
  select.top.chain <- function(tcr){
    tcr <- tcr[,c("X1st_cdr_aa","X1st_v","X1st_j","X2nd_cdr_aa","X2nd_v","X2nd_j")]
    tcr$cdr <- ifelse(!grepl("\\?",tcr$X1st_cdr_aa)&tcr$X1st_cdr_aa!="", tcr$X1st_cdr_aa,
                      ifelse(!grepl("\\?",tcr$X2nd_cdr_aa)&tcr$X2nd_cdr_aa!="",tcr$X2nd_cdr_aa,""))
    tcr$v <- ifelse(!grepl("\\?",tcr$X1st_cdr_aa)&tcr$X1st_cdr_aa!="", tcr$X1st_v,
                    ifelse(!grepl("\\?",tcr$X2nd_cdr_aa)&tcr$X2nd_cdr_aa!="",tcr$X2nd_v,""))
    tcr$j <- ifelse(!grepl("\\?",tcr$X1st_cdr_aa)&tcr$X1st_cdr_aa!="", tcr$X1st_j,
                    ifelse(!grepl("\\?",tcr$X2nd_cdr_aa)&tcr$X2nd_cdr_aa!="",tcr$X2nd_j,""))
    tcr <- tcr[,c("cdr","v","j")]
    
    return(tcr)
  }
  tcra <- select.top.chain(tcra)
  tcra <- tcra %>%
    dplyr::rename(
      alpha = cdr,
      trav = v,
      traj = j
    )
  tcrb <- select.top.chain(tcrb)
  tcrb <- tcrb %>%
    dplyr::rename(
      beta = cdr,
      trbv = v,
      trbj = j
    )
  tcr <- cbind.data.frame(tcra,tcrb)
  tcr <- tcr[rownames(meta.data),]
  meta.data <- cbind.data.frame(meta.data,tcr)
  meta.data$clone <- ifelse(meta.data$alpha!=""&meta.data$beta!="",
                            paste(paste(meta.data$beta,meta.data$alpha,sep='-'),meta.data$donor.id,sep='_'),"noPair")
  meta.data$cloneSize <- ifelse(meta.data$clone=='noPair',0,table(meta.data$clone)[meta.data$clone])
  meta.data$clonality <- ifelse(meta.data$cloneSize==0,'noPairing',
                                ifelse(meta.data$cloneSize==1,'not_clonal','clonal'))
  meta.data$beta.clone <- ifelse(meta.data$beta!='', paste0(meta.data$beta,'-',meta.data$trbv,'_',meta.data$donor.id),'noTCRb')
  meta.data$beta.cloneSize <- ifelse(meta.data$beta.clone=='noTCRb',0,table(meta.data$beta.clone)[meta.data$beta.clone])
  meta.data$beta.clonality <- ifelse(meta.data$beta.cloneSize==0,'noTCRb',
                                ifelse(meta.data$beta.cloneSize==1,'not_clonal','clonal'))
  
  return(meta.data)
}

tcr.clones.rhapsody <- function(meta_data,tcr){
  tcr$trav <- tcr$TCR_Alpha_Gamma_V_gene_Dominant
  tcr$traj <- tcr$TCR_Alpha_Gamma_J_gene_Dominant
  tcr$tcra <- tcr$TCR_Alpha_Gamma_CDR3_Translation_Dominant
  tcr$trbv <- tcr$TCR_Beta_Delta_V_gene_Dominant
  tcr$trbj <- tcr$TCR_Beta_Delta_J_gene_Dominant
  tcr$tcrb <- tcr$TCR_Beta_Delta_CDR3_Translation_Dominant
  tcr$tcr <- paste(tcr$TCR_Alpha_Gamma_CDR3_Translation_Dominant,tcr$TCR_Beta_Delta_CDR3_Translation_Dominant,sep = "_")
  tcr$tcr <- tcr %>%
    mutate(tcr = replace(tcr, tcr == "_",""),
           tcr = replace(tcr, startsWith(tcr,"_"),""),
           tcr = replace(tcr, endsWith(tcr,"_"),"")) %>%
    pull(tcr)
  tcr <- tcr[,c("trav","traj","tcra","trbv","trbj","tcrb","tcr")]
  md <- cbind.data.frame(meta_data,tcr)
  md <- md %>%
    group_by(sample.id,tcr) %>%
    mutate(clone.counts = n()) %>%
    group_by(sample.id) %>%
    mutate(clone.counts = replace(clone.counts, tcr == "",0), 
           f = (clone.counts / sum(clone.counts)),
           spike.freq = f * flow.freq * sort.freq)
  md$clonality <- ifelse(md$clone.counts == 0,"unpaired",
                      ifelse(md$clone.counts == 1,"singleton","clonal"))
  md <- as.data.frame(md)
  md$tcr <- ifelse(md$tcr != "",paste(md$tcr,md$donor.id,sep='-'),"")
  rownames(md) <- md$cell
  
  return(md)
}

tet.boxplots <- function(seurat.obj,tet){
  ### plots antigen and clone boxplots
  ### each cell is a point, x-axis = antigens on panel, y-axis = umi counts
  ### tet is list of antigens to add to plot (x-axis)
  
  df <- seurat.obj@meta.data
  shapes <- c(17,19)
  names(shapes) <- c("Pre","Post")
  # by clone
  pdf(file = "clone-tet_boxplots.pdf",width = 8,height = 6)
  for (clone in unique(df[df$clonality=="clonal","tcr"])){
    clone.df <- df[df$tcr==clone,]
    #get size of pre and post time points
    if (nrow(clone.df[clone.df$time.point=="Pre",])>0){
      size.pre <- unique(clone.df[clone.df$time.point=="Pre","clone.counts"])
    }else{size.pre <- 0}
    if (nrow(clone.df[clone.df$time.point=="Post",])>0){
      size.post <- unique(clone.df[clone.df$time.point=="Post","clone.counts"])
    }else{size.post <- 0}
    size <- paste(size.pre,size.post,sep="|")
    p <- df %>%
      filter(tcr == clone) %>%
      dplyr::select(any_of(tet),cell,time.point) %>%
      tidyr::gather(variable,value,-c(cell,time.point)) %>%
      ggplot(aes(x = reorder(as.character(variable),-value,FUN = median),y=value))+
      geom_boxplot(show.legend = FALSE,outlier.shape = NA)+
      geom_jitter(aes(color=cell,shape=time.point),show.legend = FALSE)+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90,hjust = 1),
            plot.title = element_text(hjust = .5))+
      labs(title = paste(clone,size,sep=": "),x="Antigen",y="Tetramer UMI")+
      scale_shape_manual(values = shapes)
    
    plot(p)
  }
  dev.off()
  
  # by specificity
  pdf(file = "antigen-tet_boxplots.pdf",width = 8,height = 6)
  antigens <- names(sort(table(seurat.all$combined@meta.data$tet.call.refined)[table(seurat.all$combined@meta.data$tet.call.refined)>1 & names(table(seurat.all$combined@meta.data$tet.call.refined))!="Sticky" & names(table(seurat.all$combined@meta.data$tet.call.refined))!="Negative"],decreasing = TRUE))
  for (antigen in antigens){
    size <- nrow(df[df$tet.call.refined==antigen,])
    #donor colors
    n <- length(unique(seurat.obj@meta.data$donor.id))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    donor_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
    names(donor_cols) <- unique(seurat.obj@meta.data$donor.id)
    
    p <- df %>%
      filter(tet.call.refined == antigen) %>%
      dplyr::select(any_of(tet),cell,time.point,donor.id) %>%
      tidyr::gather(variable,value,-c(cell,time.point,donor.id)) %>%
      
      ggplot(aes(x = reorder(variable,-value,FUN = median),y=value))+
      geom_boxplot(show.legend = FALSE,outlier.shape = NA)+
      geom_jitter(aes(color=donor.id,shape=time.point),alpha=.6)+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90,hjust = 1),
            plot.title = element_text(hjust = .5))+
      labs(title = paste(antigen,size,sep=": "),x="Antigen",y="Tetramer UMI")+
      scale_color_manual(values = donor_cols)+
      scale_shape_manual(values = shapes)
    
    plot(p)
  }
  dev.off()
}

create.seurat <- function(expr_matrix,md,project_name,resolution){
  #split on gene/abseq
  gene <- expr_matrix[!grepl('AbSeq',rownames(expr_matrix)),]
  abseq <- expr_matrix[grepl("AbSeq",rownames(expr_matrix)),]
  
  #create object
  seurat_sct_gene <- CreateSeuratObject(counts = gene, project = project_name, min.cells = 3, min.features = 0)
  seurat_sct_gene <- SCTransform(seurat_sct_gene, verbose = FALSE)
  seurat_sct_ab <- CreateSeuratObject(counts = abseq, project = project_name, min.cells = 3, min.features = 0)
  seurat_sct_ab <- SCTransform(seurat_sct_ab, verbose = FALSE)
  seurat_sct_comb <- CreateSeuratObject(counts = expr_matrix, project = project_name, min.cells = 3, min.features = 0)
  seurat_sct_comb <- SCTransform(seurat_sct_comb, verbose = FALSE)
  
  #metadata
  seurat_sct_gene@meta.data <- cbind(seurat_sct_gene@meta.data,md[rownames(seurat_sct_gene@meta.data),])
  seurat_sct_gene@meta.data <- cbind(seurat_sct_gene@meta.data,t(as.data.frame(seurat_sct_ab@assays$SCT@scale.data)))
  seurat_sct_comb@meta.data <- cbind(seurat_sct_comb@meta.data,md[rownames(seurat_sct_comb@meta.data),])
  
  ##### Dim Reduction
  #gene
  seurat_sct_gene <- RunPCA(seurat_sct_gene, verbose = FALSE)
  seurat_sct_gene <- RunUMAP(seurat_sct_gene, dims = 1:30, verbose = FALSE)
  seurat_sct_gene <- FindNeighbors(seurat_sct_gene, dims = 1:30, verbose = FALSE)
  seurat_sct_gene <- FindClusters(seurat_sct_gene, verbose = FALSE,resolution = resolution)
  
  #abseq
  #seurat_sct_ab <- RunPCA(seurat_sct_ab, verbose = FALSE)
  #seurat_sct_ab <- RunUMAP(seurat_sct_ab, dims = 1:30, verbose = FALSE)
  #seurat_sct_ab <- FindNeighbors(seurat_sct_ab, dims = 1:30, verbose = FALSE)
  #seurat_sct_ab <- FindClusters(seurat_sct_ab, verbose = FALSE,resolution = 1)
  
  #combined
  seurat_sct_comb <- RunPCA(seurat_sct_comb, verbose = FALSE)
  seurat_sct_comb <- RunUMAP(seurat_sct_comb, dims = 1:30, verbose = FALSE)
  seurat_sct_comb <- FindNeighbors(seurat_sct_comb, dims = 1:30, verbose = FALSE)
  seurat_sct_comb <- FindClusters(seurat_sct_comb, verbose = FALSE,resolution = resolution)

  return(list("gene"=seurat_sct_gene,"abseq"=seurat_sct_ab,"combined"=seurat_sct_comb))
}

create.seurat.scvi <- function(expr,rna.counts,reduction,meta_data,res){
  # use combined rna and prot expression from total-vi output as expr
  # use previously built RNA (gene only) seurat obj as suerat.rna (for rna counts)
  # use latent matrix as reduction (total-vi output)
  # use meta_data that is used in seurat v4 meta_data pipeline
  rownames(reduction) <- rownames(expr)
  seurat.scvi <- CreateSeuratObject(t(expr),min.cells = 0,min.features = 0)
  seurat.scvi@meta.data <- cbind.data.frame(seurat.scvi@meta.data,
                                            meta_data[rownames(seurat.scvi@meta.data),])
  seurat.scvi <- FindVariableFeatures(seurat.scvi,selection.method = 'vst',
                                      nfeatures = ncol(expr))
  seurat.scvi[['scvi']] <- CreateDimReducObject(embeddings = as.matrix(reduction),
                                                key="scvi_",
                                                assay = DefaultAssay(seurat.scvi))
  seurat.scvi[['mrna']] <- CreateAssayObject(counts = t(rna.counts))
  seurat.scvi <- FindNeighbors(seurat.scvi,dims = 1:ncol(reduction),reduction = 'scvi')
  seurat.scvi <- FindClusters(seurat.scvi,resolution = res,algorithm = 2)
  seurat.scvi <- RunUMAP(seurat.scvi,dims = 1:ncol(reduction),reduction = 'scvi',
                         n.components = 2,min.dist = .2)
  #prepare assays
  seurat.scvi@assays$SCVI <- seurat.scvi@assays$RNA
  Key(seurat.scvi@assays$SCVI) <- "SCVI_"
#  seurat.scvi@assays$mRNA <- seurat.rna@assays$RNA
#  Key(seurat.scvi@assays$mRNA) <- "mRNA_"
#  seurat.scvi@assays$AbSeq <- seurat.scvi[colnames(expr)[grepl("AbSeq",colnames(expr))]]@assays$RNA
#  Key(seurat.scvi@assays$AbSeq) <- "AbSeq_"
  
  return(seurat.scvi)
}

plot.main <- function(df,plot_title){
  clust <- DimPlot(df, label = TRUE,pt.size = .5) + NoLegend()
  donor <- DimPlot(df, label = FALSE,group.by = 'donor.id',
                   cols = colorRampPalette(brewer.pal(n=8,name="Set1"))(11))
  condition <- DimPlot(df, label = FALSE,group.by = 'condition',
                       cols = c("magenta","cyan","grey30"))
  cell_type <- DimPlot(df,label = FALSE,group.by = "cell.label") +
    scale_color_manual(values = c("darkorange","steelblue","grey50"))
  p <- plot_grid(clust,donor,condition,cell_type)
  pdf(paste0(plot_title,'.pdf'),width=12,height = 8)
  plot(p)
  dev.off()
  return(p)
}

plot.main.tcrb <- function(df,plot_title){
  clust <- DimPlot(df, label = TRUE,pt.size = .5) + NoLegend()
  donor <- DimPlot(df, label = FALSE,group.by = 'donor.id',
                   cols = colorRampPalette(brewer.pal(n=8,name="Set1"))(11))
  condition <- DimPlot(df, label = FALSE,group.by = 'condition',
                       cols = c("magenta","cyan","grey30"))
  tcrb <- DimPlot(df,group.by = 'beta.clonality',
                  cols = c("skyblue","grey80","red"),order = "clonal")
  p <- plot_grid(clust,donor,condition,tcrb)
  pdf(paste0(plot_title,'.pdf'),width=12,height = 8)
  plot(p)
  dev.off()
  return(p)
}

plot.abseq <- function(df,plot_title){
  #only works for "gene" slot in seurat obj
  abseqs <- colnames(df@meta.data)[grepl('AbSeq',colnames(df@meta.data))]
  pdf(paste0(plot_title,'.pdf'),width = 9,height = 6)
  for (i in seq(1,length(abseqs),4)){
    plot(FeaturePlot(df,c(abseqs[i],abseqs[i+1],abseqs[i+2],abseqs[i+3]),
                     cols = c('blue','cyan','green','yellow','orange','red')))
  }
  dev.off()
}

plot.sig.markers <- function(df,assay,plot_title){
  #assay = type of feature (i.e. "RNA"/"SCT" or "AbSeq")
  if (assay == "mRNA"){
    markers.all <- FindAllMarkers(df,assay = assay,test.use = "negbinom")
    markers.all <- markers.all[markers.all$gene%in%df@assays$SCVI@var.features,]
  } else {
    markers.all <- FindAllMarkers(df,assay = assay)
  }
  markers.all.sig <- markers.all[markers.all$p_val_adj<=.01,]
  sig_genes <- unique(markers.all.sig$gene)
  if (assay=="RNA"){
    sig_genes <- paste0('rna_',sig_genes)
  } else if(assay=='SCT'){
    sig_genes <- paste0("sct_",sig_genes)
  } else if(assay=='SCVI'|assay=="mRNA"){
    sig_genes <- paste0("SCVI_",sig_genes)
  }
  pdf(paste0(plot_title,'.pdf'),width = 9,height = 6)
  for (i in seq(1,length(sig_genes),4)){
    plot(FeaturePlot(df, c(sig_genes[i],sig_genes[i+1],sig_genes[i+2],sig_genes[i+3]),
                     cols = c('blue','cyan','green','yellow','orange','red')))
  }
  dev.off()
  
  return(markers.all.sig)
}

label.cutoff <- function(MIDs, breaks = 50){
  cat = colnames(MIDs)
  summary <- MIDs
  #determine negative threshold for each category with normal mixture model
  neg <- c()
  cat.name <- list()
  pdf(file="./QC_negative_threshold.pdf",width = 4,height=4)
  for (j in seq(1:length(cat))){
    cat.name[[j]] <- as.character(colnames(MIDs)[j])
    hist <- hist(asinh(summary[,cat[j]]),breaks = breaks,plot=FALSE)
    df <- data.frame(mid=hist$mids,cou=hist$counts)
    #iterate mixture model
    for(i in seq(0.1,5,by=0.1)){
      res <- try(mix.model <- mix(as.mixdata(df),mixparam(mu=c(0,2),sigma=i)))
      if(inherits(res, "try-error")){
        next
      }else{
        break
      }
    }
    #negative threshold
    f <- function(x){dnorm(x,res$parameters$mu[1],res$parameters$sigma[1]) * res$parameters$pi[1] - dnorm(x,res$parameters$mu[2],res$parameters$sigma[2]) * res$parameters$pi[2]}
    thr <- uniroot(f,interval = c(res$parameters$mu[1],res$parameters$mu[2]))$root
    #threshold plot
    plot(res,main = cat[j],xlab="Normalized SCT")
    abline(v=thr, col=c("blue"), lty=2,lwd=2)
    neg <- c(neg,thr)
  }
  names(neg) <- cat
  dev.off()
  
  #list of cells for each marker
  labels <- c()
  for (cell in rownames(MIDs)){
    cell.label <- c()
    for (marker in markers){
      if (MIDs[cell,marker]>neg[marker]){
        cell.label <- c(cell.label,"pos")
      } else {cell.label <- c(cell.label,"neg")}
    } 
    tab <- table(cell.label)
    if (TRUE %in% grepl("pos",names(tab)) & TRUE %in% tab["pos"]>0){
      ind <- which(cell.label %in% "pos")
      label <- paste(markers[ind],collapse="|")
    } else{
      label = "Undetermined"
    }
    labels <- c(labels,label)
  }
  MIDs$label <- labels
  
  return(list("cutoff"=neg,"labels"=MIDs))
}

plot.cluster.hm <- function(seurat.obj,markers,type){
  # plots cluster heatmap
  # markers = df output of markers from FindAllMarkers (can be rna + abseq combined)
  # type = type of analysis (choose "seurat" or "scvi")
  #   if type is "seurat", SCT and AbSeq assays will be called
  #   if type is "scvi", 
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  if (type == "seurat"){
    rna <- as.data.frame(seurat.obj@assays$SCT@data)
    abseq <- as.data.frame(seurat.obj@assays$AbSeq@data)
    expr <- rbind.data.frame(rna,abseq)
  } else if (type == "scvi"){
    expr <- as.data.frame(seurat.obj@assays$SCVI@data)
  } else{stop("Error! Invalid argument. Please use 'seurat' or 'scvi' for type argument.")}
  meta.data <- seurat.obj@meta.data
  
  genes <- c()
  for (cluster in unique(markers$cluster)){
    clust.df <- markers[markers$cluster==cluster,]
    up <- as.character(head(clust.df[order(-clust.df$avg_log2FC),]$gene,n = 10))
    dn <- as.character(tail(clust.df[order(-clust.df$avg_log2FC),]$gene,n = 10))
    genes <- c(genes,up,dn)
  }
  genes <- unique(genes)
  
  subset <- as.data.frame(t(expr[genes,]))
  #subset$cluster <- paste('cluster',meta.data$seurat_clusters,sep='_')
  subset$cluster <- meta.data$seurat_clusters
  d = plyr::ddply(subset,'cluster', function(x) colMeans(x[genes]))
  rownames(d) <- d[,1]
  d <- d[,-1]

  n = length(rownames(d))
  cols = gg_color_hue(n)
  
  # Data frame with column annotations.
  mat_col <- data.frame(cluster = rownames(d))
  rownames(mat_col) <- rownames(d)
  
  # List with colors for each annotation.
  mat_colors <- list(cluster = cols)
  names(mat_colors$cluster) <- rownames(d)
  
  hm <- pheatmap(mat = as.matrix(t(d)),scale='row',color = colorRampPalette(c('blue','white','red'))(100),
                 annotation_col = mat_col,annotation_colors = mat_colors,fontsize_row = 8,
                 border_color = NA)
  
  return(hm)
}

cluster.distribution <- function(seurat.obj,groups){
  cluster <- seurat.obj$seurat_clusters
  donor <- seurat.obj$donor.id
  condition <- seurat.obj$condition
  df <- cbind.data.frame(cluster,donor,condition)
  df <- df %>% 
    group_by(donor,cluster,condition) %>% 
    summarise(N=n()) %>%
    group_by(donor) %>%
    mutate(Freq = N / sum(N))
  for (donor in unique(df$donor)){
    all.clust <- table(df[df$donor==donor,"cluster"])
    zero.clust <- names(all.clust[which(all.clust %in% 0)])
    if (length(zero.clust)>0){
      cond <- unique(df[df$donor==donor,] %>% pull(condition))
      for (clust in zero.clust){
        row.add <- c(donor,clust,cond,0,0)
        df <- rbind.data.frame(df,row.add)
      }
    }
  }
  df$Freq <- as.numeric(as.character(df$Freq))
  comparisons <- compare_means(Freq ~ condition,data=df,
                               group.by = "cluster")
  df <- df[complete.cases(df),]
  p <- ggplot(df,aes(x=cluster,y=Freq,color=condition))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_dodge(width=.75),size=1)+
    theme_bw()+
    scale_color_manual(values = c("magenta","cyan"),na.value='grey30')+
    stat_compare_means(label="p.signif",method = "wilcox")
  
  return(p)
}

gather.marker.data <- function(marker.df,markers.select,p_thresh=.01,fc_thresh=.25,n_genes=6){
  marker.df$log2.padj <- -log2(marker.df$p_val_adj)
  marker.df$gene <- rownames(marker.df)
  sig.markers <- rownames(marker.df[marker.df$p_val_adj<=p_thresh & abs(marker.df$avg_log2FC)>=fc_thresh,])
  marker.df$signif <- ifelse(marker.df$gene %in% sig.markers,"signif",NA)
  marker.df <- marker.df %>% arrange(marker.df,p_val_adj)
  pos.markers <- head(marker.df[marker.df$signif=="signif" & marker.df$avg_log2FC>0,"gene"],n=n_genes)
  neg.markers <- head(marker.df[marker.df$signif=="signif" & marker.df$avg_log2FC<0,"gene"],n=n_genes)
  lab.markers <- c(pos.markers,neg.markers,markers.select)
  marker.df$labs <- ifelse(marker.df$gene %in% lab.markers,marker.df$gene,NA)
  
  return(marker.df)
}

plot.diff.cluster <- function(seurat.obj,cluster,group,conditions,markers.select,type){
  # if type = "seurat", will use combined sct as assay (maybe update for wnn??)
  # if type = "scvi", will use rna counts with negbinom and abseq scvi with wilcox
  if (type == "seurat"){
    markers <- FindMarkers(seurat.obj,
                           group.by = group,subset.ident = cluster,
                           ident.1 = conditions[1],ident.2 = conditions[2],
                           logfc.threshold = 0,test.use = "negbinom")
    markers <- gather.marker.data(markers,markers.select = markers.select)
    p <- ggplot(markers,aes(x=avg_log2FC,y=log2.padj,color=signif,label=labs))+
      geom_point(alpha=.5)+
      theme_bw()+labs(subtitle=paste(conditions[1],conditions[2],sep = ' vs. '),
                      title = paste0("Cluster: ",cluster),
                      x="Avg. Log2FC",y="-log2(padj)")+
      geom_text_repel(segment.color = "grey50",segment.alpha = .6,size=2)+
      theme(legend.position = "none")+
      scale_color_manual(values = c("red"),na.valu="grey50")+
      geom_vline(xintercept = 0,linetype="dashed",alpha=.3,color="steelblue")
    
  }else if(type == "scvi"){
    markers.rna <- FindMarkers(seurat.obj,
                               group.by = group,subset.ident = cluster,
                               ident.1 = conditions[1],ident.2 = conditions[2],
                               logfc.threshold = 0,test.use = "negbinom",
                               assay = "mRNA")
    markers.abseq <- FindMarkers(seurat.obj,
                                 group.by = group,subset.ident = cluster,
                                 ident.1 = conditions[1],ident.2 = conditions[2],
                                 logfc.threshold = 0,
                                 assay = "AbSeq")
    markers.rna <- gather.marker.data(markers.rna,markers.select = markers.select)
    markers.abseq <- gather.marker.data(markers.abseq,markers.select = markers.select)
    p.rna <- ggplot(markers.rna,aes(x=avg_log2FC,y=log2.padj,color=signif,label=labs))+
      geom_point(alpha=.5)+
      theme_bw()+labs(subtitle=paste(conditions[1],conditions[2],sep = ' vs. '),
                      title = paste0("Cluster: ",cluster),
                      x="Avg. Log2FC",y="-log2(padj)")+
      geom_text_repel(segment.color = "grey50",segment.alpha = .6,size=2)+
      theme(legend.position = "none")+
      scale_color_manual(values = c("red"),na.valu="grey50")+
      geom_vline(xintercept = 0,linetype="dashed",alpha=.3,color="steelblue")
    p.abseq <- ggplot(markers.abseq,aes(x=avg_log2FC,y=log2.padj,color=signif,label=labs))+
      geom_point(alpha=.5)+
      theme_bw()+labs(subtitle=paste(conditions[1],conditions[2],sep = ' vs. '),
                      title = paste0("Cluster: ",cluster),
                      x="Avg. Log2FC",y="-log2(padj)")+
      geom_text_repel(segment.color = "grey50",segment.alpha = .6,size=2)+
      theme(legend.position = "none")+
      scale_color_manual(values = c("red"),na.valu="grey50")+
      geom_vline(xintercept = 0,linetype="dashed",alpha=.3,color="steelblue")
    p <- list("rna.plot"=p.rna,"abseq.plot"=p.abseq)
  } else {stop("Error! Invalid argument. Please use 'seurat' or 'scvi' for type argument.")}
  
  return(p)
}

plot.diff.subsets <- function(seurat.obj,clusters,group,conditions,title,markers.select,type){
  seurat.subset <- subset(seurat.obj,subset=seurat_clusters%in%clusters)
  if (type=="seurat"){
    markers <- FindMarkers(seurat.subset,
                           group.by=group,
                           ident.1=conditions[1],ident.2=conditions[2],
                           logfc.threshold = 0,test.use="negbinom")
    markers <- gather.marker.data(markers,markers.select)
    p <- ggplot(markers,aes(x=avg_log2FC,y=log2.padj,color=signif,label=labs))+
      geom_point(alpha=.5)+
      theme_bw()+labs(subtitle=paste(conditions[1],conditions[2],sep = ' vs. '),
                      title = title,
                      x="Average LogFC",y="-log2(padj)")+
      geom_text_repel(segment.color = "grey50",segment.alpha = .6,size=2)+
      theme(legend.position = "none")+
      scale_color_manual(values = c("red"),na.valu="grey50")+
      geom_vline(xintercept = 0,linetype="dashed",alpha=.3,color="steelblue")
  } else if(type=="scvi"){
    markers.rna <- FindMarkers(seurat.subset,
                               group.by = group,
                               ident.1 = conditions[1],ident.2 = conditions[2],
                               logfc.threshold = 0,test.use = "negbinom",
                               assay = "mRNA")
    markers.abseq <- FindMarkers(seurat.subset,
                                 group.by = group,
                                 ident.1 = conditions[1],ident.2 = conditions[2],
                                 logfc.threshold = 0,
                                 assay = "AbSeq")
    markers.rna <- gather.marker.data(markers.rna,markers.select = markers.select)
    markers.abseq <- gather.marker.data(markers.abseq,markers.select = markers.select)
    p.rna <- ggplot(markers.rna,aes(x=avg_log2FC,y=log2.padj,color=signif,label=labs))+
      geom_point(alpha=.5)+
      theme_bw()+labs(subtitle=paste(conditions[1],conditions[2],sep = ' vs. '),
                      title = title,
                      x="Avg. Log2FC",y="-log2(padj)")+
      geom_text_repel(segment.color = "grey50",segment.alpha = .6,size=2)+
      theme(legend.position = "none")+
      scale_color_manual(values = c("red"),na.valu="grey50")+
      geom_vline(xintercept = 0,linetype="dashed",alpha=.3,color="steelblue")
    p.abseq <- ggplot(markers.abseq,aes(x=avg_log2FC,y=log2.padj,color=signif,label=labs))+
      geom_point(alpha=.5)+
      theme_bw()+labs(subtitle=paste(conditions[1],conditions[2],sep = ' vs. '),
                      title = title,
                      x="Avg. Log2FC",y="-log2(padj)")+
      geom_text_repel(segment.color = "grey50",segment.alpha = .6,size=2)+
      theme(legend.position = "none")+
      scale_color_manual(values = c("red"),na.valu="grey50")+
      geom_vline(xintercept = 0,linetype="dashed",alpha=.3,color="steelblue")+
      xlim(-50,50)
    p <- list("rna.plot"=p.rna,"abseq.plot"=p.abseq)
  }
  
  return(p)
}

plot.vln.condition <- function(seurat.obj,markers){
  vln.plots <- list()
  for (marker in markers){
    vln.plots[[marker]] <- VlnPlot(subset(seurat.obj,subset = condition!="NA"),
                                 group.by = "donor.id",
                                 features = marker,
                                 pt.size = 0,cols=c("magenta","cyan"))+
      aes(fill=subset(seurat.obj,subset = condition!="NA")$condition)+
      scale_color_manual(values = c("magenta","cyan"))+NoLegend()
  }
  p <- plot_grid(plotlist = vln.plots,ncol = 4)
  
  return(p)
}

run.wnn <- function(seurat.obj,res){
  # need to have run SCT on mRNA stored in `seurat.obj$gene`
  # need to have loaded AbSeq counts into `seurat.obj$abseq`
  # res = clustering resolution
  seurat <- seurat.obj$gene
  seurat@assays$AbSeq <- seurat.obj$abseq@assays$RNA
  DefaultAssay(seurat) <- "AbSeq"
  VariableFeatures(seurat) <- rownames(seurat[["AbSeq"]])
  seurat <- NormalizeData(seurat,normalization.method = "CLR",margin=2) %>%
    ScaleData() %>% RunPCA(reduction.name = "abpca")
  seurat <- FindMultiModalNeighbors(
    seurat, reduction.list = list("pca","abpca"),
    dims.list = list(1:30,1:18), modality.weight.name = "RNA.weight"
  )
  seurat <- RunUMAP(seurat, nn.name = "weighted.nn",reduction.name = "wnn.umap",
                 reduction.key = "wnnUMAP_")
  seurat <- FindClusters(seurat, graph.name = "wsnn", algorithm = 3, resolution = res,
                      verbose = FALSE)
  seurat.obj$wnn <- seurat
  
  return(seurat.obj)
}

nse <- function(p){
  n <- length(p)
  X <- c()
  for (i in seq(1:length(p))){
    x <- (p[i] * log(p[i])) / log(n)
    X <- c(X,x)
  }
  X <- -1 * (sum(X))
  
  return(X)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

bhatta_dist_clusters <- function(md){
  t <- md %>%
    filter(tcrb != "",
           seurat_clusters %in% 0:9) %>%
    mutate(beta_clone = paste(tcrb,donor.id,sep = "_")) %>%
    select(seurat_clusters,donor.id,beta_clone) %>%
    group_by(seurat_clusters,donor.id)  %>%
    mutate(tot_cells_cluster = n()) %>%
    group_by(seurat_clusters,donor.id,beta_clone) %>%
    summarize(cluster_freq = n() / tot_cells_cluster) %>%
    distinct()
  
  pairs <- apply(combn(c(0:(length(unique(t$seurat_clusters))-1)),2),2,paste,collapse='_')
  all.donors <- unique(t$donor.id)
  battacharya.df <- data.frame(all.donors)
  for (pair in pairs){
    
    p.id <- strsplit(pair,"_")[[1]][1]
    q.id <- strsplit(pair,"_")[[1]][2]
    
    df <- t %>%
      filter(seurat_clusters %in% c(p.id,q.id))
    donor_b_dists <- c()
    for (donor in all.donors){
      
      overlap <- intersect(df %>% filter(seurat_clusters==p.id,donor.id == donor) %>% pull(beta_clone), 
                           df %>% filter(seurat_clusters==q.id,donor.id == donor) %>% pull(beta_clone))
      if (length(overlap)>0){
        b_dists <- c()
        for (clone in overlap){
          p <- df %>% filter(beta_clone == clone,seurat_clusters == p.id) %>% pull(cluster_freq)
          q <- df %>% filter(beta_clone == clone,seurat_clusters == q.id) %>% pull(cluster_freq)
          b <- sqrt(p * q)
          b_dists <- c(b_dists,b)
        }
        donor_b_dists <- c(donor_b_dists,-log(sum(b_dists)))
      }else{
        donor_b_dists <- c(donor_b_dists,NA)
      }
    }
    battacharya.df[,pair] <- donor_b_dists
  }
  #plot heatmap
  x <- as.data.frame(melt(battacharya.df) %>%
                       na.omit() %>%
                       group_by(variable) %>%
                       summarize(mean = mean(value)) %>%
                       separate(variable,into = c("id_1","id_2"),sep = "_") %>%
                       spread(id_1,mean))
  rownames(x) <- x$id_2
  x <- as.matrix(subset(x, select = - id_2))
  p <- pheatmap::pheatmap(x,cluster_rows = F,cluster_cols = F,na_col = "white",display_numbers = TRUE,
                          color = colorRampPalette(c('red',"orange",'yellow',"lightyellow"))(20))
  
  return(list(df=battacharya.df,plot=p))
}

#END
