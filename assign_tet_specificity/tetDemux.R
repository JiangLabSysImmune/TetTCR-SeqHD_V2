#' @importFrom cluster clara
#' @importFrom Matrix colSums
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pnbinom kmeans
library(cluster)
library(Matrix)
library(fitdistrplus)
library(stats)
library(cgwtools)
library(tibble)
library(gt)
library(ggridges)

tetDemux <- function (object, assay = "tet", positive.quantile = 0.99, 
          init = NULL, nstarts = 100, kfunc = "clara", nsamples = 500, dist_cutoff = 0.1, SNR_cutoff = 3,
          seed = 42, verbose = TRUE) 
{
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(object = object, assay = assay, slot = "counts")[, 
                                                                          colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- init %||% (nrow(x = data) + 1)
  switch(EXPR = kfunc, kmeans = {
    init.clusters <- kmeans(x = t(x = GetAssayData(object = object, 
                                                   assay = assay)), centers = ncenters, nstart = nstarts)
    Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
  }, clara = {
    init.clusters <- clara(x = t(x = GetAssayData(object = object, 
                                                  assay = assay)), k = ncenters, samples = nsamples)
    Idents(object = object, cells = names(x = init.clusters$clustering), 
           drop = TRUE) <- init.clusters$clustering
  }, stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'"))
  average.expression <- AverageExpression(object = object, 
                                          assays = assay, verbose = FALSE)[[assay]]
  if (sum(average.expression == 0) > 0) {
    stop("Cells with zero counts exist as a cluster.")
  }
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    
    
    lst <- sort(average.expression[iter, ], index.return=TRUE, decreasing=FALSE)
    idx <- lapply(lst, `[`, lst$x %in% head(unique(lst$x),round(ncenters*dist_cutoff, 0)))$ix
    
    
    
    values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[idx])]
    
    
    
    
    fit <- suppressWarnings(expr = fitdist(data = values.use, 
                                           distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    #if no cells are over cutoff
    if (length(which(x = values > cutoff)) > 0){
      discrete[iter, names(x = which(x = values > cutoff))] <- 1
    }
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", 
                     cutoff, " reads"))
    }
  }
  
  cell_ID <- list()
  SNR <- list()
  for (cell in colnames(discrete)) {
    pos_counts <- mean(counts[discrete[, cell] == 1, cell])
    
    #Considering only top 3 of non-called tetramers
    neg_counts <- mean(counts[discrete[, cell] == 0, cell][order(counts[discrete[, cell] == 0, cell], decreasing = TRUE)][1:3])
    sigtonoise <- pos_counts/neg_counts
    
    cell_ID <- append(cell_ID, cell)
    SNR <- append(SNR, sigtonoise)
    
    #set cell as sticky if signal to noise less than 3x
    
    if ((sigtonoise < SNR_cutoff)&(sum(discrete[, cell]) >= 1)&(sum(discrete[, cell]) < 4)) {
      discrete[, cell] <- 1
    }
    
  }
  ID_SNR <- do.call(rbind, Map(data.frame, cellID=cell_ID, SNR=SNR))
  ID_SNR <- ID_SNR %>% remove_rownames %>% column_to_rownames(var="cellID")
  object <- AddMetaData(object = object, metadata = ID_SNR, col.name = "SNR")
  
  
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive == 2] <- "Doublet"
  classification.global[npositive == 3] <- "Triplet"
  classification.global[npositive > 3] <- "Sticky"
  
  
  donor.id = rownames(x = data)
  hash.max <- apply(X = counts, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = counts, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = counts, MARGIN = 2, FUN = maxn, nth = 2)
  hash.third <- apply(X = counts, MARGIN = 2, FUN = maxn, nth = 3)
  
  hash.maxID <- as.character(x = donor.id[sapply(X = 1:ncol(x = counts), 
                                                 FUN = function(x) {
                                                   return(which(x = counts[, x] == hash.max[x])[1])
                                                 })])
  hash.secondID <- as.character(x = donor.id[sapply(X = 1:ncol(x = counts), 
                                                    FUN = function(x) {
                                                      return(which(x = counts[, x] == hash.second[x])[1])
                                                    })])
  hash.thirdID <- as.character(x = donor.id[sapply(X = 1:ncol(x = counts), 
                                                    FUN = function(x) {
                                                      return(which(x = counts[, x] == hash.third[x])[1])
                                                    })])
  
  
  
  hash.margin <- hash.max - hash.second
  
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                 collapse = "_"))
  })
  
  triplet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x], hash.thirdID[x])), 
                 collapse = "_"))
  })
  
  
  
  
  
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == 
                                                                           "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == 
                                                                           "Doublet")]
  classification[classification.global == "Triplet"] <- triplet_id[which(x = classification.global == 
                                                                           "Triplet")]
  
  
  
  classification.metadata <- data.frame(hash.maxID, hash.secondID, hash.thirdID,
                                        hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, c("maxID", 
                                                          "secondID", "thirdID","margin", "classification", 
                                                          "classification.global"), sep = "_")
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, "_classification")
  
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, 
                                                            "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  
  triplets <- rownames(x = object[[]])[which(object[[paste0(assay, 
                                                            "_classification.global")]] == "Triplet")]
  Idents(object = object, cells = triplets) <- "Triplet"
  
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, 
                                                            "_classification.global")]] == "Sticky")]
  Idents(object = object, cells = doublets) <- "Sticky"
  object$hash.ID <- Idents(object = object)
  return(object)
}

tet.call.refined <- function(seurat.obj,SNR_cutoff = 3){
  ### refined tetramer classification 
  ### takes seurat object that has been run with tetDemux.R, outputs a seurat obj with a new tet.call.refined coulumn
  ### uses `tet_classification` column from tetDemux.R and attempts to rescue false "Sticky" or Multiplets classifications
  ### sets pre-defined cutoff for SNR (tet_umi / (total_umi - tet_umi)), rescues cells with SNR > cutoff
  ### Iterates over doublets/triplets and ensures 2nd/3rd ranked tet pass SNR cutoff to be labeled as cross-reactive
  
  #Set up arrays
  tet_calls <- c()
  SNR_list <- c()
  max_list <- c()
  md <- seurat.obj@meta.data
  
  #iterate over cells
  for (cell in rownames(md)){
    tet_class <- md[cell,"tet_classification.global"]
    tet_maxID <- md[cell,"tet_maxID"]
    tet_max_counts <- md[cell,tet_maxID]
    total_counts <- md[cell,'nCount_tet'] - md[cell,"nFeature_tet"]
    SNR <- tet_max_counts / (total_counts - tet_max_counts)
    #append to arrays
    SNR_list <- c(SNR_list,SNR)
    max_list <- c(max_list,tet_max_counts)
    
    if (tet_class == "Sticky" | tet_class == "Doublet" | tet_class == "Triplet"){
      if (tet_class == "Sticky" & SNR > SNR_cutoff){
        #pass, call is maxID
        tet_call <- tet_maxID
      }else if (tet_class == "Doublet" | tet_class == "Triplet"){
        tet_secondID <- md[cell,"tet_secondID"]
        tet_thirdID <- md[cell,"tet_thirdID"]
        tet_second_counts <- md[cell,tet_secondID]
        tet_third_counts <- md[cell,tet_thirdID]
        SNR_second <- tet_second_counts / (total_counts - tet_max_counts - tet_second_counts)
        SNR_third <- tet_third_counts / (total_counts - tet_max_counts - tet_second_counts - tet_third_counts)
        if (tet_class == "Doublet" & SNR > SNR_cutoff & SNR_second > SNR_cutoff){
          #pass, CR
          tet_call <- md[cell,"tet_classification"]
        }else if(tet_class == "Triplet" & SNR > SNR_cutoff & SNR_second > SNR_cutoff & SNR_third > SNR_cutoff){
          #pass, CR
          tet_call <- md[cell,"tet_classification"]
        }else if(SNR > SNR_cutoff){
          #pass
          tet_call <- tet_maxID
        }else{
          #fail
          tet_call <- "Sticky"
        }
      }else{
        tet_call <- md[cell,"tet_classification"]
      }
    }else{
      tet_call <- md[cell,"tet_classification"]
    }
    tet_calls <- c(tet_calls,tet_call)
  }
  seurat.obj@meta.data$tet.call.refined <- tet_calls
  seurat.obj@meta.data$SNR.refined <- SNR_list
  seurat.obj@meta.data$tet.max.counts <- max_list
  
  return(seurat.obj)
}

#END
