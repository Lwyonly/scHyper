scHyper_function <- function(fpath){
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  suppressMessages(library(lobstr))
  set.seed(123)

  result <- read.csv(fpath, header = TRUE)
  result <- separate(result, 'ligand', c('Sender', 'Ligand'), sep = '\\.')
  result <- separate(result, 'receptor', c('Receiver', 'Receptor'), sep = '\\.')
  result <- result[which(result$Sender != result$Receiver),]
  result$Ligand <- gsub('_', '&', result$Ligand)
  result$Receptor <- gsub('_', '&', result$Receptor)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  return(result)
}

CellChat_function <- function(data, meta, species = 'human'){
  suppressMessages(library(Seurat))
  suppressMessages(library(CellChat))
  suppressMessages(library(tidyverse))
  suppressMessages(library(lobstr))
  set.seed(123)

  data <- as.matrix(data)
  meta$celltype <- as.character(meta$celltype)
  
  cellchat <- createCellChat(object = data, meta = meta, group.by = "celltype")
  if(species == 'human'){
    cellchat@DB <- CellChatDB.human
  }else{
    cellchat@DB <- CellChatDB.mouse
  }
  
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  result <- computeCommunProb(cellchat) %>% filterCommunication(.) %>% 
    subsetCommunication(.)
  
  result <- result[,1:5]
  colnames(result) <- c('Sender','Receiver','Ligand', 'Receptor', 'LRscore')
  result <- result[which(result$Sender != result$Receiver),]
  result$Ligand <- gsub('_', '&', result$Ligand)
  result$Receptor <- gsub('_', '&', result$Receptor)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  return(result)
}

iTALK_function <- function(data, meta, top_genes){
  suppressMessages(library(iTALK))
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
 
  data <- as.matrix(data)
  data <- as.data.frame(t(data))
  data <- cbind(data, cell_type = meta[rownames(data), "celltype"])

  # find top 50 percent highly expressed genes
  highly_exprs_genes <- rawParse(data, top_genes=top_genes, stats='mean')#top_genes=50
  # find the ligand-receptor pairs from highly expressed genes
  comm.list<-c('growth factor','other','cytokine','checkpoint')
  
  result <- NULL
  for(comm.type in comm.list){
    res.tmp <- FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm.type)
    res.tmp <- res.tmp[order(res.tmp$cell_from_mean_exprs*res.tmp$cell_to_mean_exprs,decreasing=T),]
    result<-rbind(result,res.tmp)
  }
  
 
  result$LRscore <- result$cell_from_mean_exprs*result$cell_to_mean_exprs
  result[,c(3,5,7)] <- NULL
  colnames(result) <- c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore')
  result <- result[which(result$Receiver != result$Sender), ]
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  return(result)
}

SCSR_function <- function(data, meta, s.score, species = 'human'){
  suppressMessages(library(SingleCellSignalR))
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyr))
  suppressMessages(library(dplyr))
  suppressMessages(library(lobstr))
  set.seed(123)
 
  if(species == 'human'){
    species <- "homo sapiens"
  }else{
    species <- "mus musculus"
  }
  
  data <- as.matrix(data)
  meta <- as.data.frame(meta)
  #meta$celltype <- factor(meta$celltype)
  #meta$celltype <- as.integer(meta$celltype)
  
  # Digitize Labels
  i <- 1
  for(ct in unique(meta$celltype)){
    meta[which(meta$celltype == ct), "ct_num"] <- i
    i <- i+1
  }
  c.names <- as.character(unique(meta$celltype))
  
  #clust.ana <- cluster_analysis(data = data, genes = rownames(data), cluster = meta$ct_num, c.names = c.names, dif.exp=TRUE)
  signal <- cell_signaling(data = data, genes = rownames(data), int.type = "paracrine",
                           species = species, cluster = meta$ct_num, c.names = c.names, s.score = s.score, write = FALSE)
  inter.net <- inter_network(data = data, signal = signal, genes = rownames(data), 
                             cluster = meta$ct_num, c.names = c.names, species = species,write = FALSE)
  
  result <- inter.net$`full-network`
  result$interaction.type <- NULL
  result <- separate(result, 'ligand', c('Sender', 'Ligand'), sep = '\\.')
  result <- separate(result, 'receptor', c('Receiver','Receptor'), sep = '\\.')
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- distinct(result, all, .keep_all = TRUE)
  rownames(result) <- NULL
  return(result)
}

CellPhoneDB2_function <- function(fpath){
  suppressMessages(library(tidyr))
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  file.path <- paste(fpath, 'significant_means.txt', sep = '/')
  result <- read.table(file.path, header = TRUE, sep = "\t",check.names = FALSE)
  result <- result[, c(2, 15:dim(result)[2])]
  
  result <- result %>% pivot_longer(cols = -interacting_pair, names_to = "sr", values_to = "LRscore")
  result <- dplyr::filter(result,  !is.na(LRscore))
  result <- tidyr::separate(data = result, col = sr, into = c("Sender", "Receiver"), sep = "\\|")
  
  # handle the complex information
  if(T){
    cpdb.complex <-read.csv("./cpdb/db/v5.0.0/complex_input.csv")
    rec.complex <- cpdb.complex$complex_name[which(cpdb.complex$receptor == TRUE)]
    lig.complex <- cpdb.complex$complex_name[which(cpdb.complex$receptor != TRUE)]
    cpdb.complex <- cpdb.complex[,1:5]
    cpdb.gene <- read.csv("./cpdb/db/v5.0.0/gene_input.csv")
    cpdb.gene <- dplyr::distinct(cpdb.gene, gene_name, uniprot, hgnc_symbol, .keep_all = TRUE)
    cpdb.gene <- cpdb.gene[,1:2]
    
    tmp.cpdb.complex <- merge(cpdb.gene, cpdb.complex, by.y = "uniprot_1", by.x = "uniprot")
    tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
    colnames(tmp.cpdb.complex)[1] <- "gene_1"
    tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_2", by.x = "uniprot")
    tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
    colnames(tmp.cpdb.complex)[1] <- "gene_2"
    tmp.cpdb.complex <- tmp.cpdb.complex[,-5]
    tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_3", by.x = "uniprot", all.y = TRUE)
    tmp.cpdb.complex <- tmp.cpdb.complex[-117,]
    tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
    colnames(tmp.cpdb.complex)[1] <- "gene_3"
    cpdb.complex <- tmp.cpdb.complex[,4:1]
    rm(cpdb.gene, tmp.cpdb.complex)
    
    cpdb.complex$gene <- paste(cpdb.complex$gene_1, cpdb.complex$gene_2, cpdb.complex$gene_3, sep = "&")
    cpdb.complex$gene <- gsub("&NA", "", cpdb.complex$gene)
    cpdb.complex <- cpdb.complex[,c("complex_name", "gene")]
  }
  
  # combine the complex information and result of cellphonedb
  if(T){
    complexes <- cpdb.complex$complex_name[which(stringr::str_count(cpdb.complex$complex_name, pattern = '_')>=1)]
    for (complex in complexes) {
      change.pair <- which(grepl(complex, result$interacting_pair))
      if(length(change.pair)>0){
        change.complex <- gsub('_', '*', complex)
        result$interacting_pair <- gsub(complex, change.complex, result$interacting_pair)
      }
    }
    
    result <- tidyr::separate(data = result, col = interacting_pair, into = c("Ligand", "Receptor"), sep = "_")
    result$Ligand <- gsub('\\*', '_', result$Ligand)
    result$Receptor <- gsub('\\*', '_', result$Receptor)
    result <- merge(result, cpdb.complex, by.x = "Ligand", by.y = "complex_name", all.x = TRUE)
    result$Ligand[!is.na(result$gene)] <- result$gene[!is.na(result$gene)]
    result <- result[,-6]
    result <- merge(result, cpdb.complex, by.x = "Receptor", by.y = "complex_name", all.x = TRUE)
    result$Receptor[!is.na(result$gene)] <- result$gene[!is.na(result$gene)]
    result <- result[,-6]
    
    result.rec <- which(result$Receptor %in% rec.complex)
    result.lig <- which(result$Ligand %in% lig.complex)
    
    if(length(result.lig)>0){
      result <- result[-result.lig,]
      print(paste0('Delate：', length(result.lig)))
    }else if(length(result.rec)>0){
      result <- result[-result.rec,]
      print(paste0('Delate：', length(result.rec)))
    }
    
    result <- result[which(result$Sender!=result$Receiver),]
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
    result <- dplyr::distinct(result, all, .keep_all = TRUE)
  }
  return(result)
}

MDIC3_function <- function(fpath){
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  suppressMessages(library(lobstr))
  set.seed(123)

  result <- read.csv(fpath, header = TRUE)
  result <- separate(result, 'L_R_pair', c('Ligand','Receptor'), sep = '\\ - ')
  result$Ligand <- gsub("\\(", "", result$Ligand)
  result$Ligand <- gsub("\\)", "", result$Ligand)
  result$Receptor <- gsub("\\(", "", result$Receptor)
  result$Receptor <- gsub("\\)", "", result$Receptor)
  
  result$Ligand <- gsub('\\+', '&', result$Ligand)
  result$Receptor <- gsub('\\+', '&', result$Receptor)
  result <- result[,1:5]
  colnames(result) <- c('Sender','Receiver','Ligand', 'Receptor', 'LRscore')
  result <- result[which(result$Sender != result$Receiver),]
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  return(result)
}

scTenifoldXct_function <- function(fpath){
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  suppressMessages(library(lobstr))
  set.seed(123)

  result <- read.csv(fpath, header = TRUE)
  result <- result[,c(1,2,4,7,8)]
  colnames(result) <- c('Ligand', 'Receptor', 'LRscore', 'Sender','Receiver')
  result <- result[which(result$Sender != result$Receiver),]
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  return(result)
}


