subsetCommunication_internal <- function(net, LR, cells.level, slot.name = "net",
                                         sources.use = NULL, targets.use = NULL,
                                         signaling = NULL,
                                         pairLR.use = NULL,
                                         thresh = 0.05,
                                         datasets = NULL, ligand.pvalues = NULL, ligand.logFC = NULL, ligand.pct.1 = NULL, ligand.pct.2 = NULL,
                                         receptor.pvalues = NULL, receptor.logFC = NULL, receptor.pct.1 = NULL, receptor.pct.2 = NULL) {
  if (!is.data.frame(net)) {
    prob <- net$prob
    pval <- net$pval
    prob[pval >= thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source","target","interaction_name")
    net.pval <- reshape2::melt(pval, value.name = "pval")
    net$pval <- net.pval$pval
    # remove the interactions with zero values
    net <- subset(net, prob > 0)
  }
  if (!("ligand" %in% colnames(net))) {
    pairLR <- dplyr::select(LR, c("interaction_name_2", "pathway_name", "ligand",  "receptor"))
    idx <- match(net$interaction_name, rownames(pairLR))
    net <- cbind(net, pairLR[idx,])
  }
  
  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], pairLR.use = LR, key = "pathway_name", matching.exact = T, pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }
  
  if (!is.null(pairLR.use)){
    net <- tryCatch({
      subset(net,interaction_name %in% pairLR.use$interaction_name)
    }, error = function(e) {
      subset(net, pathway_name %in% pairLR.use$pathway_name)
    })
  }
  
  if (!is.null(datasets)) {
    if (!("datasets" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before selecting 'datasets'")
    }
    net <- net[net$datasets %in% datasets, , drop = FALSE]
  }
  if (!is.null(ligand.pvalues)){
    if (!("ligand.pvalues" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pvalues'")
    }
    net <- net[net$ligand.pvalues <= ligand.pvalues, , drop = FALSE]
  }
  if (!is.null(ligand.logFC)){
    if (!("ligand.logFC" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.logFC'")
    }
    if (ligand.logFC >= 0) {
      net <- net[net$ligand.logFC >= ligand.logFC, , drop = FALSE]
    } else {
      net <- net[net$ligand.logFC <= ligand.logFC, , drop = FALSE]
    }
  }
  if (!is.null(ligand.pct.1)){
    if (!("ligand.pct.1" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pct.1'")
    }
    net <- net[net$ligand.pct.1 >= ligand.pct.1, , drop = FALSE]
  }
  if (!is.null(ligand.pct.2)){
    if (!("ligand.pct.2" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pct.2'")
    }
    net <- net[net$ligand.pct.2 >= ligand.pct.2, , drop = FALSE]
  }
  
  if (!is.null(receptor.pvalues)){
    if (!("receptor.pvalues" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pvalues'")
    }
    net <- net[net$receptor.pvalues <= receptor.pvalues, , drop = FALSE]
  }
  if (!is.null(receptor.logFC)){
    if (!("receptor.logFC" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.logFC'")
    }
    if (receptor.logFC >= 0) {
      net <- net[net$receptor.logFC >= receptor.logFC, , drop = FALSE]
    } else {
      net <- net[net$receptor.logFC <= receptor.logFC, , drop = FALSE]
    }
  }
  if (!is.null(receptor.pct.1)){
    if (!("receptor.pct.1" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pct.1'")
    }
    net <- net[net$receptor.pct.1 >= receptor.pct.1, , drop = FALSE]
  }
  if (!is.null(receptor.pct.2)){
    if (!("receptor.pct.2" %in% colnames(net))) {
      stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'receptor.pct.2'")
    }
    net <- net[net$receptor.pct.2 >= receptor.pct.2, , drop = FALSE]
  }
  
  net <- net[rowSums(is.na(net)) != ncol(net), , drop = FALSE]
  
  if (nrow(net) == 0) {
    stop("No significant signaling interactions are inferred based on the input!")
  }
  
  
  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source","target","pathway_name","prob", "pval"))
    net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
    # net$source_target_pathway <- paste(paste(net$source, net$target, sep = "_"), net$pathway_name, sep = "_")
    net.pval <- net %>% group_by(source_target, pathway_name) %>% summarize(pval = mean(pval), .groups = 'drop')
    net <- net %>% group_by(source_target, pathway_name) %>% summarize(prob = sum(prob), .groups = 'drop')
    a <- stringr::str_split(net$source_target, "sourceTotarget", simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net <- dplyr::select(net, -source_target)
    net$pval <- net.pval$pval
  }
  
  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  
  net <- BiocGenerics::as.data.frame(net, stringsAsFactors=FALSE)
  
  if (nrow(net) == 0) {
    warning("No significant signaling interactions are inferred!")
  } else {
    rownames(net) <- 1:nrow(net)
  }
  
  if (slot.name == "net") {
    if (("ligand.logFC" %in% colnames(net)) & ("datasets" %in% colnames(net))) {
      net <- net[,c("source", "target", "ligand", "receptor",  "prob", "pval", "interaction_name", "interaction_name_2", "pathway_name",
                    "datasets","ligand.logFC", "ligand.pct.1", "ligand.pct.2", "ligand.pvalues",  "receptor.logFC", "receptor.pct.1", "receptor.pct.2", "receptor.pvalues")]
    } else if ("ligand.logFC" %in% colnames(net)) {
      net <- net[,c("source", "target", "ligand", "receptor",  "prob", "pval", "interaction_name", "interaction_name_2", "pathway_name",
                    "ligand.logFC", "ligand.pct.1", "ligand.pct.2", "ligand.pvalues",  "receptor.logFC", "receptor.pct.1", "receptor.pct.2", "receptor.pvalues")]
    } else {
      net <- net[,c("source", "target", "ligand", "receptor",  "prob", "pval", "interaction_name", "interaction_name_2", "pathway_name")]
    }
  } else if (slot.name == "netP") {
    net <- net[,c("source", "target", "pathway_name", "prob", "pval")]
  }
  
  return(net)
  
}

subsetCommunication <- function(validlrs, celltypes, net = net, slot.name = "net",
                                sources.use = NULL, targets.use = NULL,
                                signaling = NULL,
                                pairLR.use = NULL,
                                thresh = 0.05,
                                datasets = NULL, ligand.pvalues = NULL, ligand.logFC = NULL, ligand.pct.1 = NULL, ligand.pct.2 = NULL,
                                receptor.pvalues = NULL, receptor.logFC = NULL, receptor.pct.1 = NULL, receptor.pct.2 = NULL) {
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name = "netP"
    }
  }
  
  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }
  
  
  if (is.null(net)) {
    net <- net
  }
  LR <- validlrs
  cells.level <- celltypes
  df.net <- subsetCommunication_internal(net, LR, cells.level, slot.name = slot.name,
                                         sources.use = sources.use, targets.use = targets.use,
                                         signaling = signaling,
                                         pairLR.use = pairLR.use,
                                         thresh = thresh,
                                         datasets = datasets, ligand.pvalues = ligand.pvalues, ligand.logFC = ligand.logFC, ligand.pct.1 = ligand.pct.1, ligand.pct.2 = ligand.pct.2,
                                         receptor.pvalues = receptor.pvalues, receptor.logFC = receptor.logFC, receptor.pct.1 = receptor.pct.1, receptor.pct.2 =receptor.pct.2)
  
  return(df.net)
  
}
aggregateNet <- function(net, validlrs, celltypes, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05) {
  net <- net
  if (is.null(sources.use) & is.null(targets.use) & is.null(signaling) & is.null(pairLR.use)) {
    prob <- net$prob
    pval <- net$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net$count <- apply(prob > 0, c(1,2), sum)
    net$weight <- apply(prob, c(1,2), sum)
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  } else {
    
    df.net <- subsetCommunication(validlrs, celltypes, net = net, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source_target <- paste(df.net$source, df.net$target, sep = "_")
    df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = 'drop')
    df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = 'drop')
    df.net2$prob <- df.net3$prob
    a <- stringr::str_split(df.net2$source_target, "_", simplify = T)
    df.net2$source <- as.character(a[, 1])
    df.net2$target <- as.character(a[, 2])
    cells.level <- celltypes
    if (remove.isolate) {
      message("Isolate cell groups without any interactions are removed. To block it, set `remove.isolate = FALSE`")
      df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% unique(df.net2$source)])
      df.net2$target <- factor(df.net2$target, levels = cells.level[cells.level %in% unique(df.net2$target)])
    } else {
      df.net2$source <- factor(df.net2$source, levels = cells.level)
      df.net2$target <- factor(df.net2$target, levels = cells.level)
    }
    
    count <- tapply(df.net2[["count"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    net$count <- count
    net$weight <- prob
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  }
  
  return(net)
  
  
}


scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}
netVisual_circle <-function(net, color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, top = 1,
                            weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= "black",
                            edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                            edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2, vertex.size = NULL,
                            arrow.width=1,arrow.size = 0.2){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0
  
  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use)) ) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  
  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


searchPair <- function (signaling = c(), pairLR.use, key = c("pathway_name", "ligand"), matching.exact = FALSE, pair.only = TRUE){
  key <- match.arg(key)
  pairLR = future.apply::future_sapply(X = 1:length(signaling), 
                                       FUN = function(x) {
                                         if (!matching.exact) {
                                           index <- grep(signaling[x], pairLR.use[[key]])
                                         }
                                         else {
                                           index <- which(pairLR.use[[key]] %in% signaling[x])
                                         }
                                         if (length(index) > 0) {
                                           if (pair.only) {
                                             pairLR <- dplyr::select(pairLR.use[index, ], 
                                                                     interaction_name, pathway_name, ligand, receptor)
                                           }
                                           else {
                                             pairLR <- pairLR.use[index, ]
                                           }
                                           return(pairLR)
                                         }
                                         else {
                                           stop(cat(paste("Cannot find ", signaling[x], 
                                                          ".", "Please input a correct name!"), "\n"))
                                         }
                                       })
  if (pair.only) {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[c(4 * i - 3, 
                                             4 * i - 2, 4 * i - 1, 4 * i)]), ncol = 4, byrow = F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]][1:4]
    rownames(pairLR) <- pairLR[, 1]
  }
  else {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[(i * ncol(pairLR.use) - 
                                              (ncol(pairLR.use) - 1)):(i * ncol(pairLR.use))]), 
                             ncol = ncol(pairLR.use), byrow = F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]]
    rownames(pairLR) <- pairLR[, 1]
  }
  return(as.data.frame(pairLR, stringsAsFactors = FALSE))
}



netVisual_hierarchy1 <- function(net, vertex.receiver, color.use = NULL, title.name = NULL,  sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                 weight.scale = FALSE, vertex.weight=20, vertex.weight.max = NULL, vertex.size.max = NULL,
                                 edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6,
                                 label.dist = 2.8, space.v = 1.5, space.h = 1.6, shape= NULL, label.edge=FALSE,edge.curved=0, margin=0.2,
                                 vertex.label.cex=0.6,vertex.label.color= "black",arrow.width=1,arrow.size = 0.2,edge.label.color='black',edge.label.cex=0.5, vertex.size = NULL){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0
  cells.level <- rownames(net)
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  
  if (is.null(color.use)) {
    color.use <- scPalette(nrow(net))
  }
  
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+6
  
  m <- length(vertex.receiver)
  net2 <- net
  reorder.row <- c(vertex.receiver, setdiff(1:nrow(net),vertex.receiver))
  net2 <- net2[reorder.row,vertex.receiver]
  # Expand out to symmetric (M+N)x(M+N) matrix
  m1 <- nrow(net2); n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, m1+n1))
  
  row.names(net3) <- c(row.names(net)[vertex.receiver], row.names(net)[setdiff(1:m1,vertex.receiver)], rep("",m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,vertex.receiver)], rep("#FFFFFF",length(vertex.receiver)))
  color.use3.frame <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,vertex.receiver)], color.use[vertex.receiver])
  
  if (length(vertex.weight) != 1) {
    vertex.weight = c(vertex.weight[vertex.receiver], vertex.weight[setdiff(1:m1,vertex.receiver)],vertex.weight[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- c(rep("circle",m), rep("circle", m1-m), rep("circle",m))
  }
  
  g <- graph_from_adjacency_matrix(net3, mode = "directed", weighted = T)
  edge.start <- ends(g, es=E(g), names=FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m,1] <- 0; coords[(m+1):m1,1] <- space.h; coords[(m1+1):nrow(net3),1] <- space.h/2;
  coords[1:m,2] <- seq(space.v, 0, by = -space.v/(m-1)); coords[(m+1):m1,2] <- seq(space.v, 0, by = -space.v/(m1-m-1));coords[(m1+1):nrow(net3),2] <- seq(space.v, 0, by = -space.v/(n1-1));
  coords_scale<-coords
  
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    E(g)$label<-E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    # E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    E(g)$width<- 0.3+E(g)$weight/edge.weight.max*edge.width.max
  }else{
    E(g)$width<-0.3+edge.width.max*E(g)$weight
  }
  
  E(g)$arrow.width<-arrow.width
  E(g)$arrow.size<-arrow.size
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
  E(g)$color<-adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  
  label.dist <- c(rep(space.h*label.dist,m), rep(space.h*label.dist, m1-m),rep(0, nrow(net3)-m1))
  label.locs <- c(rep(-pi, m), rep(0, m1-m),rep(-pi, nrow(net3)-m1))
  # text.pos <- cbind(c(-space.h/1.5, space.h/10, space.h/1.2), space.v-space.v/10)
  text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), space.v-space.v/7)
  igraph::add.vertex.shape("fcircle", clip=igraph::igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin,rescale=T,vertex.shape="fcircle", vertex.frame.width = c(rep(1,m1), rep(2,nrow(net3)-m1)),
       vertex.label.degree=label.locs, vertex.label.dist=label.dist, vertex.label.family="Helvetica")
  text(text.pos, c("Source","Target","Source"), cex = 0.8, col = c("#c51b7d","#c51b7d","#2f6661"))
  arrow.pos1 <- c(-space.h/1.5, space.v-space.v/4, space.h/100000, space.v-space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v-space.v/4, space.h/20, space.v-space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], arrow.pos1[4], col = "#c51b7d",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], arrow.pos2[4], col = "#2f6661",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  if (!is.null(title.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1],title.pos[2],paste0(title.name, " signaling network"), cex = 1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}

scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}
netVisual_circle <-function(net, color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, top = 1,
                            weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= "black",
                            edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                            edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2, vertex.size = NULL,
                            arrow.width=1,arrow.size = 0.2){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0
  
  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use)) ) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  
  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}



searchPair <- function (signaling = c(), pairLR.use, key = c("pathway_name", "ligand"), matching.exact = FALSE, pair.only = TRUE){
  key <- match.arg(key)
  pairLR = future.apply::future_sapply(X = 1:length(signaling), 
                                       FUN = function(x) {
                                         if (!matching.exact) {
                                           index <- grep(signaling[x], pairLR.use[[key]])
                                         }
                                         else {
                                           index <- which(pairLR.use[[key]] %in% signaling[x])
                                         }
                                         if (length(index) > 0) {
                                           if (pair.only) {
                                             pairLR <- dplyr::select(pairLR.use[index, ], 
                                                                     interaction_name, pathway_name, ligand, receptor)
                                           }
                                           else {
                                             pairLR <- pairLR.use[index, ]
                                           }
                                           return(pairLR)
                                         }
                                         else {
                                           stop(cat(paste("Cannot find ", signaling[x], 
                                                          ".", "Please input a correct name!"), "\n"))
                                         }
                                       })
  if (pair.only) {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[c(4 * i - 3, 
                                             4 * i - 2, 4 * i - 1, 4 * i)]), ncol = 4, byrow = F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]][1:4]
    rownames(pairLR) <- pairLR[, 1]
  }
  else {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[(i * ncol(pairLR.use) - 
                                              (ncol(pairLR.use) - 1)):(i * ncol(pairLR.use))]), 
                             ncol = ncol(pairLR.use), byrow = F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]]
    rownames(pairLR) <- pairLR[, 1]
  }
  return(as.data.frame(pairLR, stringsAsFactors = FALSE))
}


mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}
netVisual_hierarchy2 <-function(net, vertex.receiver, color.use = NULL, title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                weight.scale = FALSE, vertex.weight=20, vertex.weight.max = NULL, vertex.size.max = NULL,
                                edge.weight.max = NULL, edge.width.max=8,alpha.edge = 0.6,
                                label.dist = 2.8, space.v = 1.5, space.h = 1.6, shape= NULL, label.edge=FALSE,edge.curved=0, margin=0.2,
                                vertex.label.cex=0.6,vertex.label.color= "black",arrow.width=1,arrow.size = 0.2,edge.label.color='black',edge.label.cex=0.5, vertex.size = NULL){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- levels(object@idents)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- levels(object@idents)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- levels(object@idents)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  
  
  if (is.null(color.use)) {
    color.use <- scPalette(nrow(net))
  }
  
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+6
  
  m <- length(vertex.receiver)
  m0 <- nrow(net)-length(vertex.receiver)
  net2 <- net
  reorder.row <- c(setdiff(1:nrow(net),vertex.receiver), vertex.receiver)
  net2 <- net2[reorder.row,vertex.receiver]
  # Expand out to symmetric (M+N)x(M+N) matrix
  m1 <- nrow(net2); n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, m1+n1))
  row.names(net3) <- c(row.names(net)[setdiff(1:m1,vertex.receiver)],row.names(net)[vertex.receiver],  rep("",m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[setdiff(1:m1,vertex.receiver)],color.use[vertex.receiver],  rep("#FFFFFF",length(vertex.receiver)))
  color.use3.frame <- c(color.use[setdiff(1:m1,vertex.receiver)], color.use[vertex.receiver], color.use[vertex.receiver])
  
  
  if (length(vertex.weight) != 1) {
    vertex.weight = c(vertex.weight[setdiff(1:m1,vertex.receiver)], vertex.weight[vertex.receiver], vertex.weight[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- rep("circle",nrow(net3))
  }
  
  g <- graph_from_adjacency_matrix(net3, mode = "directed", weighted = T)
  edge.start <- ends(g, es=igraph::E(g), names=FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m0,1] <- 0; coords[(m0+1):m1,1] <- space.h; coords[(m1+1):nrow(net3),1] <- space.h/2;
  coords[1:m0,2] <- seq(space.v, 0, by = -space.v/(m0-1)); coords[(m0+1):m1,2] <- seq(space.v, 0, by = -space.v/(m1-m0-1));coords[(m1+1):nrow(net3),2] <- seq(space.v, 0, by = -space.v/(n1-1));
  coords_scale<-coords
  
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    # E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<-adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  
  label.dist <- c(rep(space.h*label.dist,m), rep(space.h*label.dist, m1-m),rep(0, nrow(net3)-m1))
  label.locs <- c(rep(-pi, m0), rep(0, m1-m0),rep(-pi, nrow(net3)-m1))
  #text.pos <- cbind(c(-space.h/1.5, space.h/10, space.h/1.2), space.v-space.v/10)
  text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), space.v-space.v/7)
  igraph::add.vertex.shape("fcircle", clip=igraph::igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin,rescale=T,vertex.shape="fcircle", vertex.frame.width = c(rep(1,m1), rep(2,nrow(net3)-m1)),
       vertex.label.degree=label.locs, vertex.label.dist=label.dist, vertex.label.family="Helvetica")
  text(text.pos, c("Source","Target","Source"), cex = 0.8, col = c("#c51b7d","#2f6661","#2f6661"))
  
  arrow.pos1 <- c(-space.h/1.5, space.v-space.v/4, space.h/100000, space.v-space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v-space.v/4, space.h/20, space.v-space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], arrow.pos1[4], col = "#c51b7d",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], arrow.pos2[4], col = "#2f6661",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  
  if (!is.null(title.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1],title.pos[2],paste0(title.name, " signaling network"), cex = 1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}

netVisual_chord_cell_internal <- function(net, color.use = NULL, group = NULL, cell.order = NULL,
                                          sources.use = NULL, targets.use = NULL,
                                          lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                          remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                          transparency = 0.4, link.border = NA,
                                          title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,...){
  if (inherits(x = net, what = c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source","target")
  } else if (is.data.frame(net)) {
    if (all(c("source","target", "prob") %in% colnames(net)) == FALSE) {
      stop("The input data frame must contain three columns named as source, target, prob")
    }
    cell.levels <- as.character(union(net$source,net$target))
  }
  if (!is.null(cell.order)) {
    cell.levels <- cell.order
  }
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)
  
  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- cell.levels[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- cell.levels[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  # remove the interactions with zero values
  net <- subset(net, prob > 0)
  if(dim(net)[1]<=0){message("No interaction between those cells")}
  # create a fake data if keeping the cell types (i.e., sectors) without any interactions
  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels, as.character(union(net$source,net$target)))
    if (length(cells.removed) > 0) {
      net.fake <- data.frame(cells.removed, cells.removed, 1e-10*sample(length(cells.removed), length(cells.removed)))
      colnames(net.fake) <- colnames(net)
      net <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      if(nrow(net) > nrow(net.fake)){
        link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
      }
      # directional <- net[, 1:2]
      # directional$plot <- 0
      # directional$plot[1:(nrow(net) - nrow(net.fake))] <- 1
      # link.arr.type = "big.arrow"
      # message("Set scale = TRUE when remove.isolate = FALSE")
      scale = TRUE
    }
  }
  
  df <- net
  cells.use <- union(df$source,df$target)
  
  # define grid order
  order.sector <- cell.levels[cell.levels %in% cells.use]
  
  # define grid color
  if (is.null(color.use)){
    color.use = scPalette(length(cell.levels))
    names(color.use) <- cell.levels
  } else if (is.null(names(color.use))) {
    names(color.use) <- cell.levels
  }
  grid.col <- color.use[order.sector]
  names(grid.col) <- order.sector
  
  # set grouping information
  if (!is.null(group)) {
    group <- group[names(group) %in% order.sector]
  }
  
  # define edge color
  edge.color <- color.use[as.character(df$source)]
  
  if (directional == 0 | directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }
  
  circos.clear()
  chordDiagram(df,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
               link.arr.type = link.arr.type, # link.border = "white",
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = scale,
               group = group,
               link.target.prop = link.target.prop,
               reduce = reduce,
               ...)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)
  
  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(grid.col), type = "grid", legend_gp = grid::gpar(fill = grid.col), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }
  
  if(!is.null(title.name)){
    # title(title.name, cex = 1)
    text(-0, 1.02, title.name, cex=1)
  }
  circos.clear()
  gg <- recordPlot()
  return(gg)
}
netVisual_aggregate <- function(groupSize, pairLR.use, prob, pval, signaling, signaling.name = NULL, color.use = NULL, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE,
                                vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
                                weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8,
                                layout = c("circle","hierarchy","chord"), thresh = 0.05, from = NULL, to = NULL, bidirection = NULL, vertex.size = NULL,
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8,
                                group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20,
                                ...) {
  layout <- match.arg(layout)
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.weight)) {
    vertex.weight <- groupSize
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = pairLR.use, key = "pathway_name", matching.exact = T, pair.only = T)
  
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  
  
  pairLR.use.name <- dimnames(prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- prob
  pval <- pval
  
  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }
  
  
  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  nRow <- length(pairLR.name.use)
  
  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]
  
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
  # prob <-(prob-min(prob))/(max(prob)-min(prob))
  
  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.null(edge.weight.max)) {
      edge.weight.max = max(prob.sum)
    }
    par(mfrow=c(1,2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
    # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
    # grid.echo()
    # gg <-  grid.grab()
    gg <- recordPlot()
  } else if (layout == "circle") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    gg <- netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
  } else if (layout == "chord") {
    prob.sum <- apply(prob, c(1,2), sum)
    gg <- netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                        group = group, cell.order = cell.order,
                                        lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                        scale = scale, reduce = reduce,
                                        title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y= legend.pos.y)
  }
  
  return(gg)
  
}



netAnalysis_contribution <- function(pairLR.use, prob, pval, signaling, signaling.name = NULL, width = 0.1, vertex.receiver = NULL, thresh = 0.05, return.data = FALSE,
                                     x.rotation = 0, title = "Contribution of each L-R pair",
                                     font.size = 10, font.size.title = 10) {
  pairLR <- searchPair(signaling = signaling, pairLR.use = pairLR.use, key = "pathway_name", matching.exact = T, pair.only = T)
  pair.name.use = select(pairLR.use[rownames(pairLR),],"interaction_name_2")
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  
  pairLR.use.name <- dimnames(prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- prob
  pval <- pval
  
  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }
  
  
  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  
  prob <- prob[,,pairLR.name.use]
  
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    dimnames(prob)[3] <- pairLR.name.use
  }
  prob <-(prob-min(prob))/(max(prob)-min(prob))
  
  if (is.null(vertex.receiver)) {
    pSum <- apply(prob, 3, sum)
    pSum.max <- sum(prob)
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    y.lim <- max(pSum)
    
    pair.name <- unlist(dimnames(prob)[3])
    pair.name <- factor(pair.name, levels = unique(pair.name))
    if (!is.null(pairLR.name.use)) {
      pair.name <- pair.name.use[as.character(pair.name),1]
      pair.name <- factor(pair.name, levels = unique(pair.name))
    }
    mat <- pSum
    df1 <- data.frame(name = pair.name, contribution = mat)
    if(nrow(df1) < 10) {
      df2 <- data.frame(name = as.character(1:(10-nrow(df1))), contribution = rep(0, 10-nrow(df1)))
      df <- rbind(df1, df2)
    } else {
      df <- df1
    }
    df <- df[order(df$contribution, decreasing = TRUE), ]
    # df$name <- factor(df$name, levels = unique(df$name))
    df$name <- factor(df$name,levels=df$name[order(df$contribution, decreasing = TRUE)])
    df1$name <- factor(df1$name,levels=df1$name[order(df1$contribution, decreasing = TRUE)])
    gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity", width = 0.7) +
      theme_classic() + theme(axis.text.y = element_text(angle = x.rotation, hjust = 1,size=font.size, colour = 'black'), axis.text=element_text(size=font.size),
                              axis.title.y = element_text(size= font.size), axis.text.x = element_blank(), axis.ticks = element_blank()) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim) + coord_flip() + theme(legend.position="none") +
      scale_x_discrete(limits = rev(levels(df$name)), labels = c(rep("", max(0, 10-nlevels(df1$name))),rev(levels(df1$name))))
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5, size = font.size.title))
    }
    gg
    
  } else {
    pair.name <- factor(unlist(dimnames(prob)[3]), levels = unique(unlist(dimnames(prob)[3])))
    # show all the communications
    pSum <- apply(prob, 3, sum)
    pSum.max <- sum(prob)
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    y.lim <- max(pSum)
    
    df<- data.frame(name = pair.name, contribution = pSum)
    gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity",width = 0.2) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8),
                              axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("All")+ theme(plot.title = element_text(hjust = 0.5))#+
    
    # show the communications in Hierarchy1
    if (dim(prob)[3] > 1) {
      pSum <- apply(prob[,vertex.receiver,], 3, sum)
    } else {
      pSum <- sum(prob[,vertex.receiver,])
    }
    
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    
    df<- data.frame(name = pair.name, contribution = pSum)
    gg1 <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity",width = 0.2) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8), axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("Hierarchy1") + theme(plot.title = element_text(hjust = 0.5))#+
    #scale_x_discrete(limits = c(0,1))
    
    # show the communications in Hierarchy2
    
    if (dim(prob)[3] > 1) {
      pSum <- apply(prob[,setdiff(1:dim(prob)[1],vertex.receiver),], 3, sum)
    } else {
      pSum <- sum(prob[,setdiff(1:dim(prob)[1],vertex.receiver),])
    }
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    
    df<- data.frame(name = pair.name, contribution = pSum)
    gg2 <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity", width=0.9) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8), axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("Hierarchy2")+ theme(plot.title = element_text(hjust = 0.5))#+
    #scale_x_discrete(limits = c(0,1))
    title <- cowplot::ggdraw() + cowplot::draw_label(paste0("Contribution of each signaling in ", signaling.name, " pathway"), fontface='bold', size = 10)
    gg.combined <- cowplot::plot_grid(gg, gg1, gg2, nrow = 1)
    gg.combined <- cowplot::plot_grid(title, gg.combined, ncol = 1, rel_heights=c(0.1, 1))
    gg <- gg.combined
    gg
  }
  if (return.data) {
    df <- subset(df, contribution > 0)
    return(list(LR.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}


library(ComplexHeatmap)
netVisual_heatmap <- function(net, measure = c("count", "weight"), signaling = NULL, slot.name = c("net"), color.use = NULL, color.heatmap = c("#2166ac","#b2182b"),
                              title.name = NULL, width = NULL, height = NULL, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                              sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL){
  # obj1 <- object.list[[comparison[1]]]
  # obj2 <- object.list[[comparison[2]]]
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  # slot.name <- match.arg(slot.name)
  
  message("Do heatmap based on a single object \n")
  if (!is.null(signaling)) {
    net.diff <- net$prob[,,signaling]
    if (is.null(title.name)) {
      title.name = paste0(signaling, "network")
    }
    legend.name <- "Communication Prob."
  } else if (!is.null(measure)) {
    net.diff <- net[[measure]]
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Interaction strength"
      }
    }
    legend.name <- title.name
  }
  
  
  net <- net.diff
  
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show]
    color.use <- color.use[col.show]
  }
  
  
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), 0, round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
    # color.heatmap.use = colorRamp3(c(seq(min(mat), -(max(mat)-min(max(mat)))/9, length.out = 4), 0, seq((max(mat)-min(max(mat)))/9, max(mat), length.out = 4)), RColorBrewer::brewer.pal(n = 9, name = color.heatmap))
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
  }
  # col_fun(as.vector(mat))
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size+4),column_names_gp = gpar(fontsize = font.size+4),
                # width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title.name,column_title_gp = gpar(fontsize = font.size.title+4),column_names_rot = 90,
                row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title+4),row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, #at = colorbar.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}


netVisual_bubble <- function(net, celltypes, validlrs, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, color.heatmap = c("Spectral","viridis"), n.colors = 10, direction = -1, thresh = 0.05,
                             comparison = NULL, group = NULL, remove.isolate = FALSE, max.dataset = NULL, min.dataset = NULL,
                             min.quantile = 0, max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, color.text = NULL,
                             title.name = NULL, font.size = 10, font.size.title = 10, show.legend = TRUE,
                             grid.on = TRUE, color.grid = "grey90", angle.x = 90, vjust.x = NULL, hjust.x = NULL,
                             return.data = FALSE){
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(net[[1]])) {
    message("Comparing communications on a merged object \n")
  } else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 1, 1)
    vjust=c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
    })
  } else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  
  if (is.null(comparison)) {
    cells.level <- celltypes
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(validlrs, celltypes, net = net, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }
    
    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)
    
    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T)*1.1, max(df.net$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    # rownames(df.net) <- df.net$interaction_name_2
    
    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")
    
    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2),])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, levels = cells.order)
    df <- df.net
  }
  
  min.cutoff <- quantile(df$prob, min.quantile,na.rm= T)
  max.cutoff <- quantile(df$prob, max.quantile,na.rm= T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  
  
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    # line.on <- FALSE
    # df <- df[!is.na(df$prob),]
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, ,drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        #idx.na <- c(which(is.na(values)), which(!(dataset.name[comparison] %in% df.i.j$dataset)))
        dataset.na <- c(df.i.j$dataset[is.na(values)], setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          } else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            } else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
    #df <- df[!is.na(df$prob), ]
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target),unique(df$source.target)))
  
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
    geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")
  
  values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
  g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
  #g <- g + scale_radius(range = c(1,3), breaks = values,labels = names(values), name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  } else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }
  
  g <- g + theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  
  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept=seq(1.5, length(unique(df$source.target))-0.5, 1),lwd=0.1,colour=color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour=color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5+length(dataset.name[comparison]), length(group.names0)*length(dataset.name[comparison]), by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype="dashed", color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      } else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      #names(color) <- dataset.name[comparison]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order)-1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  } else {
    return(g)
  }
  
}






