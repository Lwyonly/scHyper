CalSI_1 <- function(a, b) {
  intersection = length(intersect(unique(a), unique(b)))
  #union = length(a) + length(b) - intersection
  union <- min(length(unique(a)), length(unique(b)))
  return (intersection/union)
}

CalSI_2 <- function(result.path){
  methods <- c('scHyper', 'CellChat', 'CellPhoneDB', 'iTALK', 'SingleCellSignalR', 'MDIC3', 'scTenifoldXct')
  
  SI_res <- lapply(methods, function(method1){
    print(method1)
    result1 <- readRDS(paste0(result.path, '/', method1, '/result.rds'))
    #result1 <- IDTransform(data, method1, result1)
    temp.si <- lapply(methods, function(method2){
      result2 <- readRDS(paste0(result.path, '/', method2, '/result.rds'))
      #result2 <- IDTransform(data, method2, result2)
      if(dim(result1)[1]==0 | dim(result2)[1]==0){
        tmp.si <- NA
      }else{
        tmp.si <- CalSI_1(as.character(result1$all), as.character(result2$all))
      }
      
      tmp.si
    })
    temp.si <- data.frame(SI = unlist(temp.si), row.names = methods)
    temp.si
  })
  SI_res <- do.call(cbind, SI_res)
  colnames(SI_res) <- methods
  #SI_res <- SI_res[!apply(SI_res, 1, function(x){all(is.na(x))}), ]
  #SI_res <- SI_res[, !apply(SI_res, 2, function(x){all(is.na(x))})]
  return(SI_res)
}

CalRSI <- function(result.path){
  methods <- c('scHyper', 'CellChat', 'CellPhoneDB', 'iTALK', 'SingleCellSignalR', 'MDIC3', 'scTenifoldXct')
  
  RSI_res <- lapply(methods, function(method1){
    print(paste0('Method1: ', method1))
    result1 <- readRDS(paste0(result.path, '/', method1, '/result.rds')) # recode the LR score
    result1 <- as.data.frame(result1)
    result1$LRscore <- as.numeric(result1$LRscore)
    #result1 <- IDTransform(data, method1, result1)
    result1$rank <- rank(-result1$LRscore)
    
    temp.rsi <- lapply(methods, function(method2){
      print(paste0('Method2: ', method2))
      result2 <- readRDS(paste0(result.path, '/', method2, '/result.rds'))
      result2 <- as.data.frame(result2)
      result2$LRscore <- as.numeric(result2$LRscore)
      #result2 <- IDTransform(data, method2, result2)
      result2$rank <- rank(-result2$LRscore)
      
      if(dim(result1)[1]==0 | dim(result2)[1]==0){
        tmp.rsi <- NA
      }else{
        overlap.lr <- intersect(result1$all, result2$all)
        if(length(overlap.lr)!=0){
          rank1 <- result1[which(result1$all %in% overlap.lr),]$rank/dim(result1)[1]
          rank2 <- result2[which(result2$all %in% overlap.lr), ]$rank/dim(result2)[1]
          mean.rank <- mean(abs(rank1-rank2))
          tmp.rsi <- 1-mean.rank
        }else{tmp.rsi <- 0}
      }
      tmp.rsi
    })
    temp.rsi <- data.frame(RSI = unlist(temp.rsi), row.names = methods)
    temp.rsi
  })
  RSI_res <- do.call(cbind, RSI_res)
  colnames(RSI_res) <- methods
  #RSI_res <- RSI_res[!apply(RSI_res, 1, function(x){all(is.na(x))}), ]
  #RSI_res <- RSI_res[, !apply(RSI_res, 2, function(x){all(is.na(x))})]
  
  return(RSI_res)
}
