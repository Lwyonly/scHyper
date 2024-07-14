rm(list = ls());gc()
source('./benchmark/methods.R')
source('./benchmark/RSIfunction.R')
set.seed(123)


data <- read.csv('./methods_result/pbmc3k/dataset/data.csv', header=TRUE, row.names = 1)
meta <- read.csv('./methods_result/pbmc3k/dataset/meta.csv', header=TRUE, row.names = 1)
result_CellChat <- CellChat_function(data, meta, species = 'human')
saveRDS(result_CellChat, file = paste0('./methods_result/pbmc3k/CellChat', '/result.rds'))
result_iTALK <- iTALK_function(data, meta, top_genes=50)
saveRDS(result_iTALK, file = paste0('./methods_result/pbmc3k/iTALK', '/result.rds'))
result_SCSR <- SCSR_function(data, meta, s.score = 0.5, species = 'human')
saveRDS(result_SCSR, file = paste0('./methods_result/pbmc3k/SingleCellSignalR', '/result.rds'))
result_CellPhoneDB <- CellPhoneDB2_function(fpath = './methods_result/pbmc3k/dataset')
saveRDS(result_CellPhoneDB, file = paste0('./methods_result/pbmc3k/CellPhoneDB', '/result.rds'))
result_scHyper <- scHyper_function(fpath = './methods_result/pbmc3k/dataset/df_enriched_copy.csv')
saveRDS(result_scHyper, file = paste0('./methods_result/pbmc3k/scHyper', '/result.rds'))
result_MDIC3 <- MDIC3_function(fpath="./methods_result/pbmc3k/MDIC3.csv")
saveRDS(result_MDIC3, file = paste0('./methods_result/pbmc3k/MDIC3', '/result.rds'))
result_scTenifoldXct <-scTenifoldXct_function(fpath="./methods_result/pbmc3k/dataset/scTenifoldXct.csv")
saveRDS(result_scTenifoldXct, file = paste0('./methods_result/pbmc3k/scTenifoldXct', '/result.rds'))


#Calculate RSI
RSIResult <- CalRSI(result.path = './methods_result/pbmc3k')
saveRDS(RSIResult, file = paste0('./methods_result/pbmc3k', 'RSIResult.rds'))

