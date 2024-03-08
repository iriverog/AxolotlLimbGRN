################################################################################
#                                                                              #
#   PREPARING EXPRESSION DATA TO FIND LEF1 AND TCF2L2 CANDIDATE TARGET GENES   #
#                                                                              #
################################################################################

## Date: 2023/02/01
## Author: Ines Rivero Garcia
## This script extracts the expression of LEF1, TCF2L2 and other genes with > 1 
## log-norm counts in > 5% cells to predict interactions among them using CNNC.


## 1. LOAD DATA ################################################################
library(Seurat)
setwd("~/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/")

BL <- readRDS("Raw.data/BLct_CSDct_filter_Harmony_SeuratObj.RDS")
CSD <- readRDS("Raw.data/BLct_CSDct_filter_Harmony_SeuratObj.RDS")
Idents(BL)<-"orig.ident"
Idents(CSD) <- "orig.ident"
BL <- subset(BL, idents="BL_3_5dpa")
CSD <- subset(CSD, idents = "CSD_3dpa")
BLn <- as.matrix(GetAssayData(BL, assay = "RNA", slot = "data")) # For filtering
CSDn <- as.matrix(GetAssayData(CSD, assay = "RNA", slot = "data")) # For filtering


## 2. FIND GENES THAT PASS THE THE EXPRESSION FILTER? ##########################
countFilter <- 1
cellFilterBL <- round(ncol(BLn) * 0.05)
cellfilterCSD <- round(ncol(CSDn) * 0.05)

keepBL <- c()
keepCSD <- c()

for(i in 1:nrow(BLn)){
  print(i)
  k1 <- BLn[i,] > countFilter
  k2 <- sum(k1) > cellFilterBL
  keepBL <- c(keepBL,k2)
}
table(keepBL)

genesBL <- rownames(BLn)[keepBL]

for(i in 1:nrow(CSDn)){
  print(i)
  k1 <- CSDn[i,] > countFilter
  k2 <- sum(k1) > cellfilterCSD
  keepCSD <- c(keepCSD, k2)
}
table(keepCSD)

genesCSD <- rownames(CSDn)[keepCSD]


## 3. PREPARE TABLES OF GENE PAIRS FOR PREDICTION ##############################
# LEF1 only passes the expression filter in BL and TCF2L2 only passes the expression
# filter in CSD.
LEF1 <- "AMEX60DD044502"
TCF2L2 <- "AMEX60DD053152"

genePairsBL <- data.frame(Var1 = genesBL, Var2 = rep(LEF1, length(genesBL)))
genePairsCSD <- data.frame(Var1 = genesCSD, Var2 = rep(TCF2L2, length(genesCSD)))

# Remove the interaction of LEF1 with itself and TCF7L2 with itself.
keepBL <- 1:nrow(genePairsBL)
remove <- c()
for(i in 1:nrow(genePairsBL)){
  if(genePairsBL[i,1] == genePairsBL[i, 2]){
    remove <- c(remove, i)
  }
}
keepBL <- setdiff(keepBL, remove)
genePairsBL <- genePairsBL[keepBL,]

keepCSD <- 1:nrow(genePairsCSD)
remove <- c()
for(i in 1:nrow(genePairsCSD)){
  if(genePairsCSD[i,1] == genePairsCSD[i,2]){
    remove <- c(remove, i)
  }
}
keepCSD <- setdiff(keepCSD, remove)
genePairsCSD <- genePairsCSD[keepCSD,]


## 4. SAVE TABLES FOR PREDICTION ###############################################
write.table(genePairsBL, "Results/CNNC/LEF1_TCF7L2_Modeling/LEF1_BL3dpi_gene_pairs.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(genePairsCSD, "Results/CNNC/LEF1_TCF7L2_Modeling/TCF7L2_CSD3dpi_gene_pairs.txt",
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(0:nrow(genePairsBL), "Results/CNNC/LEF1_TCF7L2_Modeling/LEF1_BL3dpi_gene_pairs_num.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(0:nrow(genePairsCSD), "Results/CNNC/LEF1_TCF7L2_Modeling/TCF7L2_CSD3dpi_gene_pairs_num.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


## 5. SAVE EXPRESSION MATRIX FOR PREDICTION ####################################
BL <- NormalizeData(BL, normalization.method = "RC")
BLn <- as.matrix(GetAssayData(BL, assay = "RNA", slot = "data"))
BLn <- BLn[genesBL,]
BLn <- t(BLn)
write.table(BLn, "Results/CNNC/LEF1_TCF7L2_Modeling/LEF1_BL3dpi_NormalizedExpressionMatrix.tsv",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

CSD <- NormalizeData(CSD, normalization.method = "RC")
CSDn <- as.matrix(GetAssayData(CSD, assay = "RNA", slot = "data"))
CSDn <- CSDn[genesCSD,]
CSDn <- t(CSDn)
write.table(CSDn, "Results/CNNC/LEF1_TCF7L2_Modeling/TCF7L2_CSD3dpi_NormalizedExpressionMatrix.tsv",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
