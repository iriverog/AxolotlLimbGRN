################################################################################
#                                                                              #
#          PREPARE EXPRESSION MATRICES TO MODEL WNT NETWORK WITH CNNC          #
#                                                                              #
################################################################################

setwd("~/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/CNNC/Wnt_Modeling/")

# Load Wnt network
Wnt <- read.csv("FullWntNetwork_GeneSymbol_AmexID_SignOfInteractions.csv",
                header = TRUE, sep = ",")

# Load genes that pass expression threshold in BL and CSD.
BLgenes3 <- readRDS("../../../Raw.data/BL_3-5dpa_Genes_More_1logCount_5PercentCells.rds")
CSDgenes3 <- readRDS("../../../Raw.data/CSD_3dpa_Genes_More_1logCount_5PercentCells.rds")

# Filter for BL network
keepBL <- c()
for(i in 1:nrow(Wnt)){
  if(Wnt[i, "SourceAMEX"] %in% BLgenes3){
    if(Wnt[i, "TargetAMEX"] %in% BLgenes3){
      keepBL <- c(keepBL, i)
    }
  }
}

keepCSD <- c()
for(i in 1:nrow(Wnt)){
  if(Wnt[i, "SourceAMEX"] %in% CSDgenes3){
    if(Wnt[i, "TargetAMEX"] %in% CSDgenes3){
      keepCSD <- c(keepCSD, i)
    }
  }
}

WntBL3 = Wnt[keepBL,]
WntCSD3 = Wnt[keepCSD,]

# Prepare data in the CNNC input format:
# We need 3 things: pairs of genes file, matrix with expression, numpairs
genePairsBL3 <- WntBL3[,c("SourceAMEX", "TargetAMEX")]
genePairsCSD3 <- WntCSD3[,c("SourceAMEX", "TargetAMEX")]

geneNumBL3 <- 0:nrow(genePairsBL3)+1
geneNumCSD3 <- 0:nrow(genePairsCSD3)+1

BL <- NormalizeData(BL, normalization.method = "RC")
BLn <- as.matrix(GetAssayData(BL, assay = "RNA", slot = "data"))
BLn <- BLn[unique(c(WntBL3$SourceAMEX, WntBL3$TargetAMEX)),]
BLn <- t(BLn)

CSD <- NormalizeData(CSD, normalization.method = "RC")
CSDn <- as.matrix(GetAssayData(CSD, assay = "RNA", slot = "data"))
CSDn <- CSDn[unique(c(WntCSD3$SourceAMEX, WntCSD3$TargetAMEX)),]
CSDn <- t(CSDn)

write.table(genePairsBL3, "BL3dpo_IntegratedWntNetwork_ExpressionFilter_GenePairs.tsv", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(genePairsCSD3, "CSD3dpi_IntegratedWntNetwork_ExpressionFilter_GenePairs.tsv",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(geneNumBL3, "BL3dpi_IntegratedWntNetwork_ExpressionFilter_GenePairsNum.txt", 
            col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(geneNumCSD3, "CSD3dpi_IntegratedWntNetwork_ExpressionFilter_GenePairsNum.txt", 
            col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

write.table(BLn, "BL3dpi_IntegratedWntNetwork_ExpressionFilter_ExpressionMatrix.tsv", 
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
write.table(CSDn, "CSD3dpi_IntegratedWntNetwork_ExpressionFilter_ExpressionMatrix.tsv", 
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")