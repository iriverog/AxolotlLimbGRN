################################################################################
#                                                                              #
#     PROCESS WNT PATHWAY PREDICTIONS BEFORE VISUALIZATION WITH CYTOSCAPE      #
#                                                                              #
################################################################################

setwd("~/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/CNNC/Wnt_Modeling/")

# Load Wnt network and predictions for each injury
Wnt = read.csv("FullWntNetwork_GeneSymbol_AmexID_SignOfInteractions.csv",
               header = TRUE, sep = ",")
WntBL3 = read.csv("BL3dpi_IntegratedWntNetwork_ExpressionFilter_PredictionDF.tsv", 
                  header = TRUE, sep = "\t")
WntCSD3 = read.csv("CSD3dpi_IntegratedWntNetwork_ExpressionFilter_PredictionDF.tsv", 
                   header = TRUE, sep = "\t")

# Add injury information
WntBL3$Input = "BL"
WntCSD3$Input = "CSD"

colnames(WntBL3) = c("SourceAMEX", "TargetAMEX", "BL3_Label", "BL3_Prob", "InputBL")
colnames(WntCSD3) = c("SourceAMEX", "TargetAMEX", "CSD3_Label", "CSD3_Prob", "InputCSD")

Wnt = merge(Wnt, WntBL3, all.x = TRUE, by = c("SourceAMEX", "TargetAMEX"))
Wnt = merge(Wnt, WntCSD3, all.x = TRUE, by = c("SourceAMEX", "TargetAMEX"))

# Tengo que cambiar los NAs por zeros si no no puedo comparar probabilidades.
Wnt[is.na(Wnt)] <- 0
Wnt$ProbBL_ProbCSD = Wnt$BL3_Prob - Wnt$CSD3_Prob

# Add cor in BL3 and CSD3
BL <- readRDS("../../../Raw.data/BLct_CSDct_filter_Harmony_SeuratObj.RDS")
CSD <- readRDS("../../../Raw.data/BLct_CSDct_filter_Harmony_SeuratObj.RDS")
Idents(BL) <- "orig.ident"
Idents(CSD) <- "orig.ident"
BL <- subset(BL, idents = "BL_3_5dpa")
CSD <- subset(CSD, idents = "CSD_3dpa")
BLn <- as.matrix(GetAssayData(BL, assay = "RNA", slot = "data"))
BLn <- BLn[unique(c(WntBL3$SourceAMEX, WntBL3$TargetAMEX)),]
CSDn <- as.matrix(GetAssayData(CSD, assay = "RNA", slot = "data"))
CSDn <- CSDn[unique(c(WntCSD3$SourceAMEX, WntCSD3$TargetAMEX)),]

BLn<-t(BLn)
CSDn <- t(CSDn)

corBL <- cor(BLn, method = "spearman")
corCSD <- cor(CSDn, method = "spearman")

Wnt$BL3_Cor <- NA
Wnt$CSD3_Cor <- NA

for(i in 1:nrow(Wnt)){
  gene1 = Wnt[i, "SourceAMEX"]
  gene2 = Wnt[i, "TargetAMEX"]
  
  if(gene1 %in% rownames(corBL)){
    if(gene2 %in% colnames(corBL)){
      Wnt[i, "BL3_Cor"] <- corBL[gene1, gene2]
    }
  }
  
  if(gene1 %in% rownames(corCSD)){
    if(gene2 %in% colnames(corCSD)){
      Wnt[i, "CSD3_Cor"] <- corCSD[gene1, gene2]
    }
  }
}


Wnt$Source_BL3_AvgCPM <- NA
Wnt$Target_BL3_AvgCPM  <- NA
Wnt$Source_CSD3_AvgCPM  <- NA
Wnt$Target_CSD3_AvgCPM  <- NA

BLn <- as.matrix(GetAssayData(BL, assay = "RNA", slot = "data"))
CSDn <- as.matrix(GetAssayData(CSD, assay = "RNA", slot = "data"))

for(i in 1:nrow(Wnt)){
  gene1 = Wnt[i,"SourceAMEX"]
  gene2 = Wnt[i,"TargetAMEX"]
  
  if(gene1 %in% rownames(BLn)){
    Wnt[i, "Source_BL3_AvgCPM"] <- mean(BLn[gene1,])
  }
  if(gene2 %in% rownames(BLn)){
    Wnt[i, "Target_BL3_AvgCPM"] <- mean(BLn[gene2,])
  }
  if(gene1 %in% rownames(CSDn)){
    Wnt[i, "Source_CSD3_AvgCPM"] <- mean(CSDn[gene1,])
  }
  if(gene2 %in% rownames(CSDn)){
    Wnt[i, "Target_CSD3_AvgCPM"] <- mean(CSDn[gene2,])
  }
}

write.table(Wnt[, c(1:7,9,10,12:14)], "IntegratedWntNetwork_BL3dpi_CSD3dpi_EdgeTable.tsv", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Get a mini DF with only those interactions that involve genes that pass the expression threshold in either BL or CSD.
a = rownames(Wnt[Wnt$InputBL == "BL",])
b = rownames(Wnt[Wnt$InputCSD == "CSD",])
keep=unique(c(a,b))
miniWnt = Wnt[keep,]
write.table(miniWnt[, c(1:7,9,10,12:14)], "IntegratedWntNetwork_BL3dpi_CSD3dpi_EdgeTable_Filtered.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


# Prepare node table
WntNode = data.frame(AMEX = unique(c(Wnt$SourceAMEX, Wnt$TargetAMEX)), Symbol = NA, BL3_AvgLogCPM = NA, CSD3_AvgLogCPM = NA)
for(i in 1:nrow(WntNode)){
  amex=WntNode[i,"AMEX"][1]
  if(amex %in% Wnt$SourceAMEX){
    WntNode[i, "Symbol"] = Wnt[Wnt$SourceAMEX == amex, "Source"][1]
    WntNode[i, "BL3_AvgLogCPM"] = Wnt[Wnt$SourceAMEX == amex, "Source_BL3_AvgCPM"][1]
    WntNode[i, "CSD3_AvgLogCPM"] = Wnt[Wnt$SourceAMEX == amex, "Source_CSD3_AvgCPM"][1]
  }else if(amex %in% Wnt$TargetAMEX){
    WntNode[i, "Target"] = Wnt[Wnt$TargetAMEX == amex, "Target"][1]
    WntNode[i, "BL3_AvgLogCPM"] = Wnt[Wnt$TargetAMEX == amex, "Target_BL3_AvgCPM"][1]
    WntNode[i, "CSD3_AvgLogCPM"] = Wnt[Wnt$TargetAMEX == amex, "Target_CSD3_AvgCPM"][1]
  }
}



write.table(WntNode, "IntegratedWntNetwork_BL3dpi_CSD3dpi_NodeTable.tsv", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
