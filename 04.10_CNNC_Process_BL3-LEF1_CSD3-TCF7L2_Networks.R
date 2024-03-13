################################################################################
#                                                                              #
#           PARSING BL3-5dpa LEF1 NETWORK AND CSD3dpa TCF7L2 NETWORK           #
#                                                                              #
################################################################################
library(ggplot2)
library(Seurat)

setwd("~/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/CNNC/LEF1_TCF7L2_Modeling/")
BL<-read.csv("BL3-LEF1_PredictionResults/BL3-LEF1_PredictionResults.tsv", 
             header = TRUE, sep = "\t")
CSD <- read.csv("CSD3-TCF7L2_PredictionResults/CSD3-TCF2L2_PredictionResult.tsv",
                header = TRUE, sep = "\t")
dictionary <- read.csv("../../../Raw.data/AnnotationFile_AmexG_v6_chr_unscaffolded_CherryGFP_Ines.tsv",
                       header = TRUE, sep = "\t")

#pdf("20230201_BL3-5-LEF1_and_CSD3-TCF7L2_histogramOfProbabilitiesZoomed.pdf")
ggplot() +
  geom_histogram(data = BL, aes(x = Prob_Label), fill = "forestgreen", alpha = 0.5, binwidth = 0.01) +
  geom_histogram(data = CSD, aes(x = Prob_Label), fill = "pink", alpha = 0.5, binwidth = 0.01) +
  theme_classic() +
  ggtitle("Distribution of predicted probabilities of interaction")
#dev.off()

BLconfident <- BL[BL$Pred_Label == 1,]
CSDconfident <-CSD[CSD$Pred_Label == 1,]

# Add correlations
BLgenes <- unique(c(BLconfident$To, BLconfident$From))
CSDgenes <-unique(c(CSDconfident$To, CSDconfident$From))

CT <- readRDS("../../../Raw.data/BLct_CSDct_filter_Harmony_SeuratObj.RDS")
Idents(CT) = "orig.ident"
BL <- subset(CT, idents = "BL_3_5dpa")
CSD <- subset(CT, idents = "CSD_3dpa")

CalculateTScore <- function(data){
  dataMEAN <- mean(data)
  dataSD <- sd(data)
  n <- length(data)
  return(sapply(data, function(x) (x-dataMEAN)/(dataSD/sqrt(n))))
}

BLn <- as.matrix(GetAssayData(BL, assay = "RNA", slot = "data"))
BLn <- BLn[BLgenes,]
BLt <- t(apply(BLn, 1, CalculateTScore))
CSDn <- as.matrix(GetAssayData(CSD, assay = "RNA", slot = "data"))
CSDn <- CSDn[CSDgenes,]
CSDt <- t(apply(CSDn, 1, CalculateTScore))

BLt_complete <- BLt[complete.cases(BLt),]
CSDt_complete <- CSDt[complete.cases(CSDt),]

BLconfident$BL3_Cor <- NA
CSDconfident$CSD3_Cor <- NA

for(i in 1:nrow(BLconfident)){
  g1 <- BLconfident[i, "From"]
  g2 <- BLconfident[i, "To"]
  BLconfident[i, "BL3_Cor"] <- cor(BLt_complete[g1,], BLt_complete[g2,], method = "spearman")
}

for(i in 1:nrow(CSDconfident)){
  g1 <- CSDconfident[i, "From"]
  g2 <- CSDconfident[i, "To"]
  CSDconfident[i, "CSD3_Cor"] <- cor(CSDt_complete[g1,], CSDt_complete[g2,], method = "spearman")
}

# Add mean CPM at CSD0, CSD3, CSD5, CSD8, CSD11
nodeDF <- data.frame(AMEXID = unique(c(BLconfident$From, BLconfident$To, CSDconfident$From, CSDconfident$To)),
                     Symbol = NA,
                     meanCPM_CSD0 = NA,
                     meanCPM_CSD3 = NA,
                     meanCPM_CSD5 = NA,
                     meanCPM_CSD8 = NA,
                     meanCPM_CSD11 = NA,
                     meanCPM_BL3 = NA,
                     meanCPM_BL8 = NA,
                     meanCPM_BL11 = NA)
# Add mean CPM at CSD0, BL3, BL8, BL11

# get genesets: TCFonly, LEFonly, common
LEFgenes <- BLconfident$From
TCF7L2genes <- CSDconfident$From
LEFonly <- setdiff(LEFgenes, TCF7L2genes)
LEFonlySymbol <- LEFonly
TCF7L2only <- setdiff(TCF7L2genes, LEFgenes)
TCF7L2onlySymbol <- TCF7L2only
commongenes <- intersect(TCF7L2genes, LEFgenes)
commongenesSymbol <- commongenes

for(i in 1:length(LEFonlySymbol)){
  a = LEFonlySymbol[i]
  g = dictionary[dictionary$AMEXID == a, "Symbol"]
  LEFonlySymbol[i] = g
}

for(i in 1:length(TCF7L2onlySymbol)){
  a = TCF7L2onlySymbol[i]
  g = dictionary[dictionary$AMEXID == a, "Symbol"]
  TCF7L2onlySymbol[i] = g
}

# GSEA
# DE in BL3 vs CSD3
keepBL <- read.csv("LEF1_BL3dpi_gene_pairs.txt", header = FALSE, sep = "\t")
keepCSD <- read.csv("TCF7L2_CSD3dpi_gene_pairs.txt", header = FALSE, sep = "\t")
keepGenes <- unique(c(keepBL$V1, keepBL$V2, keepCSD$V1, keepCSD$V2))

BL3_vs_CSD3_DEGs = FindMarkers(CT, ident.1 = "BL_3_5dpa", ident.2 = "CSD_3dpa", 
                        logfc.threshold = 0, test.use = "wilcox", only.pos = FALSE)
BL3_vs_CSD3_DEGs$LEF1only <- sapply(rownames(BL3_vs_CSD3_DEGs), function(x) is.element(x, LEFonly))
BL3_vs_CSD3_DEGs$TCF7L2only <- sapply(rownames(BL3_vs_CSD3_DEGs), function(x) is.element(x, TCF7L2only))
BL3_vs_CSD3_DEGs$LEF1_TCF7L2_commonGenes <- sapply(rownames(BL3_vs_CSD3_DEGs), function(x) is.element(x, commongenes))

# Add symbols!
BL3_vs_CSD3_DEGs$Symbol <- NA
for(i in 1:nrow(BL3_vs_CSD3_DEGs)){
  a <- rownames(BL3_vs_CSD3_DEGs)[i]
  g <- dictionary[dictionary$AMEXID == a, "Symbol"]
  BL3_vs_CSD3_DEGs[i, "Symbol"] <- g
  a <- NA
  g <- NA
}

# Add other information from biomaRt (human gene symbols)
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_info <- getBM(attributes=c("hgnc_symbol", "description", "gene_biotype"),
                      filters = 'hgnc_symbol', 
                      values = BL3_vs_CSD3_DEGs$Symbol, 
                      mart = ensembl)
BL3_vs_CSD3_DEGs$AMEXID <- rownames(BL3_vs_CSD3_DEGs)
BL3_vs_CSD3_DEGs <- merge(BL3_vs_CSD3_DEGs, gene_info, by.x = "Symbol", 
                          by.y = "hgnc_symbol", all.x = TRUE)

# Save results
write.table(BL3_vs_CSD3_DEGs, "DEGs_BL3_CSD3_LEF1-TCF7L2_annotation.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


# Venn diagram of target genes
library(ggvenn)
svg("../../../Figures/VennDiagram_LEF1_TCF7L2_Targets.svg", height = 8, width = 8)
ggvenn(data = list(LEF1 = LEFgenes, TCF7L2 = TCF7L2genes), show_percentage = FALSE, 
       fill_color = c("forestgreen", "hotpink3"), text_size = 8, set_name_size = 8)
dev.off()
