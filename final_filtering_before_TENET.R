################################################################################
#                                                                              #
#                            Filtering before TENET                            #
#                                                                              #
################################################################################

  ## Author: Ines Rivero-Garcia
  ## Date: 02/08/2022
  ## Last update: 15/08/2022
  ## This script filters the BL and CSD seurat object before running TENET

## 1. LOAD LIBRARIES ###########################################################
library(Seurat)
library(ggplot2)
library(matrixStats)
library(ggvenn)
library(UpSetR)
#library(pheatmap)
#library(RColorBrewer)

## 2. LOAD SEURAT OBJECTS AND PREPARE THEM FOR FILTERING #######################
setwd("Desktop/Axolotl_Limb_scGRN/Raw.data/Seurat_Objects/")
# Load seurat objects
BL <- readRDS("BL_SeuratObj.RDS")
CSD <- readRDS("CSD_SeuratObj.RDS")

# Focus on connective tissue
Idents(BL) <- "celltype"
Idents(CSD) <- "celltype"
BL <- subset(BL, idents = "connective tissue")
CSD <- subset(CSD, idents = "connective tissue")

# Remove non-injured cells from CSD
Idents(CSD) <- "orig.ident"
CSD <- subset(CSD, idents = "CSD_0dpa", invert = TRUE)

# Get matrix of raw and normalized counts
BLr <- GetAssayData(BL, assay = "RNA", slot = "counts")
CSDr <- GetAssayData(CSD, assay = "RNA", slot = "counts")
BLn <- GetAssayData(BL, assay = "RNA", slot = "data")
CSDn <- GetAssayData(CSD, assay = "RNA", slot = "data")


## 3. FILTER 1: KEEP GENES WITH > 1 LOG-NORM COUNTS IN > 5% CELLS ##############
# Which BL genes fulfill this condition?
countFilter <- 1
cellFilter <- round(ncol(BLn) * 0.05)
keepBL5 <- c()
for(i in 1:nrow(BLn)){
  k1 <- BLn[i,] > countFilter
  k2 <- sum(k1) > cellFilter
  keepBL5 <- c(keepBL5,k2)
}
table(keepBL5)

# Which CSD genes fulfill this condition?
cellFilter <- round(ncol(CSDn) * 0.05)
keepCSD5 <- c()
for(i in 1:nrow(CSDn)){
  k1 <- CSDn[i,] > countFilter
  k2 <- sum(k1) > cellFilter
  keepCSD5 <- c(keepCSD5,k2)
}
table(keepCSD5)

# Subset matrices 
BLr <- BLr[keepBL5, ]
BLn <- BLn[keepBL5, ]
CSDr <- CSDr[keepCSD5, ]
CSDn <- CSDn[keepCSD5, ]

# Create new seurat objects with the subseted matrices.
#BLmini <- CreateSeuratObject(BLr)
#BLmini@assays$RNA@data <- BLn
#CSDmini <- CreateSeuratObject(CSDr)
#CSDmini@assays$RNA@data <- CSDn
#BLmini@meta.data <- BL@meta.data
#CSDmini@meta.data <- CSD@meta.data
#limbmini <- merge(BLmini, CSDmini)
#limbmini$Injury <- ifelse(limbmini$orig.ident %in% c("BL_3_5dpa", "BL_8dpa", "BL_11dpa"),
#                          "BL", "CSD")

## 4. FILTER 2: REMOVE INVARIANT GENES #########################################
# Calculate z-scores.
    # NOTE:The scale function does column-wise scaling. 
    # We transpose the matrix to do row (gene)-wise scaling.
    # Then we transpose again to get the gene x cell matrices.
BLnz <- t(scale(t(as.matrix(BLn)), center = TRUE, scale = TRUE))
CSDnz <- t(scale(t(as.matrix(CSDn)), center = TRUE, scale = TRUE))

# Remove genes with |z-score| < 0.5 in all cells.
# Which BL genes fulfill this condition?
zFilter <- 0.5
cellFilter <- ncol(BLnz)
removeBL <- c()
for(i in 1:nrow(BLnz)){
  print(i)
  k1 <- abs(BLnz[i,]) < zFilter
  k2 <- sum(k1) == cellFilter
  removeBL <- c(removeBL,k2)
}
table(removeBL) # TRUE is genes to remove, FALSE is genes to keep.

# Which CSD genes fulfill this condition?
zFilter <- 0.5
cellFilter <- ncol(CSDnz)
removeCSD <- c()
for(i in 1:nrow(CSDnz)){
  k1 <- abs(CSDnz[i,]) < zFilter
  k2 <- sum(k1) == cellFilter
  removeCSD <- c(removeCSD,k2)
}
table(removeCSD) # TRUE is genes to remove, FALSE is genes to keep.

# Set difference between all genes and the ones which don't pass the z-score filter.
#keepBLz <- setdiff(rownames(BLnz, removeBL))
#keepCSDz <- setdiff(rownames(CSDnz, removeCSD))


# NO GENES TO BE REMOVED ACCORDING TO THIS CRITERIA!.
# Subset matrices 
#BLr <- BLr[keepBLz, ]
#BLnz <- BLn[keepBLz, ]
#CSDr <- CSDr[keepCSDz, ]
#CSDnz <- CSDn[keepCSDz, ]

# Save RDS with everything ?

## 5. SAVE DATA FOR TENET ######################################################
# Save raw count matrix as a csd with cells as rows and genes as columns.
BLnz <- t(BLnz)
CSDnz <- t(CSDnz)

write.csv(BLnz, "../TENET_input/2022-08-15_ConnectiveTissue_Injured_ExpFilter1logUMI-5perCell_zFilter0.5-allCell/2022-08-15_BL_TENET_ConnectiveTissueOnly_Injured_countMatrix.csv",
          row.names = TRUE, quote = FALSE)
write.csv(CSDnz, "../TENET_input/2022-08-15_ConnectiveTissue_Injured_ExpFilter1logUMI-5perCell_zFilter0.5-allCell/2022-08-15_CSD_TENET_ConnectiveTissueOnly_Injured_countMatrix.csv",
          row.names = TRUE, quote = FALSE)

# Save trajectory file: txt with time point of each cell.
timeBL <- BL@meta.data[rownames(BLnz), "orig.ident"]
timeCSD <- CSD@meta.data[rownames(CSDnz), "orig.ident"]

timeBL <- as.character(timeBL)
timeCSD <- as.character(timeCSD)

timeBL[timeBL == "BL_3_5dpa"] <- "1"
timeBL[timeBL == "BL_8dpa"] <- "2"
timeBL[timeBL == "BL_11dpa"] <- "3"
timeCSD[timeCSD == "CSD_3dpa"]<- "1"
timeCSD[timeCSD == "CSD_5dpa"] <- "2"
timeCSD[timeCSD == "CSD_8dpa"] <- "3"
timeCSD[timeCSD == "CSD_11dpa"] <- "4"

timeBL <- as.numeric(timeBL)
timeCSD <- as.numeric(timeCSD)

write.table(timeBL, "../TENET_input/2022-08-15_ConnectiveTissue_Injured_ExpFilter1logUMI-5perCell_zFilter0.5-allCell/2022-08-15_BL_TENET_ConnectiveTissueOnly_Injured_timeVector.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(timeCSD, "../TENET_input/2022-08-15_ConnectiveTissue_Injured_ExpFilter1logUMI-5perCell_zFilter0.5-allCell/2022-08-15_CSD_TENET_ConnectiveTissueOnly_Injured_timeVector.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save cell file: txt with a 1 for all cells (all cells should be considered).
cellBL <- rep(1, nrow(BLnz))
cellCSD <- rep(1, nrow(CSDnz))

write.table(cellBL, "../TENET_input/2022-08-15_ConnectiveTissue_Injured_ExpFilter1logUMI-5perCell_zFilter0.5-allCell/2022-08-15_BL_TENET_ConnectiveTissueOnly_Injured_cellVector.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cellCSD, "../TENET_input/2022-08-15_ConnectiveTissue_Injured_ExpFilter1logUMI-5perCell_zFilter0.5-allCell/2022-08-15_CSD_TENET_ConnectiveTissueOnly_Injured_cellVector.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

################################################################################
## 4. FILTER 2: REMOVE INVARIANT GENES #########################################
# Explanation: we discard those genes that are invariant over time in BL and CSD
# (separately), and that aren't differentially expressed in both injuries.
# Important gene sets
  # BLcandidates = genes that aren't significantly DE between all BL timepoints
  # CSDcandidates = genes that aren't significantly DE between all CSD timepoints
  # injury.candidates = genes that aren't significantly DE between CSD and BL
  # delete = intersect of BLcandidates, CSDcandidates and injury.candidates


# Find invariant genes over time in BL
#Idents(BLmini)<-"orig.ident"
#BL35_8.markers <- FindMarkers(BLmini, ident.1 = "BL_3_5dpa", ident.2 = "BL_3_5dpa",
#                              logfc.threshold = 0)
#BL8_11.markers<- FindMarkers(BLmini, ident.1 = "BL_8dpa", ident.2 = "BL_11dpa",
#                             logfc.threshold = 0)
#BL35_8.candidates <- rownames(BL35_8.markers[BL35_8.markers$p_val_adj > 0.05,])
#BL8_11.candidates <- rownames(BL8_11.markers[BL8_11.markers$p_val_adj > 0.05,])
#BLcandidates <- intersect(BL35_8.candidates, BL8_11.candidates)

# Find invariant genes over time in CSD
#Idents(CSDmini) <- "orig.ident"
#CSD3_5.markers <- FindMarkers(CSDmini, ident.1 = "CSD_3dpa", ident.2 = "CSD_5dpa",
#                              logfc.threshold = 0)
#CSD5_8.markers <- FindMarkers(CSDmini, ident.1 = "CSD_5dpa", ident.2 = "CSD_8dpa", 
#                              logfc.threshold = 0)
#CSD8_11.markers <- FindMarkers(CSDmini, ident.1 = "CSD_8dpa", ident.2 = "CSD_11dpa", 
#                               logfc.threshold = 0)
#CSD3_5.candidates <- rownames(CSD3_5.markers[CSD3_5.markers$p_val_adj > 0.05,])
#CSD5_8.candidates <- rownames(CSD5_8.markers[CSD5_8.markers$p_val_adj > 0.05,])
#CSD8_11.candidates <- rownames(CSD8_11.markers[CSD8_11.markers$p_val_adj > 0.05,])
#CSDcandidates <- intersect(CSD3_5.candidates, 
#                           intersect(CSD5_8.candidates, CSD8_11.candidates))

# Find not-DEGs in both injuries
#Idents(limbmini) <- "Injury"
#injury.markers <- FindMarkers(limbmini, ident.1 = "BL", ident.2 = "CSD", 
#                              logfc.threshold = 0)

# Find genes that (a) don't change over time in BL, (b) don't change over time
# in CSD and (c) aren't DE in BL vs CSD.
#delete<-intersect(BLcandidates, intersect(CSDcandidates, injury.candidates))

# Subset matrices


## 5. GENE EXPRESSION VISUALIZATION ############################################
#Limb <- readRDS("ConnectiveTissue_CSDandBL_Ines.RDS")
#Idents(Limb) <- "orig.ident"
#dictionary <- read.csv("../AnnotationFile_AmexG_v6_chr_unscaffolded_CherryGFP_v1.2.csv",
#                       header = FALSE, sep = ";")
#dictionary <- dictionary[,1:2]

# Visualize genes that pass the expression filter in BL (if a gene belongs to 
# delete genes, it's indicated in the plot title with a "(d")
#BLlist <- vector(mode = "list")

#for(i in 1:length(blgenes)){
#  print(i)
#  if(blgenes[i] %in% delete){
#    d <- " (d)"
#  }else{
#    d <- ""
#  }
#  gene <- paste0(dictionary[dictionary$V1 == blgenes[i], "V2"], d)
#  a <- FeaturePlot(Limb, features = blgenes[i], cols = c("grey", "red")) +
#    ggtitle(gene)
#  b<- VlnPlot(Limb, features = blgenes[i], split.by = "orig.ident") +
#    ggtitle(gene)
#  BLlist[[i]] <- a + b
#}

#BLlist2<-BLlist[2001:length(BLlist)]
#pdf("~/Desktop/BLcandidates.pdf", onefile = TRUE)
#BLlist2
#dev.off()

# Visualize genes that pass the expression filter in CSD (if a gene belongs to 
# delete genes, it's indicated in the plot title with a "(d")
#CSDlist <- vector(mode = "list")

#for(i in 1:length(csdgenes)){
#  print(i)
#  if(csdgenes[i] %in% delete){
#    d <- " (d)"
#  }else{
#    d <- ""
#  }
#  gene <- paste0(dictionary[dictionary$V1 == csdgenes[i], "V2"], d)
#  a <- FeaturePlot(Limb, features = csdgenes[i], cols = c("grey", "red")) +
#    ggtitle(gene)
#  b<- VlnPlot(Limb, features = csdgenes[i], split.by = "orig.ident") +
#    ggtitle(gene)
#  CSDlist[[i]] <- a + b
#}

#pdf("~/Desktop/CSDcandidates.pdf", onefile = TRUE)
#CSDlist
#dev.off()

# Visualize delete genes
#DELlist <- vector(mode = "list")

#for(i in 1:length(delete)){
#  print(i)
#  gene <- paste(dictionary[dictionary$V1 == delete[i], "V2"], "(d)")
#  a <- FeaturePlot(Limb, features = delete[i], cols = c("grey", "red")) +
#    ggtitle(gene)
#  b<- VlnPlot(Limb, features = delete[i], split.by = "orig.ident") +
#    ggtitle(gene)
#  DELlist[[i]] <- a + b
#}

#pdf("~/Desktop/delete.pdf", onefile = TRUE)
#DELlist
#dev.off()

# Intersect between conditions
#pdf("~/Desktop/Upset.pdf")
#upset(fromList(list(BL_3.5_8 = BL35_8.candidates, BL_8_11 = BL8_11.candidates,
#               CSD_3_5 = CSD3_5.candidates, CSD_5_8 = CSD5_8.candidates,
#               CSD_8_11 = CSD8_11.candidates, Injury = injury.candidates)),
#      nintersects = NA, order.by = "freq")
#dev.off()

# Timeline plot
#timeDF <- rbind(reshape2::melt(as.matrix(BLn[delete,])),
#                reshape2::melt(as.matrix(CSDn[delete,])))
#timeDF$Injury <- ifelse(grepl('^BL_', timeDF$Var2), "BL", "CSD")
#timeDF$time <- timeDF$Var2
#timeDF$time <- as.character(timeDF$time)
#timeDF[grepl('^BL_3', timeDF$time), "time"]<-4
#timeDF[grepl('BL_8', timeDF$time), "time"]<-8
#timeDF[grepl('^BL_11', timeDF$time), "time"]<- 11
#timeDF[grepl('^CSD_3', timeDF$time), "time"]<-3
#timeDF[grepl('^CSD_5', timeDF$time), "time"]<-5
#timeDF[grepl('^CSD_8', timeDF$time), "time"]<-8
#timeDF[grepl('^CSD_11', timeDF$time), "time"]<-11
#timeDF$time <- as.numeric(timeDF$time)
#colnames(timeDF) <- c("Gene", "Cell", "NormExpr", "Injury", "Time")

#ggplot(timeDF[timeDF$Gene %in% delete[61:70], ], aes(x = Time, y = NormExpr, color=Gene)) + 
#  stat_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=3)) +
#  theme_classic() + 
#  theme(legend.position = "none") + 
#  scale_color_manual(values = c("red", "yellow", "green", "blue", "orange", "cyan", "pink", "black", "purple", "brown")) + 
#  facet_grid(. ~ Injury) + ylim(c(0,7.5))+xlim(c(3,12))




################################################################################
#Mast with model ~ injury + time + injury:time

#Idents(biglimb)<- "orig.ident"


# DE
#Idents(limb)<-""
#BL35_8.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono")
