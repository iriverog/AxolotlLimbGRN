################################################################################
#                                                                              #
#         TRANSFORMING JUNGKUI'S GRN MATRIX INTO CNNC GENE_PAIRS_TABLE         #
#                                                                              #
################################################################################

## Date: 2023/01/20
## Author: Ines Rivero Garcia
## This script transforms Jingkui's GRN PKN (mUA day5, day9, day13 proximal and
## day13 distal) into the files that CNNC requires for training

## 1. LOAD DATA ################################################################
setwd("/data_lab_MT/Ines/Tanaka/BoneRegeneration/")
PKN <- readRDS("Raw.data/TFdf_withAMEXIDs_someHoxMissing.rds")
BL <- readRDS("Raw.data/ConnectiveTissue_CSDandBL_Ines.RDS")
Idents(BL) <- "orig.ident"
BL <- subset(BL, idents = c("BL_3_5dpa", "BL_8dpa", "BL_11dpa"))
Dictionary <- read.csv("Raw.data/AnnotationFile_AmexG_v6_chr_unscaffolded_CherryGFP_Ines.tsv",
                       header = TRUE, sep = "\t")

# Check matrix dimensions: ~ 20,000 genes and ~ 200 TFs
dim(PKN)

# check annotation
head(colnames(PKN))
head(PKN)

## 2. MODIFYING THE TABLE FOR CNNC #############################################
## 2.1 GENERAL STRUCTURE AND CONCERNS.
# It should follow this structure: GeneA  GeneB p (\t separated).
# p is 0 if there's no relationship between the genes.
# p is 1 if GeneA regulates GeneB
# p is 2 if GeneB regulates GeneA.

# We want to make the training as balanced as possible. Therefore, for each TF
# we will write half of its targets as GeneA = Tf, GeneB = target, p = 1. The
# other half will be written as GeneA = target, GeneB = TF and p = 2. Random pairs
# of the TF with non-target genes will be generated in the same proportions.
TFs <- unique(PKN$Regulator)
nrTargets <- sapply(TFs, function(x) length(unique(PKN[PKN$Regulator == x, "Gene"])))
## 2.6 DOWNSAMPLING
# I have realized that I have almost 1M images for training and that's sooo much.
nrTargetsShort <- round(nrTargets/5)
names(nrTargets) <- TFs
names(nrTargetsShort) <- TFs

targets <- unique(rownames(BL))

## 2.2 TABLE CONSTRUCTION
CNNCtable <- as.data.frame(matrix(NA, nrow = 1, ncol = 3))
colnames(CNNCtable) <- c("GeneA", "GeneB", "p")

## 2.3 ADDING POSITIVE TRAINING EXAMPLES
for(tf in TFs){
  print(tf)
  # Extract PKN specific for that TF
  miniPKN <- PKN[PKN$Regulator == tf,]
  r <- sample(1:nrow(miniPKN), nrTargetsShort[tf], replace = FALSE)
  # Count nr of interactions and split in half for balance.
  rr <- sample(r, round(length(r)/2), replace = FALSE)
  rrr <- setdiff(r, rr)
  miniCNNC1 <- data.frame(GeneA = tf, GeneB = miniPKN[rr, 2], p = 1)
  miniCNNC2 <- data.frame(GeneA = miniPKN[rrr, 2], GeneB = tf, p = 2)
  
  CNNCtable <- rbind(CNNCtable, miniCNNC1)
  CNNCtable <- rbind(CNNCtable, miniCNNC2)
}

## 2.4 ADDING NEGATIVE EXAMPLES
# For each TF, randomly sample n genes. n is the number of + interactions for that TF
# in PKN.
# Check that none of those n genes is actually a + interaction. If so, remove it 
# and re-sample as needed.
# Split n in two. Write one half as TF target 0 and the other half as target TF 0

for(tf in TFs){
  print(tf)
  l <- 0
  n <- nrTargetsShort[tf]
  realTargets <- unique(PKN[PKN$Regulator == tf, "Gene"])
  candidates <- sample(targets, n, replace = FALSE)
  while (l != nrTargetsShort[tf]) {
    keep <- 1:length(candidates)
    for(g in candidates){
      if(g %in% realTargets){
        keep <- setdiff(keep, which(candidates == g))
      }else{
        next()
      }
    }
    candidates <- candidates[keep]
    l <- length(candidates)
    if(l < nrTargetsShort[tf]){
      candidates <- c(candidates, sample(targets, nrTargetsShort[tf]-l, 
                                         replace = FALSE))
    }
  }
  rr <- round(length(candidates)/2)
  rrr <- rr + 1
  
  miniCNNC1 <- data.frame(GeneA = tf, GeneB = candidates[1:rr], p = 0)
  miniCNNC2 <- data.frame(GeneA = candidates[rrr:length(candidates)], GeneB = tf, p = 0)
  
  CNNCtable <- rbind(CNNCtable, miniCNNC1)
  CNNCtable <- rbind(CNNCtable, miniCNNC2)
}


## 2.5 POLISHING
# Remove empty row used for initialization
CNNCtable <- CNNCtable[complete.cases(CNNCtable),]

# Order table by TF name 
TFrownames <- vector(mode = "list", length = length(TFs))
names(TFrownames) <- TFs
for(tf in TFs){
  print(tf)
  TFrownames[[tf]] <- c(rownames(CNNCtable[CNNCtable$GeneA == tf & CNNCtable$p == 1,]), 
                        rownames(CNNCtable[CNNCtable$GeneB == tf & CNNCtable$p == 2,]),
                        rownames(CNNCtable[CNNCtable$GeneA == tf & CNNCtable$p == 0,]), 
                        rownames(CNNCtable[CNNCtable$GeneB == tf & CNNCtable$p == 0,]))
  nrofrows <- length(TFrownames[[tf]])
  TFrownames[[tf]] <- sample(TFrownames[[tf]], nrofrows, replace = FALSE)
  print(length(TFrownames[[tf]]))
  print(length(unique(TFrownames[[tf]])))
}
CNNCtable <- CNNCtable[unlist(TFrownames),]

# Rename rows.
rownames(CNNCtable) <- 1:nrow(CNNCtable)

# Change TF symbols to AMEXIDs!!!!
columnas <- c(1,2)
for(c in columnas){
  for(r in 1:nrow(CNNCtable)){
    if(CNNCtable[r,c] %in% c("ZNF148", "NRF1", "RXRA", "ZNF454", "HOXD9", "VDR",
                             "ZNF774", "MEF2D", "HOXA9", "POU2F2")){
      next()
    }else{
      if(CNNCtable[r,c] %in% Dictionary$AMEXID){
        next()
      }else{
        id <- CNNCtable[r, c]
        amex <- Dictionary[Dictionary$Symbol == id, "AMEXID"][1]
        CNNCtable[r,c] <- amex
      }
    }
  }
}

# Remove rows that lack AMEXID
keep <- 1:nrow(CNNCtable)
remove <- c()
for(i in 1:nrow(CNNCtable)){
  if(CNNCtable[i, "GeneA"] %in% Dictionary$AMEXID){
    next()
  }else{
    remove <- c(remove, i)
  }
  if(CNNCtable[i, "GeneB"] %in% Dictionary$AMEXID){
    next()
  }else{
    remove <- c(remove, i)
  }
}

keep <- setdiff(keep, remove)
CNNCtable <- CNNCtable[keep,]
rownames(CNNCtable) <- 1:nrow(CNNCtable)

# Remove rows with genes that are not in the expression data 
goodgenes <- rownames(subset(BL, idents = c("BL_8dpa", "BL_11dpa")))
keep <- 1:nrow(CNNCtable)
remove <- c()
for(i in 1:nrow(CNNCtable)){
  if(CNNCtable[i, "GeneA"] %in% goodgenes){
    print("okA")
  }else{
    remove <- c(remove,i)
  }
  if(CNNCtable[i, "GeneB"] %in% goodgenes){
    print("okB")
  }else{
    remove <- c(remove,i)
  }
}

keep <- setdiff(keep, remove)
CNNCtable <- CNNCtable[keep,]
rownames(CNNCtable) <- 1:nrow(CNNCtable)

# Save as .txt tabseparated no headers.
write.table(CNNCtable, "CNNC_cluster2_network/BL8-11_gene_pairs.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# Save splits as txt for training (last row is 153043)
splits <- c(0, 605, 1389, 3063, 3531, 3960, 5325, 8210, 8726, 9273, 11320, 11932,
            12999, 16051, 17435, 20183, 20712, 21285, 21830, 23354, 24533, 24946,
            25390, 27286, 30088, 31824, 33877, 34646, 35260, 38616, 41127, 44263,
            46623, 47927, 49189, 50598, 53758, 55103, 56456, 57618, 57777, 59176,
            60078, 62002, 62139, 62332, 62775, 63552, 64040, 65385, 66064, 66863,
            67662, 68034, 68344, 69080, 69356, 69749, 70302, 70947, 71532, 73474,
            74384, 75005, 75489, 76081, 76561, 77448, 78431, 78702, 78976, 79961,
            80507, 82176, 83229, 83835, 85901, 87725, 89678, 90442, 91092, 92259,
            92685, 93542, 94083, 95605, 96071, 97185, 98212, 98840, 100408, 
            101088, 102149, 102590, 105198, 105259, 105672, 106345, 107331, 
            108113, 110302, 110812, 111216, 112610, 112950, 113232, 113495,
            113895, 114584, 115370, 115857, 116723, 117445, 118189, 118873,
            119812, 120187, 120939, 121441, 122570, 122895, 123662, 124211, 
            124907, 125492, 126751, 127260, 127420, 127551, 128294, 129056,
            129686, 129927, 130676, 131693, 132328, 132875, 133239, 133696,
            134506, 135186, 135737, 135902, 136529, 137058, 139370, 139652,
            140562, 141294, 141700, 141993, 142531, 143188, 143813, 144408,
            144561, 145143, 145277, 145850, 146176, 146565, 146877, 147097,
            147195, 147376, 148057, 148682, 148848, 149106, 149484, 149793,
            149916, 150398, 151131, 151313, 151385, 151703, 151811, 151883,
            152020, 152169, 152535, 152819, 152921)
write.table(splits, "CNNC_cluster2_network/BL_gene_pairs_num.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")