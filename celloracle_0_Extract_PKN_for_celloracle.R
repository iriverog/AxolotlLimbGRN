################################################################################
#                                                                              #
#                        EXTRACT GRN TO RUN CELLORACLE                         #
#                                                                              #
################################################################################

# Load TENET network (merged) and prioiritization of candidates
setwd("/data_lab_MT/Ines/Tanaka/BoneRegeneration/UpdateBatch2/")
net <- read.csv("20230704_BLct_CSDct_edgeDF_TENET_ForCytoscape.tsv",
                header = TRUE, sep = "\t")
genes <- read.csv("20230704_BLct_CSDct_CandidateSelection_TENET_Annotated_AllGenes_WithNetworkInfo_ForCytoscape.tsv",
                  header = TRUE, sep = "\t")

# Network filtering 1: find the genes with a |z-score| > 1 that are also differentially expressed
candidates <- genes[abs(genes$CompareRanks_Zscore) > 1 & genes$p_val_adj < 0.05, "AMEXID"]

# Network filtering 2: get all genes that are TFs. Keep only candidates that are TFs.
TFs <- genes[genes$IsTF == TRUE, "AMEXID"]
candidates <- intersect(candidates, TFs)

# Network filtering 3: Expand the gene set to include 1st neighbors
vecinos <- c(net[net$Source_AMEXID %in% candidates, "Target_AMEXID"],
             net[net$Target_AMEXID %in% intersect(candidates, TFs), "Source_AMEXID"],
             candidates)
vecinos <- unique(vecinos)

# Network filtering 4: Select interactions candidate -> candidate, candidate -> neighbor,
# neighbor -> candidate and neighbor -> neighbor
net2 <- net[net$Source_AMEXID %in% vecinos & net$Target_AMEXID %in% vecinos,]
genes2 <- genes[genes$AMEXID %in% vecinos,]
# Prepare data for cell oracle
