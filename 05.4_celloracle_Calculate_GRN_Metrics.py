#!/usr/bin/env python3

###################################################################
#                                                                 #
#               CELLORACLE STEP 4: NETWORK ANALYSIS               #
#                                                                 #
###################################################################

# Code created on: 2023/07/12.
# Last modified on: 2023/07/18.

## 1. IMPORT LIBRARIES ############################################
import os
import sys
from tkinter.tix import Meter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import datetime
import celloracle as co
co.__version__

print("Running script celloracle_4_NetworkAnalysis.py")
e = datetime.datetime.now()
print ("Today's date:  = %s/%s/%s" % (e.day, e.month, e.year))
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))


## 2. SET VARIABLES ###############################################
scorePlots = '/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_ranked_score_plots'
scoreComparison = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_score_comparison_plots"
networkScorePerGene = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_network_score_per_gene_plots"
NetworkScoreDistribution = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_network_score_distribution_plots"
GRN_processed = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_links_processed.celloracle.links"
cluster_colors = ["#56b4e9", "#cc79a7", "#6a3d9a", "orange", "#0072b2", "#b15928", "#ff7f00"]


print(">>> Setting variables")
print("Path to GRN object: ", GRN_processed)
print("The Genes with highest scores will be plotted and saved at: ",  scorePlots)
print("The score comparison between networks will be saved at: ", scoreComparison)
print("The network score per gene will be saved at: ", networkScorePerGene)
print("The network score distribution will be saved at: ", NetworkScoreDistribution)

## 3. LOAD NETWORK ################################################
print(">>> Loading network")
links = co.load_hdf5(file_path=GRN_processed)


## 4. PLOT TOP GENES FOR EACH METRIC IN EACH CLUSTER ##############
print(">>> Plotting top genes for each cluster")
for c in links.cluster:
    print("Plotting top30 genes for cluster ", c)
    plt.figure(figsize = (10, 6))
    links.plot_scores_as_rank(cluster=c, n_gene=30, save=scorePlots)
    plt.close()

## 5. PLOTTING A GENE's NETWORK SCORE DURING "DIFFERENTIATION" #####

print(">>> Plotting gene's network score changes during differentiation")
df = pd.concat(links.filtered_links, axis=0)
genes1=list(df["source"])
genes2=list(df["target"])
genes = genes1 + genes2
genes=list(set(genes))
i = 1
for g in genes:
    print(i,"/",len(genes))
    plt.figure(figsize=(10,6))
    links.plot_score_per_cluster(goi=g, 
                                save=networkScorePerGene)
    plt.close()
    i= i + 1


## 6. PLOTTING NETWORK SCORE DISTRIBUTIONS ########################
print(">>> Plotting score distributions")
metrics = ["degree_centrality_all", "degree_centrality_in", 
            "degree_centrality_out", "betweenness_centrality", 
            "eigenvector_centrality"]
for m in metrics:
    plt.figure(figsize = (6, 4.5))
    maxval = round(max(links.merged_score[m]),1)
    plt.ylim([0, maxval])
    links.plot_score_discributions(values=[m],
                               method="boxplot",
                               save=NetworkScoreDistribution)
    plt.close()


## 7. SAVING FILES ################################################
# Save list of genes for further processing
genesDF = pd.DataFrame(genes, columns = ["AMEXID"])
genesDF.to_csv("/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_GenesInGRN.tsv", 
                header = True, index = False, sep = "\t")

# Save df with all metrics
links.merged_score.to_csv("/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_NetworkScore_AllClusters.tsv",
                            header = True, index = True, sep = "\t")

# Save the GRN of each cluster
links.filtered_links["c0"].to_csv("/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_GRN_c0.tsv",
                                    header = True, index = False, sep = "\t")
links.filtered_links["c1"].to_csv("/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_GRN_c1.tsv",
                                    header = True, index = False, sep = "\t")
links.filtered_links["c2"].to_csv("/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_GRN_c2.tsv",
                                    header = True, index = False, sep = "\t")
links.filtered_links["c3"].to_csv("/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_GRN_c3.tsv",
                                    header = True, index = False, sep = "\t")
links.filtered_links["c4"].to_csv("/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_GRN_c4.tsv",
                                    header = True, index = False, sep = "\t")
links.filtered_links["c5"].to_csv("/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_GRN_c5.tsv",
                                    header = True, index = False, sep = "\t")
links.filtered_links["c6"].to_csv("/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_GRN_c6.tsv",
                                    header = True, index = False, sep = "\t")                                                                                                                                                                                                                        
e = datetime.datetime.now()
print ("Today's date:  = %s/%s/%s" % (e.day, e.month, e.year))
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))
print("DONE.")
exit()