#!/usr/bin/env python3

###################################################################
#                                                                 #
#    CELLORACLE STEP 1: LOAD DATA AND CREATE CELLORACLE OBJECT    #
#                                                                 #
###################################################################

# Code created on: 2023/07/11.
# Last modified on: 2023/07/18.

## 1. IMPORT LIBRARIES ############################################
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import datetime
import celloracle as co
co.__version__

print("Running script celloracle_1_CreateOracleObject.py")
e = datetime.datetime.now()
print ("Today's date:  = %s/%s/%s" % (e.day, e.month, e.year))
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))


## 2. SET VARIABLES ###############################################
scdata = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_filter_Harmony_v2.h5ad"
GRN = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/20230704_BLct_CSDct_edgeDF_TENET_celloracle.tsv"
cluster = "seurat_clusters"
reduction = "X_umap"
save_oracle = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct.celloracle.oracle"
pcplot = '/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_SelectImportantPCs.pdf'
umapplot = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_UMAP.pdf"
cluster_colors = ["#56b4e9", "#cc79a7", "#6a3d9a", "orange", "#0072b2", "#b15928", "#ff7f00"]

print(">>> Setting variables")
print("Path to scRNA-seq data: ", scdata)
print("Path to GRN: ", GRN)
print("Cluster information: ", cluster)
print("Dimensionality reduction used: ", reduction)
print("The oracle object will be saved at: ", save_oracle)


## 3. LOAD scRNA-SEQ DATA AS SCANPY OBJECT ########################
print(">>> Loading scRNA-seq data")
adata = sc.read_h5ad(scdata)
print(adata)
print(f"Cell number is :{adata.shape[0]}")
print(f"Gene number is :{adata.shape[1]}") 

## 4. INITIALITE ORACLE OBJECT ####################################
print(">>> Initializing oracle object and checking scRNAseq metadata")
oracle = co.Oracle()

# check anndata
print("Metadata columns :", list(adata.obs.columns))
print("Dimensional reduction: ", list(adata.obsm.keys()))

# In this notebook, we use the unscaled mRNA count for the input of 
# Oracle object.
adata.X = adata.raw.X

# Instantiate Oracle object.
# First we need to save seurat clusters as categorical.
# Linea para cambiar los colores y que coincidan con los de R??
plt.rcParams["figure.figsize"] = [8, 6]
sc.pl.umap(adata, color="seurat_clusters", palette = cluster_colors)
plt.savefig(umapplot, format = "pdf")
plt.close()
oracle.import_anndata_as_raw_count(adata=adata, 
    cluster_column_name=cluster, 
    embedding_name=reduction)


## 5. LOAD GRN INFORMATION ########################################
# Load the TF and target gene information.
print(">>> Loading GRN and re-formatting for cell oracle.")
GRN_data = pd.read_csv(GRN,sep = "\t")
# Make dictionary: dictionary key is TF and dictionary value is list of target genes.
TF_to_TG_dictionary = {}

for TF, TGs in zip(GRN_data.TF, GRN_data.Target_genes):
    # convert target gene to list
    TG_list = TGs.replace(" ", "").split(",")
    # store target gene list in a dictionary
    TF_to_TG_dictionary[TF] = TG_list

# We invert the dictionary above using a utility function in celloracle.
TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)

# Add GRN to oracle object
oracle.import_TF_data(TFdict=TG_to_TF_dictionary)


## 6. IMPUTE KNN TO CALCULATE TRANSITIONS #########################
# Perform PCA
print(">>> Imputing KNN to calculate transitions")
oracle.perform_PCA()

# Select important PCs
# These next lines plot the umap but we don't need it.
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.savefig(pcplot, format = "pdf")
plt.close()

print(n_comps)
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, 
                        k=k, 
                        balanced=True, 
                        b_sight=k*8,
                        b_maxl=k*4,
                        n_jobs=4)

## 6. CHECK ORACLE OBJECT #########################################
print(">>> Checking and saving oracle object")
print(oracle)

## 7. SAVE ORACLE OBJECT AND EXIT #################################
# Save oracle object.
oracle.to_hdf5(save_oracle)

e = datetime.datetime.now()
print ("Today's date:  = %s/%s/%s" % (e.day, e.month, e.year))
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))
print("DONE.")
exit()