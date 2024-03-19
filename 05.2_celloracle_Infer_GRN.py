#!/usr/bin/env python3

###################################################################
#                                                                 #
#                   CELLORACLE STEP 2: INFER GRN                  #
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

print("Running script celloracle_2_InferGRN.py")
e = datetime.datetime.now()
print ("Today's date:  = %s/%s/%s" % (e.day, e.month, e.year))
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))


## 2. SET VARIABLES ###############################################
oracleObject = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct.celloracle.oracle"
cluster = "seurat_clusters"
reduction = "X_umap"
save_GRN = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_links.celloracle.links"

print(">>> Setting variables")
print("Path to oracle object: ", oracleObject)
print("Cluster information: ", cluster)
print("Dimensionality reduction used: ", reduction)
print("The GRN object will be saved at: ", save_GRN)


## 2. LOAD ORACLE OBJECT ##########################################
print(">>> Loading oracle object and plotting reduction.")
oracle = co.load_hdf5(oracleObject)
oracle


## 3. INFER GRN ###################################################
print(">>> Inferring GRN")
e = datetime.datetime.now()
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))
links = oracle.get_links(cluster_name_for_GRN_unit=cluster, 
                            alpha=10,
                            verbose_level=10)
e = datetime.datetime.now()
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))


## 4. SAVE GRN FILE ###############################################
links.to_hdf5(file_path=save_GRN)

e = datetime.datetime.now()
print ("Today's date:  = %s/%s/%s" % (e.day, e.month, e.year))
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))
print("DONE.")
exit()