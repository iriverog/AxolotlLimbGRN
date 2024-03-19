#!/usr/bin/env python3

###################################################################
#                                                                 #
#               CELLORACLE STEP 3: NETWORK PROCESSING             #
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

print("Running script celloracle_3_NetworkProcessing.py")
e = datetime.datetime.now()
print ("Today's date:  = %s/%s/%s" % (e.day, e.month, e.year))
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))


## 2. SET VARIABLES ###############################################
GRN = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_links.celloracle.links"
degreedistribution = '/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_degree_distribution.plots'
GRN_processed = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_links_processed.celloracle.links"

print(">>> Setting variables")
print("Path to GRN object: ", GRN)
print("The processed GRN object will be saved at: ", GRN_processed)
print("The degree distribution plots will be saved at: ", degreedistribution)

## 3. LOAD NETWORK ################################################
print(">>> Loading network")
links = co.load_hdf5(file_path=GRN)


## 4. FILTERING WEAK EDGES ########################################
print(">>> Filtering weak edges")
links.filter_links(p=7.671065e-06) # weight="coef_abs", threshold_number=5000)


## 5. PLOT DEGREE DISTRIBUTIONS ###################################
print(">>> Plotting degree distribution")
plt.rcParams["figure.figsize"] = [9, 4.5]
links.plot_degree_distributions(plot_model=True,
                                save=degreedistribution)


## 6. CALCULATING NETWORK METRICS #################################
print(">>> Calculating network metrics")
links.get_network_score()


## 7. SAVE PROCESSED GRN FOR IN SILICO PERTURBATION ANALYSIS ######
print(">>> Save processed GRN")
links.to_hdf5(file_path=GRN_processed)

e = datetime.datetime.now()
print ("Today's date:  = %s/%s/%s" % (e.day, e.month, e.year))
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))
print("DONE.")
exit()