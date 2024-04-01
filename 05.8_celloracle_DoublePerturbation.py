###################################################################
#                                                                 #
#             CELLORACLE STEP 8: IN SILICO KO + OE                #
#                                                                 #
###################################################################

# Code created on: 2023/07/20.
# Last modified on: 2023/07/20.

## 1. IMPORT LIBRARIES ############################################
import copy
import glob
import importlib
import time
import os
import shutil
import sys
from importlib import reload

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from tqdm.notebook import tqdm
import statsmodels.stats.multitest as multi

import time

import celloracle as co
from celloracle.applications import Oracle_development_module#, Oracle_systematic_analysis_helper
# from .systematic_analysis_helper import Oracle_systematic_analysis_helper
import systematic_analysis_helper
co.__version__

# Setting up plotting parameters
plt.rcParams["figure.figsize"] = [5,5]
plt.rcParams["savefig.dpi"] = 300


## 2. SET VARIABLES ###############################################
# Data
oracleObject = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct.celloracle.oracle"
GRN_filtered = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_links_processed.celloracle.links"
Pseudotime = "pt_rank_norm"
gradient = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct.celloracle.gradient"
cluster = "seurat_clusters"
injury = "injury"
gene1 = "AMEX60DD044502"
gene2 = "AMEX60DD053152"
value1 = "KO" # Choose from KO or MaxExpression
value2 = "MaxExpression" # Choose from KO or MaxExpression
savename = gene1 + "." + value1 + "_" + gene2 + "." + value2

# Plotting parameters
scale = 10 #12
n_grid = 40
min_mass = 7.8
scale_simulation = 4
scale_dev = 30
n_neighbors = 81
vm = 0.20
cmap = "bwr"

# Simulation parameters and hd5f to save results
n_propagation = 3
file_path = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_"  + gene1 + "-" + value1 + "_" + gene2 + "-" + value2 + "_screen.celloracle.hdf5"
ps_file = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_"  + gene1 + "-" + value1 + "_" + gene2 + "-" + value2 + "_Simulation_Result.tsv"


## 3.  LOAD DATA ##################################################
print(">>> Loading oracle object and filtered GRN")
oracle = co.load_hdf5(oracleObject)
oracle
links = co.load_hdf5(GRN_filtered)
gradient = co.load_hdf5(gradient)

# Make sure that the number of cells and dimensional reduction data is exact same between oracle object and gradient object
print(oracle.adata.shape, gradient.embedding.shape)
assert((oracle.adata.obsm[oracle.embedding_name] == gradient.embedding).all())


## 4. PROCESS GRN TO MAKE PREDICTIVE MODELS FOR SIMULATION #######
print(">>> Fitting models for simulation")
links.filter_links(p=7.671065e-06)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)


## 5. SELECT GENES FOR SIMULATION ################################
# We selected the genes that are DE between BL and CSD at 3dpi and
# have at least 1 target gene in the filtered network.
print(">>> Loading genes for simulation")
genes = ["AMEX60DD044502", "AMEX60DD053152"]

## 6. SYSTEMATIC SIMULATION OF KO ################################
print(">>> Perturbation simulation")
# We will simulate the KOs for the following conditions:
#   1. Using all cells
#   2. Focusing on the BL branch.
#   3. Focusing on the CSD branch.

# We need to get the cell indexes in order to perform simulations 2 & 3
#sc.pl.umap(oracle.adata, color=injury, legend_loc="on data")
#plt.savefig("prueba.png")
#plt.close()
#  List of cluster name
print(sorted(list(oracle.adata.obs[injury].unique())))

# Get cell_id
cell_idx_BL = np.where(oracle.adata.obs[injury].isin(["BL"]))[0]
cell_idx_CSD = np.where(oracle.adata.obs[injury].isin(["CSD"]))[0]
cell_idx_c0 = np.where(oracle.adata.obs[cluster].isin(["c0"]))[0]
cell_idx_c1 = np.where(oracle.adata.obs[cluster].isin(["c1"]))[0]
cell_idx_c2 = np.where(oracle.adata.obs[cluster].isin(["c2"]))[0]
cell_idx_c3 = np.where(oracle.adata.obs[cluster].isin(["c3"]))[0]
cell_idx_c4 = np.where(oracle.adata.obs[cluster].isin(["c4"]))[0]


# Make dictionary to store the cell index list
index_dictionary = {"Whole_cells": None,
                    "BL": cell_idx_BL,
                    "CSD": cell_idx_CSD,
                    "c0": cell_idx_c0,
                    "c1": cell_idx_c1,
                    "c2": cell_idx_c2,
                    "c3": cell_idx_c3,
                    "c4": cell_idx_c4}

# Make a custom function for the systematic KO
# We save all simulations in a single hdf5 file.
# When saving the results, we specify the gene used in the KO and
# the misc (name of the simulation, cell population or condition).

def pipeline(gene_for_OE):
    print("Calculating overexpression value")
    if value1 == "MaxExpression":
        v1 = sc.get.obs_df(oracle.adata, keys=[gene1], layer="imputed_count").max() # get max value
        # Transform value from numpy float 64 to float
        v1 = v1.tolist()
        v1 = v1[0]
    if value1 == "KO":
        v1 = 0.0
    if value2 == "MaxExpression":
        v2 = sc.get.obs_df(oracle.adata, keys=[gene2], layer="imputed_count").max() # get max value
        # Transform value from numpy float 64 to float
        v2 = v2.tolist()
        v2 = v2[0]
    if value2 == "KO":
        v2 = 0.0 

    print("Simulating Perturbation")
    # Simulate perturbation
    oracle.simulate_shift(perturb_condition={gene1: v1, gene2: v2}, ignore_warning=True, n_propagation=3)
    oracle.estimate_transition_prob(n_neighbors=n_neighbors, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    # Do simulation for all conditions.
    for lineage_name, cell_idx in index_dictionary.items():
        # Paths to save the KO plots
        PSgrid = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_OE_"  + gene1 + "-" + value1 + "_" + gene2 + "-" + value2 + "_" + lineage_name + "_PSgrid.pdf"
        PSgridvector = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_OE_" + gene1 + "-" + value1 + "_" + gene2 + "-" + value2 + "_" + lineage_name + "_PSvector.pdf"
        vectorfieldsplot = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_OE_" + gene1 + "-" + value1 + "_" + gene2 + "-" + value2 + "_" + "VectorsFieldPlot.pdf"
        vectorwithcluster = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_OE_" + gene1 + "-" + value1 + "_" + gene2 + "-" + value2 + "_" + "VectorsFieldPlot.pdf"

        dev = Oracle_development_module()
        # Load development flow
        dev.load_differentiation_reference_data(gradient_object=gradient)
        # Load simulation result
        dev.load_perturb_simulation_data(oracle_object=oracle, cell_idx_use=cell_idx, name=lineage_name)
        # Calculate inner product
        dev.calculate_inner_product()
        dev.calculate_digitized_ip(n_bins=10)
        # Visualize perturbation
        ## Plot vector fields (just to check that the scale parameters are alright and we don't see anything on random)
        oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=n_neighbors)
        oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
        fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
        oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
        ax[0].set_title(f"Simulated cell identity shift vector: {gene1} KO and {gene2} overexpression")
        oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
        ax[1].set_title(f"Randomized simulation vector")
        plt.savefig(vectorfieldsplot, format = "pdf")
        plt.close()

        ## Plot pseudoscore plot
        fig, ax = plt.subplots(1, 2, figsize=[12, 6])
        dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0], cmap = cmap)
        ax[0].set_title(f"PS")
        dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1], cmap = cmap)
        ax[1].set_title(f"PS calculated with Randomized simulation vector")
        plt.savefig(PSgrid, format = "pdf")
        plt.close()

        ## Plot pseudoscore with vectors
        fig, ax = plt.subplots(figsize=[6, 6])
        dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax, cmap = cmap)
        dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
        plt.savefig(PSgridvector, format = "pdf")
        plt.close()
        
        # Save results in a hdf5 file.
        dev.set_hdf_path(path=file_path) 
        dev.dump_hdf5(gene=gene_for_OE, misc=lineage_name)

        # Plot vector fields with clusters
        fig, ax = plt.subplots(figsize=[8, 8])
        oracle.plot_cluster_whole(ax=ax, s=10)
        oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
        plt.savefig(vectorwithcluster, format = "pdf")
        plt.close()

# Test pipeline
#pipeline(gene_for_OE="AMEX60DD000244")
#dev_test = Oracle_development_module()
#dev_test.set_hdf_path(path=file_path)
#dev_test.load_hdf5(gene="AMEX60DD000244", misc="Whole_cells")
#dev_test.visualize_development_module_layout_0(s=5, scale_for_simulation=scale_simulation, s_grid=n_grid, scale_for_pseudotime=scale_dev, vm=vm)
#plt.savefig("prueba.png")
#plt.close()

# Run pipeline
for gene in tqdm(genes):
    pipeline(gene_for_OE=gene)


## 7. EXTRACT PS AS A DATAFRAME ###################################
# Load hdf5 with KO results
dev = Oracle_development_module()
dev.set_hdf_path(path=file_path)


# See information of the saved data
info = dev.get_hdf5_info()
print("Genes\n", info["gene_list"])
print("\nSimulation conditions\n", info["misc_list"])


# Test: Check one simulation 
#dev.load_hdf5(gene="AMEX60DD053152", misc="Whole_cells")
#dev.visualize_development_module_layout_0(s=5, 
#                                            scale_for_simulation=scale_simulation,
#                                            s_grid=n_grid,
#                                            scale_for_pseudotime=scale_dev, 
#                                            vm=vm)
#plt.savefig("prueba.png")
#plt.close()                                            


# Get dataframe with perturbation scores (PS)
helper = systematic_analysis_helper.Oracle_systematic_analysis_helper(hdf5_file_path = file_path)

# The script they provided have bugs or something so here i made my own tiny script to extract ps_sum, ps_sum_random and p_value

# Create list with all genes that were KO
loop = list(info["gene_list"])

# Create lists to store values
neg_ps_BL_sums = []
neg_ps_BL_random = []
neg_pvalues_BL = []
pos_ps_BL_sums = []
pos_ps_BL_random = []
pos_pvalues_BL = []

neg_ps_CSD_sums = []
neg_ps_CSD_random = []
neg_pvalues_CSD = []
pos_ps_CSD_sums = []
pos_ps_CSD_random = []
pos_pvalues_CSD = []

neg_ps_sums = []
neg_ps_random = []
neg_pvalues = []
pos_ps_sums = []
pos_ps_random = []
pos_pvalues = []

neg_ps_c0_sums = []
neg_ps_c0_random = []
neg_pvalues_c0 = []
pos_ps_c0_sums = []
pos_ps_c0_random = []
pos_pvalues_c0 = []

neg_ps_c1_sums = []
neg_ps_c1_random = []
neg_pvalues_c1 = []
pos_ps_c1_sums = []
pos_ps_c1_random = []
pos_pvalues_c1 = []

neg_ps_c2_sums = []
neg_ps_c2_random = []
neg_pvalues_c2 = []
pos_ps_c2_sums = []
pos_ps_c2_random = []
pos_pvalues_c2 = []

neg_ps_c3_sums = []
neg_ps_c3_random = []
neg_pvalues_c3 = []
pos_ps_c3_sums = []
pos_ps_c3_random = []
pos_pvalues_c3 = []

neg_ps_c4_sums = []
neg_ps_c4_random = []
neg_pvalues_c4 = []
pos_ps_c4_sums = []
pos_ps_c4_random = []
pos_pvalues_c4 = []

# Iterate through all genes
for gene in loop:
    # BL
    helper.load_hdf5(gene = gene, misc="BL", specify_attributes=["inner_product_df"])
    # Negative scores
    p, ps_sum, ps_ran = helper.get_negative_PS_p_value(return_ps_sum=True)
    neg_ps_BL_sums.append(ps_sum)
    neg_ps_BL_random.append(ps_ran)
    neg_pvalues_BL.append(p)
    #Positive scores
    p, ps_sum, ps_ran = helper.get_positive_PS_p_value(return_ps_sum=True)
    pos_ps_BL_sums.append(ps_sum)
    pos_ps_BL_random.append(ps_ran)
    pos_pvalues_BL.append(p)

    # CSD
    helper.load_hdf5(gene = gene, misc="CSD", specify_attributes=["inner_product_df"])
    p, ps_sum, ps_ran = helper.get_negative_PS_p_value(return_ps_sum=True)
    neg_ps_CSD_sums.append(ps_sum)
    neg_ps_CSD_random.append(ps_ran)
    neg_pvalues_CSD.append(p)
    #Positive scores
    p, ps_sum, ps_ran = helper.get_positive_PS_p_value(return_ps_sum=True)
    pos_ps_CSD_sums.append(ps_sum)
    pos_ps_CSD_random.append(ps_ran)
    pos_pvalues_CSD.append(p)

    # Whole cells
    helper.load_hdf5(gene = gene, misc="Whole_cells", specify_attributes=["inner_product_df"])
    p, ps_sum, ps_ran = helper.get_negative_PS_p_value(return_ps_sum=True)
    neg_ps_sums.append(ps_sum)
    neg_ps_random.append(ps_ran)
    neg_pvalues.append(p)
    #Positive scores
    p, ps_sum, ps_ran = helper.get_positive_PS_p_value(return_ps_sum=True)
    pos_ps_sums.append(ps_sum)
    pos_ps_random.append(ps_ran)
    pos_pvalues.append(p)

    # C0 - blastema mesenchyme
    helper.load_hdf5(gene = gene, misc="c0", specify_attributes=["inner_product_df"])
    p, ps_sum, ps_ran = helper.get_negative_PS_p_value(return_ps_sum=True)
    neg_ps_c0_sums.append(ps_sum)
    neg_ps_c0_random.append(ps_ran)
    neg_pvalues_c0.append(p)
    #Positive scores
    p, ps_sum, ps_ran = helper.get_positive_PS_p_value(return_ps_sum=True)
    pos_ps_c0_sums.append(ps_sum)
    pos_ps_c0_random.append(ps_ran)
    pos_pvalues_c0.append(p)

    # C1 - CSD fibroblasts I
    helper.load_hdf5(gene = gene, misc="c1", specify_attributes=["inner_product_df"])
    p, ps_sum, ps_ran = helper.get_negative_PS_p_value(return_ps_sum=True)
    neg_ps_c1_sums.append(ps_sum)
    neg_ps_c1_random.append(ps_ran)
    neg_pvalues_c1.append(p)
    #Positive scores
    p, ps_sum, ps_ran = helper.get_positive_PS_p_value(return_ps_sum=True)
    pos_ps_c1_sums.append(ps_sum)
    pos_ps_c1_random.append(ps_ran)
    pos_pvalues_c1.append(p)

    # C2 - early injury
    helper.load_hdf5(gene = gene, misc="c2", specify_attributes=["inner_product_df"])
    p, ps_sum, ps_ran = helper.get_negative_PS_p_value(return_ps_sum=True)
    neg_ps_c2_sums.append(ps_sum)
    neg_ps_c2_random.append(ps_ran)
    neg_pvalues_c2.append(p)
    #Positive scores
    p, ps_sum, ps_ran = helper.get_positive_PS_p_value(return_ps_sum=True)
    pos_ps_c2_sums.append(ps_sum)
    pos_ps_c2_random.append(ps_ran)
    pos_pvalues_c2.append(p)

    # C3 - CSD fibroblasts II
    helper.load_hdf5(gene = gene, misc="c3", specify_attributes=["inner_product_df"])
    p, ps_sum, ps_ran = helper.get_negative_PS_p_value(return_ps_sum=True)
    neg_ps_c3_sums.append(ps_sum)
    neg_ps_c3_random.append(ps_ran)
    neg_pvalues_c3.append(p)
    #Positive scores
    p, ps_sum, ps_ran = helper.get_positive_PS_p_value(return_ps_sum=True)
    pos_ps_c3_sums.append(ps_sum)
    pos_ps_c3_random.append(ps_ran)
    pos_pvalues_c3.append(p)

    # C4 - blastema cartilage
    helper.load_hdf5(gene = gene, misc="c4", specify_attributes=["inner_product_df"])
    p, ps_sum, ps_ran = helper.get_negative_PS_p_value(return_ps_sum=True)
    neg_ps_c4_sums.append(ps_sum)
    neg_ps_c4_random.append(ps_ran)
    neg_pvalues_c4.append(p)
    #Positive scores
    p, ps_sum, ps_ran = helper.get_positive_PS_p_value(return_ps_sum=True)
    pos_ps_c4_sums.append(ps_sum)
    pos_ps_c4_random.append(ps_ran)
    pos_pvalues_c4.append(p)

# Create vectors of corrected pvalues
neg_pvalues_BL_Bonferroni = multi.multipletests(neg_pvalues_BL, alpha = 0.05, method = "bonferroni")[1] #[1] to extract vector of adjusted pvalues.
neg_pvalues_BL_FDR = multi.multipletests(neg_pvalues_BL, alpha = 0.05, method = "fdr_bh")[1]
pos_pvalues_BL_Bonferroni = multi.multipletests(pos_pvalues_BL, alpha = 0.05, method = "bonferroni")[1]
pos_pvalues_BL_FDR = multi.multipletests(pos_pvalues_BL, alpha = 0.05, method = "fdr_bh")[1]

neg_pvalues_CSD_Bonferroni = multi.multipletests(neg_pvalues_CSD, alpha = 0.05, method = "bonferroni")[1]
neg_pvalues_CSD_FDR = multi.multipletests(neg_pvalues_CSD, alpha = 0.05, method = "fdr_bh")[1]
pos_pvalues_CSD_Bonferroni = multi.multipletests(pos_pvalues_CSD, alpha = 0.05, method = "bonferroni")[1]
pos_pvalues_CSD_FDR = multi.multipletests(pos_pvalues_CSD, alpha = 0.05, method = "fdr_bh")[1]

neg_pvalues_Bonferroni = multi.multipletests(neg_pvalues, alpha = 0.05, method = "bonferroni")[1]
neg_pvales_FDR = multi.multipletests(neg_pvalues, alpha = 0.05, method = "fdr_bh")[1]
pos_pvalues_Bonferroni = multi.multipletests(pos_pvalues, alpha = 0.05, method = "bonferroni")[1]
pos_pvalues_FDR = multi.multipletests(pos_pvalues, alpha = 0.05, method = "fdr_bh")[1]

neg_pvalues_c0_Bonferroni = multi.multipletests(neg_pvalues_c0, alpha = 0.05, method = "bonferroni")[1]
neg_pvalues_c0_FDR = multi.multipletests(neg_pvalues_c0, alpha = 0.05, method = "fdr_bh")[1]
pos_pvalues_c0_Bonferroni = multi.multipletests(pos_pvalues_c0, alpha = 0.05, method = "bonferroni")[1]
pos_pvalues_c0_FDR = multi.multipletests(pos_pvalues_c0, alpha = 0.05, method = "fdr_bh")[1]

neg_pvalues_c1_Bonferroni = multi.multipletests(neg_pvalues_c1, alpha = 0.05, method = "bonferroni")[1]
neg_pvalues_c1_FDR = multi.multipletests(neg_pvalues_c1, alpha = 0.05, method = "fdr_bh")[1]
pos_pvalues_c1_Bonferroni = multi.multipletests(pos_pvalues_c1, alpha = 0.05, method = "bonferroni")[1]
pos_pvalues_c1_FDR = multi.multipletests(pos_pvalues_c1, alpha = 0.05, method = "fdr_bh")[1]

neg_pvalues_c2_Bonferroni = multi.multipletests(neg_pvalues_c2, alpha = 0.05, method = "bonferroni")[1]
neg_pvalues_c2_FDR = multi.multipletests(neg_pvalues_c2, alpha = 0.05, method = "fdr_bh")[1]
pos_pvalues_c2_Bonferroni = multi.multipletests(pos_pvalues_c2, alpha = 0.05, method = "bonferroni")[1]
pos_pvalues_c2_FDR = multi.multipletests(pos_pvalues_c2, alpha = 0.05, method = "fdr_bh")[1]

neg_pvalues_c3_Bonferroni = multi.multipletests(neg_pvalues_c3, alpha = 0.05, method = "bonferroni")[1]
neg_pvalues_c3_FDR = multi.multipletests(neg_pvalues_c3, alpha = 0.05, method = "fdr_bh")[1]
pos_pvalues_c3_Bonferroni = multi.multipletests(pos_pvalues_c3, alpha = 0.05, method = "bonferroni")[1]
pos_pvalues_c3_FDR = multi.multipletests(pos_pvalues_c3, alpha = 0.05, method = "fdr_bh")[1]

neg_pvalues_c4_Bonferroni = multi.multipletests(neg_pvalues_c4, alpha = 0.05, method = "bonferroni")[1]
neg_pvalues_c4_FDR = multi.multipletests(neg_pvalues_c4, alpha = 0.05, method = "fdr_bh")[1]
pos_pvalues_c4_Bonferroni = multi.multipletests(pos_pvalues_c4, alpha = 0.05, method = "bonferroni")[1]
pos_pvalues_c4_FDR = multi.multipletests(pos_pvalues_c4, alpha = 0.05, method = "fdr_bh")[1]

# Create log vectors
neg_ps_BL_sums_log1p = np.log1p(neg_ps_BL_sums)
neg_ps_BL_random_log1p = np.log1p(neg_ps_BL_random)
pos_ps_BL_sums_log1p = np.log1p(pos_ps_BL_sums)
pos_ps_BL_random_log1p = np.log1p(pos_ps_BL_random)

neg_ps_CSD_sums_log1p = np.log1p(neg_ps_CSD_sums)
neg_ps_CSD_random_log1p = np.log1p(neg_ps_CSD_random)
pos_ps_CSD_sums_log1p = np.log1p(pos_ps_CSD_sums)
pos_ps_CSD_random_log1p = np.log1p(pos_ps_CSD_random)

neg_ps_sums_log1p = np.log1p(neg_ps_sums)
neg_ps_random_log1p = np.log1p(neg_ps_random)
pos_ps_sums_log1p = np.log1p(pos_ps_sums)
pos_ps_random_log1p = np.log1p(pos_ps_random)

neg_ps_c0_sums_log1p = np.log1p(neg_ps_c0_sums)
neg_ps_c0_random_log1p = np.log1p(neg_ps_c0_random)
pos_ps_c0_sums_log1p = np.log1p(pos_ps_c0_sums)
pos_ps_c0_random_log1p = np.log1p(pos_ps_c0_random)

neg_ps_c1_sums_log1p = np.log1p(neg_ps_c1_sums)
neg_ps_c1_random_log1p = np.log1p(neg_ps_c1_random)
pos_ps_c1_sums_log1p = np.log1p(pos_ps_c1_sums)
pos_ps_c1_random_log1p = np.log1p(pos_ps_c1_random)

neg_ps_c2_sums_log1p = np.log1p(neg_ps_c2_sums)
neg_ps_c2_random_log1p = np.log1p(neg_ps_c2_random)
pos_ps_c2_sums_log1p = np.log1p(pos_ps_c2_sums)
pos_ps_c2_random_log1p = np.log1p(pos_ps_c2_random)

neg_ps_c3_sums_log1p = np.log1p(neg_ps_c3_sums)
neg_ps_c3_random_log1p = np.log1p(neg_ps_c3_random)
pos_ps_c3_sums_log1p = np.log1p(pos_ps_c3_sums)
pos_ps_c3_random_log1p = np.log1p(pos_ps_c3_random)

neg_ps_c4_sums_log1p = np.log1p(neg_ps_c4_sums)
neg_ps_c4_random_log1p = np.log1p(neg_ps_c4_random)
pos_ps_c4_sums_log1p = np.log1p(pos_ps_c4_sums)
pos_ps_c4_random_log1p = np.log1p(pos_ps_c4_random)

# Save everything as dataframe
result = pd.DataFrame({"AMEXID": loop,
                        "BL_negative_PS_sum":neg_ps_BL_sums, "BL_negative_PS_random": neg_ps_BL_random, "BL_negative_PS_sum_log1p":neg_ps_BL_sums_log1p, "BL_negative_PS_random_log1p": neg_ps_BL_random_log1p, "BL_negative_pvalue": neg_pvalues_BL, "BL_negative_padj_Bonferroni": neg_pvalues_BL_Bonferroni, "BL_negative_padj_FDR": neg_pvalues_BL_FDR,
                        "BL_positive_PS_sum":pos_ps_BL_sums, "BL_positive_PS_random": pos_ps_BL_random, "BL_positive_PS_sum_log1p":pos_ps_BL_sums_log1p, "BL_positive_PS_random_log1p": pos_ps_BL_random_log1p, "BL_positive_pvalue": pos_pvalues_BL, "BL_positive_padj_Bonferroni": pos_pvalues_BL_Bonferroni, "BL_positive_padj_FDR": pos_pvalues_BL_FDR,
                        "CSD_negative_PS_sum":neg_ps_CSD_sums, "CSD_negative_PS_random": neg_ps_CSD_random, "CSD_negative_PS_sum_log1p":neg_ps_CSD_sums_log1p, "CSD_negative_PS_random_log1p": neg_ps_CSD_random_log1p, "CSD_negative_pvalue": neg_pvalues_CSD, "CSD_negative_padj_Bonferroni": neg_pvalues_CSD_Bonferroni, "CSD_negative_padj_FDR": neg_pvalues_CSD_FDR,
                        "CSD_positive_PS_sum":pos_ps_CSD_sums, "CSD_positive_PS_random": pos_ps_CSD_random, "CSD_positive_PS_sum_log1p":pos_ps_CSD_sums_log1p, "CSD_positive_PS_random_log1p": pos_ps_CSD_random_log1p, "CSD_positive_pvalue": pos_pvalues_CSD, "CSD_positive_padj_Bonferroni": pos_pvalues_CSD_Bonferroni, "CSD_positive_padj_FDR": pos_pvalues_CSD_FDR,
                        "AllCells_negative_PS_sum":neg_ps_sums, "AllCells_negative_PS_random": neg_ps_random, "AllCells_negative_PS_sum_log1p":neg_ps_sums_log1p, "AllCells_positive_PS_random_log1p": neg_ps_random_log1p, "AllCells_negative_pvalue": neg_pvalues, "AllCells_negative_padj_Bonferroni": neg_pvalues_Bonferroni, "AllCells_negative_padj_FDR": neg_pvalues_Bonferroni,
                        "AllCells_positive_PS_sum":pos_ps_sums, "AllCells_positive_PS_random": pos_ps_random, "AllCells_positive_PS_sum_log1p":pos_ps_sums_log1p, "AllCells_positive_PS_random_log1p": pos_ps_random_log1p, "AllCells_positive_pvalue": pos_pvalues, "AllCells_positive_padj_Bonferroni": pos_pvalues_Bonferroni, "AllCells_positive_padj_FDR": pos_pvalues_FDR,
                        "c0_negative_PS_sum":neg_ps_c0_sums, "c0_negative_PS_random": neg_ps_c0_random, "c0_negative_PS_sum_log1p":neg_ps_c0_sums_log1p, "c0_negative_PS_random_log1p": neg_ps_c0_random_log1p, "c0_negative_pvalue": neg_pvalues_c0, "c0_negative_padj_Bonferroni": neg_pvalues_c0_Bonferroni, "c0_negative_padj_FDR": neg_pvalues_c0_FDR,
                        "c0_positive_PS_sum":pos_ps_c0_sums, "c0_positive_PS_random": pos_ps_c0_random, "c0_positive_PS_sum_log1p":pos_ps_c0_sums_log1p, "c0_positive_PS_random_log1p": pos_ps_c0_random_log1p, "c0_positive_pvalue": pos_pvalues_c0, "c0_positive_padj_Bonferroni": pos_pvalues_c0_Bonferroni, "c0_positive_padj_FDR": pos_pvalues_c0_FDR,
                        "c1_negative_PS_sum":neg_ps_c1_sums, "c1_negative_PS_random": neg_ps_c1_random, "c1_negative_PS_sum_log1p":neg_ps_c1_sums_log1p, "c1_negative_PS_random_log1p": neg_ps_c1_random_log1p, "c1_negative_pvalue": neg_pvalues_c1, "c1_negative_padj_Bonferroni": neg_pvalues_c1_Bonferroni, "c1_negative_padj_FDR": neg_pvalues_c1_FDR,
                        "c1_positive_PS_sum":pos_ps_c1_sums, "c1_positive_PS_random": pos_ps_c1_random, "c1_positive_PS_sum_log1p":pos_ps_c1_sums_log1p, "c1_positive_PS_random_log1p": pos_ps_c1_random_log1p, "c1_positive_pvalue": pos_pvalues_c1, "c1_positive_padj_Bonferroni": pos_pvalues_c1_Bonferroni, "c1_positive_padj_FDR": pos_pvalues_c1_FDR,
                        "c2_negative_PS_sum":neg_ps_c2_sums, "c2_negative_PS_random": neg_ps_c2_random, "c2_negative_PS_sum_log1p":neg_ps_c2_sums_log1p, "c2_negative_PS_random_log1p": neg_ps_c2_random_log1p, "c2_negative_pvalue": neg_pvalues_c2, "c2_negative_padj_Bonferroni": neg_pvalues_c2_Bonferroni, "c2_negative_padj_FDR": neg_pvalues_c2_FDR,
                        "c2_positive_PS_sum":pos_ps_c2_sums, "c2_positive_PS_random": pos_ps_c2_random, "c2_positive_PS_sum_log1p":pos_ps_c2_sums_log1p, "c2_positive_PS_random_log1p": pos_ps_c2_random_log1p, "c2_positive_pvalue": pos_pvalues_c2, "c2_positive_padj_Bonferroni": pos_pvalues_c2_Bonferroni, "c2_positive_padj_FDR": pos_pvalues_c2_FDR,
                        "c3_negative_PS_sum":neg_ps_c3_sums, "c3_negative_PS_random": neg_ps_c3_random, "c3_negative_PS_sum_log1p":neg_ps_c3_sums_log1p, "c3_negative_PS_random_log1p": neg_ps_c3_random_log1p, "c3_negative_pvalue": neg_pvalues_c3, "c3_negative_padj_Bonferroni": neg_pvalues_c3_Bonferroni, "c3_negative_padj_FDR": neg_pvalues_c3_FDR,
                        "c3_positive_PS_sum":pos_ps_c3_sums, "c3_positive_PS_random": pos_ps_c3_random, "c3_positive_PS_sum_log1p":pos_ps_c3_sums_log1p, "c3_positive_PS_random_log1p": pos_ps_c3_random_log1p, "c3_positive_pvalue": pos_pvalues_c3, "c3_positive_padj_Bonferroni": pos_pvalues_c3_Bonferroni, "c3_positive_padj_FDR": pos_pvalues_c3_FDR,
                        "c4_negative_PS_sum":neg_ps_c4_sums, "c4_negative_PS_random": neg_ps_c4_random, "c4_negative_PS_sum_log1p":neg_ps_c4_sums_log1p, "c4_negative_PS_random_log1p": neg_ps_c4_random_log1p, "c4_negative_pvalue": neg_pvalues_c4, "c4_negative_padj_Bonferroni": neg_pvalues_c4_Bonferroni, "c4_negative_padj_FDR": neg_pvalues_c4_FDR,
                        "c4_positive_PS_sum":pos_ps_c4_sums, "c4_positive_PS_random": pos_ps_c4_random, "c4_positive_PS_sum_log1p":pos_ps_c4_sums_log1p, "c4_positive_PS_random_log1p": pos_ps_c4_random_log1p, "c4_positive_pvalue": pos_pvalues_c4, "c4_positive_padj_Bonferroni": pos_pvalues_c4_Bonferroni, "c4_positive_padj_FDR": pos_pvalues_c4_FDR})


# Save all PS in a single data frame
result.to_csv(ps_file,   sep = "\t", header = True, index = False)
print("DONE.")
exit()