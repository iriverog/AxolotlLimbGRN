###################################################################
#                                                                 #
#        CELLORACLE STEP 6: CALCULATE PSEUDOTIME GRADIENT         #
#                                                                 #
###################################################################

# Code created on: 2023/07/12.
# Last modified on: 2023/07/20.

## 1. IMPORT LIBRARIES ############################################
from ast import expr
from cmath import exp
from genericpath import samefile
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

print("Running script 20230718_CT_celloracle_6.py")
e = datetime.datetime.now()
print ("Today's date:  = %s/%s/%s" % (e.day, e.month, e.year))
print ("Time: = %s:%s:%s" % (e.hour, e.minute, e.second))


## 2. SET VARIABLES ###############################################
oracleObject = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct.celloracle.oracle"
GRN_filtered = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_links_processed.celloracle.links"
Pseudotime = "pt_rank_norm"
goi =  "AMEX60DD044502"
symbol = "LEF1"
scale = 10
n_grid = 40
min_mass = 7.8
scale_simulation = 4
scale_dev = 40
n_neighbors = 81
vm = 0.2
cmap = "bwr"
umapgoi = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_UMAP_" + goi + "_" + symbol + ".pdf"
expressionhistogram = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_ExpressionHistogram_" + goi + "_" + symbol + ".pdf"
sanitycheck_plots = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_SanityCheck_" + goi + "_" + symbol + ".pdf"
scalesetup = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_ScaleSetUp" + str(scale) +"_" + goi + "_" + symbol + ".pdf"
minmasssuggestion = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_MinMassSuggestion_.pdf"
chosenminmass = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_ChosenMinMass_.pdf"
vectorfieldsplot = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_VectorFieldsPlot_" + goi + "_" + symbol + ".pdf"
vectorwithcluster = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_VectorsWithClusters_" + goi + "_" + symbol + ".pdf"
pseudotimeplot = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_pseudotimeUMAP.pdf"
pseudotimegrid = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_pseudotimeGrid.pdf"
pseudotimegradient = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_pseudotimeGradient.pdf"
gradient_file = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct.celloracle.gradient"
PSgrid = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_KO_" + goi + "_" + symbol + "_PSgrid.pdf"
PSgridvector = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_KO_" + goi + "_" + symbol + "_PSvector.pdf"
PSmodule = "/home/iriverog/NetVolumes/LAB_MT/RESULTADOS/Ines/Axolotl_Limb_Regeneration/Results/celloracle/BLct_CSDct_KO_" + goi + "_" + symbol + "_PSresult.pdf"


print(">>> Setting variables")
print("Path to oracle object: ", oracleObject)
print("Path to filtered GRN: ", GRN_filtered)
print("The gene of interest is: ", goi)
print("Umap plots of GOI will be saved at: ", umapgoi)
print("Histogram of GOI expression will be saved at: ", expressionhistogram)
print("Sanity check simulation plots will be saved at: ", sanitycheck_plots)
print("Variable storing pseudotime information: ", Pseudotime)
print("Scale: ", scale)

## 3.  LOAD ORACLE OBJECT AND FILTERED GRN ########################
print(">>> Loading oracle object and filtered GRN")
oracle = co.load_hdf5(oracleObject)
oracle
links = co.load_hdf5(GRN_filtered)


## 4. FIT BADGE RIDGE REGRESSION MODELS FOR PREDICTION ############
print(">>> Fitting badge regression models for prediction")
#links.filter_links(p=0.001, threshold_number = 5000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10,
                              use_cluster_specific_TFdict=True)
oracle


## 5. IN SILICO PERTURBATION ANALYSIS ##############################
print(">>> Plotting expression of ", goi, " (", symbol, ")")
plt.rcParams["figure.figsize"] = [6,6]
plt.rcParams["savefig.dpi"] = 300
sc.pl.umap(oracle.adata, 
                    color=[goi, oracle.cluster_column_name],
                    layer="imputed_count", 
                    use_raw=False, 
                    cmap="viridis")
plt.savefig(umapgoi, format = "pdf")
plt.close()

sc.get.obs_df(oracle.adata, 
                keys=[goi], 
                layer="imputed_count").hist()
plt.savefig(expressionhistogram, format = "pdf")
plt.close()


# Enter perturbation conditions to simulate signal propagation after the perturbation.
# We can enter any value but they recommend to avoid extremely
# high values far from the natural gene expression range.
# The upper limit allowed is twice the max. gene expression.
print(">>> Simulating perturbation")
oracle.simulate_shift(perturb_condition={goi: 0.0},
                      n_propagation=3)


## 6. EVALUATION OF THE DISTRIBUTION PATTERN OF SIMULATED VALUES
# The simulated values should be in the same range as real values.
# If they deviate too much, celloracle raises automatic warnings.
# Still, we evaluate the distribution to be super sure.

# Histogram op top10 genes with largest difference vs real values.
print(">>> Sanity check of simulated values")
oracle.evaluate_and_plot_simulation_value_distribution(n_genes=10,
                                                       save=sanitycheck_plots)
plt.close()

# There's a clipping optional step to fix things but I shouldn't
# need to use it.


## 7. CALCULATION OF TRANSITION PROBABILITIES BETWEEN CELL STATES
# Get transition probability
print(">>> Calculating cell state transition probabilities")
oracle.estimate_transition_prob(n_neighbors=n_neighbors,
                                knn_random=True,
                                sampled_fraction=1)

# Calculate embedding
oracle.calculate_embedding_shift(sigma_corr=0.05)


## 8. VISUALIZATION ##########################################
# CAUTION! Is important to find the optimal scale parameters
print(">>> Visualization")

# Finding the optimal scale parameters based on the data
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
# Show quiver plot
oracle.plot_quiver(scale=scale, ax=ax[0]) # Error!
ax[0].set_title(f"Simulated cell identity shift vector: {goi} {symbol} KO")
# Show quiver plot that was calculated with randomized graph.
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
plt.savefig(scalesetup, format = "pdf")
plt.close()

# Find parameters for n_grid and min_mass
oracle.calculate_p_mass(smooth=0.8, 
                       n_grid=n_grid, 
                        n_neighbors=n_neighbors)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig(minmasssuggestion)
plt.close()

oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig(chosenminmass)
plt.close()

# Plot vector fields
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
# Show quiver plot
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, 
                                    ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} {symbol} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, 
                                            ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
plt.savefig(vectorfieldsplot)
plt.close()

# Plot vector fields with cell cluster
fig, ax = plt.subplots(figsize=[8, 8])

oracle.plot_cluster_whole(ax=ax, s=10)
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, 
                                    ax=ax, 
                                    show_background=False)
plt.savefig(vectorwithcluster)
plt.close()


## 10. COMPARISON OF SIMULATION VECTORS WITH DEVELOPMENTAL VECTORS
print(">>> Comparing simulation vetors with developmental vectors")
# We have simulated how TF perturbations affect cell identity
# We can compare the simulated perturbations with devo gradient vectors
# This helps us understand how TF impact dell fate determination

# First: transfer the pseudotime data into a nxn grid
# Second: calculate the 2D gradient of pseudotime to get vector field
# Third: compare the perturbation vector field with the pseudotime vector field (inner product)

print(">>> Visualizing pseudotime")
fig, ax = plt.subplots(figsize=[6,6])
sc.pl.embedding(adata=oracle.adata, 
                    basis=oracle.embedding_name, 
                    ax=ax, 
                    cmap="rainbow",
                    color=[Pseudotime])
plt.savefig(pseudotimeplot)
plt.close()

print(">>> Transforming pseudotime to a grid")
from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, 
                                pseudotime_key=Pseudotime)
gradient.calculate_p_mass(smooth=0.8, 
                            n_grid=n_grid, 
                            n_neighbors=n_neighbors)
gradient.calculate_mass_filter(min_mass=min_mass, 
                                plot=True)
plt.savefig(chosenminmass)
plt.close()

gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, 
                                    plot=True)
plt.savefig(pseudotimegrid)
plt.close()

print(">>> Calculating the 2D vector map to represent the pseudotime gradient")
gradient.calculate_gradient()
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig(pseudotimegradient)
plt.close()

# Just to plot pseudotime gradient vectors
#fig, ax = plt.subplots(figsize=[6, 6])
#gradient.plot_dev_flow_on_grid(scale=scale_dev, ax=ax)
#plt.savefig("pruebaCT13.png")
#plt.close()

# Save gradient object if you want.
gradient.to_hdf5(gradient_file)

# Inner product to compare the 2 vectors
# The inner product will be a positive value when two vectors
#  are pointing in the same direction. The inner product will 
# be a negative value when two vectors are pointing in the 
# opposing directions.
# Example: A -> & B -> makes inner product = 1
# Example: A -> & B <- makes inner product = -1
# Example: if A and B are perpendicular inner product = 0
# The length of the vectors also affect the magnitude of the
# inner product
# Inner product gets saved as a perturbation score (PS)
# PS > 0 means that the TF perturbation promotes differentiation
# PS < 0 means that the TF perurbation blocks differentiation

print(">>> Calculating inner product")
from celloracle.applications import Oracle_development_module
# Make Oracle_development_module to compare two vector field
dev = Oracle_development_module()

# Load development flow
dev.load_differentiation_reference_data(gradient_object=gradient)

# Load simulation result
dev.load_perturb_simulation_data(oracle_object=oracle)


# Calculate inner produc scores
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

# Perturbation score visualization
# Need to adjust the vm parameter for the PS color visualization
# If no color in real plot, make vm smaller.
# If color in random plot, make vm bigger.
# Show perturbation scores

fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0], cmap = cmap)
ax[0].set_title(f"PS")

dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1], cmap = cmap)
ax[1].set_title(f"PS calculated with Randomized simulation vector")
plt.savefig(PSgrid, format = "pdf")
plt.close()

# Show perturbation scores with perturbation simulation vector field
fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax, cmap = cmap)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, 
                                    show_background=False, 
                                    ax=ax)
plt.savefig(PSgridvector)
plt.close()


#dev.visualize_development_module_layout_0(s=5,
#                                          scale_for_simulation=scale_simulation,
#                                          s_grid=50,
#                                          scale_for_pseudotime=scale_dev,
#                                          vm=vm)
#plt.savefig(PSmodule)
#plt.close()



