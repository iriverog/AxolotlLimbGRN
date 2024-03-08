######################################################################
#                                                                    #
# SAVE CNNC PREDICTIONS AS A CLEAN DF FOR CYTOSCAPE                  #
#                                                                    #
######################################################################

## Author: Ines Rivero-Garcia, based on "https://github.com/xiaoyeye/CNNC/blob/master/predict_no_y.py"
## Date: 25/01/2023
## This script takes the output from CNNC (predictions) as well as a gene pair list, cleans results 
## and saves them as a data frame (.tsv) for easier access.

## 1. LOAD LIBRARIES #################################################################################
import os,sys
import numpy as np
import pandas as pd 


## 2. SAVE COMMAND LINE ARGUMENTS AS VARIABLES #######################################################
# pred = numpy array as outputed by CNNC
# Example: pred = "CNNC_cluster2_network/25012023_BL3-5_CanonicalWnt/y_predict.npy"
pred = sys.argv[1]

# df = df with gene pairs inputed to CNNC
# Example: df = "CNNC_cluster2_network/CanonicalWnt_gene_pairs.txt"
df = sys.argv[2]

# filename = filepath where the clean DF with results should be saved
# Example: filename = 'CNNC_cluster2_network/24012023_BL8-11_2labels_TestPredictionDF.tsv'
filename = sys.argv[3]


## 3. LOAD PREDICTION RESULTS AND GENE PAIRS DATAFRAME ###############################################
pred = np.load(file = pred)
df = pd.read_csv(df, sep = '\t', header = None)


## 4. BINARIZE PREDICTION RESULTS (0 if prediction < 0.5, 1 otherwise) ###############################
prob_list = list(pred)
print('Number of interactions that were inputed for prediction: ' + str(len(prob_list)))
pred_list = []

# Transform the probs into a binary list
for i in range(0, len(prob_list)):
    if prob_list[i][0] < 0.5:
        pred_list.append(0)
    else:
        pred_list.append(1)


## 4. FIX DATAFRAME COLNAMES ##########################################################################
df.columns = ['From', 'To']


## 5. CHECK THAT WE AREN'T MISSING ANY GENES ##########################################################
sanitycheck = len(pred_list) == len(prob_list) == df.shape[0]

if sanitycheck == True:
    print('Saving predictions as ' + filename)
else:
    print('The number of predictions does not match the number of gene pairs')
    print('Number of predictions: ' + len(pred_list))
    print('Number of gene pairs: ' + len(df))
    sys.exit()

# If we passed the checks, add the predictions (both probabilities and discrete labels) to output file
df['Pred_Label'] = pred_list
df['Prob_Label'] = [item for sublist in prob_list for item in sublist]
df.to_csv(filename, sep = '\t', index=False)