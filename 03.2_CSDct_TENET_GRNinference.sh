#!/bin/bash

# Date of creation: 2023/05/19.
# Date of last modification: 2023/07/02.

echo ">>> Running TENET" || exit 100 &&
echo "source ~/miniconda3.cluster.source && conda activate TENET" || exit 100 && 
source ~/miniconda3.cluster.source && conda activate TENET || exit 100 &&

# CSD
echo "RUNNING TENET ON CSD CONNECTIVE TISSUE" || exit 100 &&
date +"%F %X" || exit 100 &&
echo "./TENET [expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length] [output_file]" || exit 100 &&
echo "./TENET 20230526_CSDct_TENET_ExpressionMatrix.csv 11 20230616_CSDct_TENET_RealTimePoints.txt 20230526_CSDct_TENET_SelectedCells.txt 1 20230616_CSDct_TENET_RealTime_TEmatrix.txt" || exit 100 &&
./TENET 20230526_CSDct_TENET_ExpressionMatrix.csv 11 20230616_CSDct_TENET_RealTimePoints.txt 20230526_CSDct_TENET_SelectedCells.txt 1 20230616_CSDct_TENET_RealTime_TEmatrix.txt || exit 100 &&
echo "DONE WITH CSD CT NETWORK" || exit 100 &&
date +"%F %X" || exit 100 &&

# mv TE_result_matrix.txt 20230526_CSDct_TENET_TEmatrix.txt

# Make GRN
echo "RECONSTRUCTING GRN ON CSD CONNECTIVE TISSUE" || exit 100 &&
date +"%F %X" || exit 100 &&
echo "python makeGRN.py [TE_matrix_file] [cutoff for FDR]" || exit 100 &&
echo "python makeGRN.py 20230616_CSDct_TENET_RealTime_TEmatrix.txt 0.05" || exit 100 &&
python makeGRN.py 20230616_CSDct_TENET_RealTime_TEmatrix.txt 0.05 || exit 100 &&
echo "RECONSTRUCTION FINISHED" || exit 100 &&

# Trim indirect edges
echo "TRIMMING INDIRECT EDGES OF CSD CONNECTIVE TISSUE GRN" || exit 100 &&
date +"%F %X" || exit 100 &&
echo "python trim_indirect.py [name of GRN] [cutoff]" || exit 100 &&
echo "python trim_indirect.py 20230616_CSDct_TENET_RealTime_TEmatrix.fdr0.05.sif 0" || exit 100 &&
python trim_indirect.py 20230616_CSDct_TENET_RealTime_TEmatrix.fdr0.05.sif 0 || exit 100 &&

echo "DONE. Exiting"
