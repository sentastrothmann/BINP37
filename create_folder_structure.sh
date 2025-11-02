#!/bin/bash
# create_folder_structure.sh

# Root folder
mkdir -p scRNAClassifiers

# 0 - Raw Data
mkdir -p scRNAClassifiers/0_RawData

# 1 - Quality Control
mkdir -p scRNAClassifiers/1_QualityControl/{Data/{Input,Inbetween,Output},Figures,Scripts}

# 2 - Dimensionality Reduction
mkdir -p scRNAClassifiers/2_DimensionalityReduction/{Data/{Input,Inbetween,Output},Figures,Scripts}

# 3 - Data Integration
mkdir -p scRNAClassifiers/3_DataIntegration/{Data/{Input,Inbetween,Output},Figures,Scripts}

# 4 - Clustering
mkdir -p scRNAClassifiers/4_Clustering/{Data/{Input,Inbetween,Output},Figures,Scripts}

# 5 - Differential Gene Expression
mkdir -p scRNAClassifiers/5_DifferentialGeneExpression/{Data/{Input,Inbetween,Output},Figures,Scripts}

# 6 - Cell Type Prediction
mkdir -p scRNAClassifiers/6_CelltypePrediction/{Data/{Input,Inbetween,Output},Figures,Scripts}

# 7 - Trajectory Inference (PAGA)
mkdir -p scRNAClassifiers/7_TrajectoryInferencePAGA/{Data/{Input,Inbetween,Output},Figures,Scripts}

# 8 - Classifiers
mkdir -p scRNAClassifiers/8_Classifiers/{Data/{Input,Inbetween,Output},InitialModels,InitialTraining,Models,Training,Validation,Preprocessing,Figures}
