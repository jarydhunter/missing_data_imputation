# missing_data_imputation
This repo houses the pipeline for the missing data imputation project.

# Data types and Terms
- Three biological data types : (1) methylation, (2) mrna, and (3) mirna
- 5 cancers: BRCA, LIHC, KIRC, LUSC, LUAD
- complete data (intersection) - the sample intersection across the three data types for each cancer. There are missing samples 
in each data type, so we only take the samples present in all three, to get a complete data set.
- original data (union) - the sample union of all data types.

# 06_Two_Thousand_Features_Dup

This folder contains 5 subfolders (1) cluster_complete_data, (2) evaluate_imputation, 
(3) evaluate_similarity, (4) evaluate_original_imputation, (5) evaluate_origingal_similarity. 

## (1) cluster_complete_data

This folder is the first step in the pipeline - two scripts (1) job.R and (2) run.sh. 

(1) Job.R contains the code for reading in the datasets into a list, normalizing the data, and 
subsetting the amount of features. Then we take the intersection of all samples across 
the three data types, so we have a "complete" data set. Then it clusters the data using all clustering methods and saves the results in the Results folder. Any errors will be sent to the Error folder. 

(2) run.sh is a shell script to execute job.R on the cluster, executing all cancer data sets and clustering methods
in parrellel. 

## evaluate_imputation 

The second step is to run Job.R within this folder - it will read in the results from cluster_complete_data as 
ground truth. In Job.R, we run 50 simulations, each one removing a random part of the "complete data" (using the same 
structure of missingness present in the original data). This script employs 4 imputation methods: (knn, random, LSA, LSS) to impute directly on the biological data. After imputation and clusters the newly imputed data and is evaluated against the ground truth using NMI, ACC, P value (survival) and concordance index (survival). 

## evaluate_similarity

This is still part of the second step and is the exact same as the evaluate_imputation method, but instead of imputing directly on the biological data, we impute on the similarity matrices computed during the SNF algorithm, prior to spectral clustering. This is a much faster approach and can be implemented off cluster.

## evaluate_original_imputation

After the second step, we have our results from the "complete data" (intersection) and now, in the 3rd step, we apply the optimal imputation/clustering combination on the original data (union) and evaluate using survival pvalue and CI (can't use ACC and NMI because there is no ground truth on the original data).

## evaluate_original_similarity

Similar to the evaluate_original_imputation folder, this is part of the third step and does the same thing, but imputes on the affinity matrices. 

# analyze_results

This folder contains two scripts (1) run.R and (2) run_original.R. Both are used to evaluate the results from steps 2 and 3 and plot the outcomes. The important thing here is to see if the results from the intersection generalize to the union. 
