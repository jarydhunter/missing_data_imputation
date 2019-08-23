#!/hpf/tools/centos6/R/3.1.0/bin/Rscript

# Script to evaluate the how each imputation method affects the
# performance of the clustering methods

argv <- as.numeric(commandArgs(T))

######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
library(impute)
library(survival)
library(tidyverse)

######################################################################
# Initialize folders
projectFolder <- getwd()
testFolder <- paste(projectFolder, "Scripts",
                    "06_Two_Thousand_Features_Dup",
                    "evaluate_imputation", sep="/")
resultsFolder <- paste(projectFolder, "Results", sep="/")

if(!file.exists(resultsFolder)){
  dir.create(resultsFolder)
}

######################################################################
### Leaving these variables as is for now, will need to adjust as I go through
### the rest of this script.
# Initialize fixed variables
jvmGBLimit <- 8
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 12
numFeat <- 2000

# Initialize variable parameters
# Data set which will be tested
dataset <- 'butterfly' #argv[1]
runType <- 1 #argv[2]
seed <- 7 #argv[3]

######################################################################
# Load functions

# Note: some functions depend on variables initialized above!
# lsaImputation:
# -imputedFile, incompleteFile, projectFolder, jvmGBLimit
# iClusterClustering:
# -numCores
# source(paste(projectFolder, "Scripts/loadFunctions.R", sep="/"))
source('./Scripts/loadFunctions.R')


######################################################################
# Load the data

# This function will essentially act as a switch for which dataset to be loaded
# in, then each function can be tailored to the unique requirements of the
# layout of those datasets. Add to the if statements for new datasets. Each
# "load" function should return a list of tibbles.
load_data <- function(dataset){
  if(dataset == 'butterfly'){
    data = load_butterfly()
  } else {
    stop('no loading function for this dataset')
  }
  data
}

tibs <- load_data(dataset = dataset)

dat <- tibs[['dat']]

######################################################################
# Subset the data into complete and incomplete observations, complete set will
# bet used to create a false missingness to evaluate the quality of the four
# imputation methods Ben used on the TCGA data, the best method will be passed to
# the next part of the analysis to impute the true missing values which will have
# SNF applied to it to generate the clusters.

# Get all observations which are complete
completeData <- dat %>% filter_all(all_vars(!is.na(.)))

# Generate a false missingness in the completeData
set.seed(seed)

addedMissingness <- generateIncompleteIntersection(dat)

# Extract all cases which appear in all of the data types in the
# incomplete data

# ######################################################################
# # Select a subset of features which differ most between cases and
# # controls.
#
# featureSubsetIndices <- function(cases, subsetSize=numFeat) {
#   numViews <- length(cases)
#   featureSubsetInd <- vector("list", numViews)
#
#   for (v in 1:numViews) {
#     # Calculate the t-test p-value for each feature, grouped by cases
#     # and controls
#     numFeatures <- nrow(cases[[v]])
#     pval <- sapply(1:numFeatures,
#                    function(i) t.test(cases[[v]][i, ],
#                                       controls[[v]][i, ])$p.value)
#
#     # Subset the data keeping the features with the smallest p-values
#     ind <- order(pval)
#     featureSubsetInd[[v]] <- ind[1:min(subsetSize, numFeatures)]
#   }
#
#   return(featureSubsetInd)
# }
#
# subsetData <- function(data, ind) {
#   for (v in 1:length(data)) {
#     data[[v]] <- data[[v]][ind[[v]], ]
#   }
#
#   return(data)
# }
#
# completeInd <- featureSubsetIndices(completeData)
# completeData <- subsetData(completeData, completeInd)
#
# incompleteInd <- featureSubsetIndices(incompleteData)
# incompleteData <- subsetData(incompleteData, incompleteInd)
# removedData <- subsetData(removedData, incompleteInd)

######################################################################
# Normalize the features in the data sets.
# Normalization is performed before imputation and we expect that the
# data will still be normalized after imputation (before clustering).

# rowStatistics <- function(cases) {
#   numViews <- length(cases)
#   rowStats <- vector("list", numViews)
#
#   for (v in 1:numViews) {
#     # Calculate the row means and standard deviations
#     rowMean <- apply(cases[[v]], 1, mean, na.rm=TRUE)
#     rowSd <- apply(cases[[v]], 1, sd, na.rm=TRUE)
#     constantInd <- rowSd==0
#     rowSd[constantInd] <- 1
#     rowStats[[v]] <- list(mean=rowMean, sd=rowSd, ind=constantInd)
#   }
#
#   return(rowStats)
# }
#
# normalizeData <- function(data, stat) {
#   for (v in 1:length(data)) {
#     data[[v]] <- (data[[v]] - stat[[v]]$mean) / stat[[v]]$sd
#     data[[v]] <- data[[v]][!stat[[v]]$ind, ]
#   }
#
#   return(data)
# }
#
# completeStat <- rowStatistics(completeData)
# completeData <- normalizeData(completeData, completeStat)
#
# incompleteStat <- rowStatistics(incompleteData)
# incompleteData <- normalizeData(incompleteData, incompleteStat)
# removedData <- normalizeData(removedData, incompleteStat)


######################################################################
# Evaluate the performance of the methods

# Flag indicating if the rows are the observations or not.
sampleRows <- FALSE
clusteringMethods <- c(hierarchicalClustering, iClusterClustering,
                       SNFClustering)
imputationMethods <- c('knnImputation', 'llsImputation', 'lsaImputation',
                       'randomImputation')

### The knn method takes weird input, it needs genes as rows and "samples" as columns,
### this is probably why Ben was making all features rows, I'd rather just
### transpose it and then untranspose it after.

# Impute the missing data and save the results
impute <- function(x) get(imputationMethods[[runType]])(x)

imputationOutput <- evaluateImputation(impute, as.matrix(addedMissingness),
                                       as.matrix(completeData))
imputedData <- imputationOutput$imputedData
imputedDataFile <- paste0(resultsFolder, '/', imputationMethods[[runType]], '_imputed_data.csv')
write_csv(as_tibble(imputedData), imputedDataFile)

imputationResults <- imputationOutput$results
imputationResultsFile <- paste0(resultsFolder, '/', imputationMethods[[runType]], "_summary_stats.csv")
write_csv(as_tibble(imputationResults), imputationResultsFile)

### In this new version I won't be comparing the different clustering methods for
### each imputation method, instead I will try each imputation method and
### choose the one that is the closest as the best imputation method, then apply
### the snf clustering method to that best imputation.
# Save the results of clustering the imputed and intersected data
# for (i in 1:length(clusteringMethods)) {
#   completeLabels <- as.matrix(tibs[['labels']])
#
#   # Initialize clustering function
#   numClus <- length(unique(completeLabels))
#   cluster <- function(x) clusteringMethods[[i]](x, numClus)
#
#   clusteringResults <- evaluateClustering(cluster, imputedData,
#                                           completeLabels)
#   writeResults(c(argv, i, j, clusteringResults), clusteringFile)
# }
