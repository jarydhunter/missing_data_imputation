# Script to evaluate the how each imputation method affects the
# performance of the clustering methods

argv <- commandArgs(T)

if (length(argv)<= 1) {
    stop("A home directory and dataset must be supplied")
} else if(length(argv)==2){
    argv[3] <- 1 # default imputation method
    argv[4] <- 42 # default seed
} else if(length(argv) == 3){
    argv[4] <- 42
}
######################################################################
# Load libraries
library(SNFtool)
library(iClusterPlus)
library(impute)
library(survival)
library(tidyverse)

######################################################################
# Initialize folders
projectFolder <- argv[1]
setwd(projectFolder)
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
cancerTypes <- c("BRCA", "KIRC", "LIHC", "LUAD", "LUSC")
dataTypes <- c("methyl", "mirna", "mrna")
numCores <- 12
numFeat <- 2000

# Initialize variable parameters
# Data set which will be tested
dataset <- argv[2]
runType <- as.numeric(argv[3])
seed <- as.numeric(argv[4])

######################################################################
# Load functions

# Note: some functions depend on variables initialized above!
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
# be used to create a false missingness to evaluate the quality of the four
# imputation methods Ben used on the TCGA data, the best method will be passed to
# the next part of the analysis to impute the true missing values which will have
# SNF applied to it to generate the clusters.

# Get all observations which are complete
completeData <- dat %>% filter_all(all_vars(!is.na(.)))

# Generate a false missingness in the completeData
set.seed(seed)

addedMissingness <- generateIncompleteIntersection(dat)

######################################################################
# Evaluate the performance of the imputation methods, and save the resulting
# imputed matrices.

# clusteringMethods <- c(hierarchicalClustering, iClusterClustering,
#                        SNFClustering)
imputationMethods <- c('knnImputation', 'llsImputation', 'lsaImputation')

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

### In this new version I won't bundle the different clustering methods with
### each imputation method, instead I will try each imputation method and
### store the resulting object, then apply each clustering method.
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
