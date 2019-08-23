# Return metrics to evaluate clustering the data using the clustering
# of the complete data and survival data.
# cluster: function used to cluster the data
# data: list of complete matrices
# completeLabels: numeric vector of labels for the data attained from
# clustering the complete data with the cluster function
# clinicalData: data frame of TCGA clinical data

evaluateClustering <- function(cluster, data, completeLabels,
                               clinicalData) {
  browser()
  # Cluster the data
  runTime <- system.time(labels <- cluster(data))[3]

  # Evaluate the new labels
  acc <- calAcc(completeLabels, labels)
  nmi <- calNMI(completeLabels, labels)
  surv <- evaluateSurvival(clinicalData, labels)
#
#   # Evaluate the new labels on individuals in the intersection
#   intLabels <- labels[intersectInd]
#   intCompleteLabels <- completeLabels[intersectInd]
#   intClinicalData <- clinicalData[intersectInd, ]
#   intAcc <- calAcc(intCompleteLabels, intLabels)
#   intNmi <- calNMI(intCompleteLabels, intLabels)
#   if (length(unique(intLabels))==1) {
#     intSurv <- c(NA, NA)
#   } else {
#     intSurv <- evaluateSurvival(intClinicalData, intLabels)
#   }
  results <- c(acc, nmi, surv) #intAcc, intNmi, intSurv, runTime)

  return(results)
}
##### Comments below was just in the script, not sure why, commenting out
##### because it causes errors, when sourcing this file, but it could be the
##### code didn't upload properly and there is more missing so it is now incorrect.
#pathway analysis
#real data set run clustering, survival analysis,
#send link