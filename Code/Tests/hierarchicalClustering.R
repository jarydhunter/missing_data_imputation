#' Hierarchical Clustering
#'
#' Given a numeric matrix, where each row is an observation and each column is a
#' feature, it will run hclust using correlation as a unit of distance, and
#' average distance as the linkage function.
#'
#' @param dat: matrix of data to be clustered
#' @param numClus: the number of clusters to be found
#'
#' @return labels assigned by hierarchical clustering on the matrix data
hierarchicalClustering <- function(dat, numClus) {
  # Transpose the data for the cor function
  # The cor function takes correlations between columns
  dat <- t(dat)

  # Calculate the correlation between samples
  # The cor function takes correlations between columns
  correlation <- cor(dat)

  # Convert the correlation to a distance object
  distance <- as.dist(1 - correlation)

  # Perform clustering
  hclustFit <- hclust(distance, method="average")
  labels <- cutree(hclustFit, k=numClus)

  # Return the resulting labels
  labels
}