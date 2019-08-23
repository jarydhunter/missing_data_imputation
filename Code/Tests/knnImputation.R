# Keeping track of the splitting conditions would be easier if there
# was a single variable which could be used to keep track of the logic

#' KNN imputation
#'
#' A wrapper around the knn.impute function, this function takes a matrix in
#' standard form, i.e. each row is an observation and each column is a feature,
#' then it transposes it so that knn.impute can be called on it, and then it
#' tranposes the results and returns it. The matrix being passed in is assumed to
#' be numeric, since knn.impute relies on euclidean distance thus all
#' assumptions for that also hold. Missingness is assumed to be less than 50%
#' for a given feature, and less than 80% for a given observation. For more
#' details about the imputation see `?impute.knn`
#'
#' @param incompleteData a matrix of entirely numeric data with some missingness.
#'
#' @return a tibble with NA's replaced with imputed values.
knnImputation <- function(incompleteData) {
  stopifnot(require(impute))
  # Transpose the data for KNN
  # Tranpose the matrix so that each column is an observation and each row is a
  # feature, this will also create an extra column feature to preserve the
  # feature names, and this will be used to undo the transposition at the end of
  # the function.
  # Note: this "tidyverse" method of transposing does not preserve the order of
  # the columns, so the column order is stored before transposing and then
  # reapplied at the end.

  incompleteData <- t(incompleteData)

  # Impute the missing data, requirements for impute.knn is that the data is a
  # matrix, and all columns must be numeric.
  imputedData <- incompleteData %>%
    impute.knn(., rowmax=1, colmax=1) %>%
    .$data

  # Untranspose the data
  t(imputedData)
}