#'local least squares imputation
#'
#' Performs imputation on missing values using the local least squares method.
#'
#' @param incompleteData a matrix of the data with some missingness, and each
#' row is an observation
#'
#' @return a matrix where all missingness has been replaced with imputed values
#' using the local least squares imputation algorithm.
#'
llsImputation <- function(incompleteData) {
  # Impute the missing data
  llsImpute(incompleteData, allVariables=TRUE)
}