#' Evaluate Imputation
#'
#' Given a function for imputation impute, incompleteData a tibble with
#' missingness artificially applied and it's corresponding complete tibble
#' completeData, calculate and return metrics for how well the imputation
#' performed.
#'
#' @param impute: function used to impute the incomplete data
#' @param incompleteData: numeric matrix with missing values artificially applied
#' @param completeData: numeric matrix with no missingness
#'
#' @return list of metrics for how well impute performed on incompleteData
#'
#'
evaluateImputation <- function(impute, incompleteData, completeData) {

  num_missing <- is.na(incompleteData) %>% sum()
  # Impute the incomplete data
  runTime <- system.time(imputedData <- impute(incompleteData))[3]

  # Evaluate the imputed data
  # Note: the metrics for each data type is calculated separately
  nrmse <- imputationErrors(completeData, imputedData, num_missing)

  results <- c(nrmse, runTime)

  list(imputedData=imputedData, results=results)
}