#' Calculate imputation errors
#'
#' For a given tibble of completeData, and imputedData a tibble with imputed
#' values, caculate the mse, rmse, nrmse, rae for the true values with the
#' imputed version.
#'
#' @param original: a tibble or the original values
#' @param imputed: a tibble where some values are imputed
#' @param num_missing: integer of how many values were imputed
#'
#' @return named list of values, mse, rmse, nrmse, and rae
imputationErrors <- function(original, imputed, num_missing) {

  # using num_missing, since we know all non-imputed entries in imputed will be
  # exactly the same as the original so those will all sum to 0.
  mse <- sum((original-imputed)^2)/num_missing
  rmse <- mse^0.5 # RMSE
  # using as definition of nrmse, rmse/sd(data), since for some datasets we
  # have a mean very very close to zero, which can cause dividing by it to
  # explode which can make this measure uninterpretable.
  nrmse <- rmse/sd(as.matrix(original))

  absError <- abs(original - imputed)
  relAbsError <- absError / abs(original)
  rae <- sum(relAbsError)/num_missing # RAE

  list(mse=mse, rmse=rmse, nrmse=nrmse, rae=rae)
}