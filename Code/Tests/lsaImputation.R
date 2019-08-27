#' least squares adaptive imputation
#'
#' Perform least squares adaptive imputation on the given matrix which already
#' has some missingness. This function calls a function written in Java to perform
#' the imputation.
#'
#' @param incompleteData a matrix with missingness to be imputed upon
#' @param jvmGBLimit the memory limit in Gigabytes for the java virtual machine
#'
#' @return a matrix with no missingness, all missingness imputed by the LSA
#' algorithm.
#'
lsaImputation <- function(incompleteData, jvmGBLimit = 8) {
  browser()
  projectFolder <- getwd()
  incompleteFile <- paste0(projectFolder, '/incompleteFile')
  imputedFile <- paste0(projectFolder, '/Results/imputedFile')

  # Concatenate the data by columns because the samples are in rows
  # Row and column names are required for the LSA code
  column_names <- colnames(incompleteData)
  row_names <- rownames(incompleteData)
  colnames(incompleteData) <- as.character(1:ncol(incompleteData))
  rownames(incompleteData) <- as.character(1:nrow(incompleteData))

  # Save the incomplete data to a file
  cat("xxx","\t", file=incompleteFile, sep="")
  write.table(incompleteData, incompleteFile,
              col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE,
              na="NULL",append=TRUE)

  # Impute the missing data
  lsaCodePath <- paste(projectFolder, "Code/Impute/LSimpute.jar",
                       sep="/")
  jvmMemLimit <- paste("-Xmx", jvmGBLimit, "G", sep="")
  javaCommand <- paste("java -jar", jvmMemLimit, "-server ",
                       lsaCodePath, incompleteFile, imputedFile, "6",
                       sep=" ")
  system(javaCommand)

  # Read the imputed file
  imputedData <- data.matrix(read.table(imputedFile, header=TRUE, sep="\t"))
  imputedData <- imputedData[, -1]
  colnames(imputedData) <- column_names
  rownames(imputedData) <- row_names

  # Remove the files used for imputation
  system(paste("rm", imputedFile, sep=" "))
  system(paste("rm", incompleteFile, sep=" "))

  imputedData
}