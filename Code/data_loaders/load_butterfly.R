library(tidyverse)

#'Load butterfly data
#'
#' Load the csv files for the butterfly dataset, will return a list tibbles. The
#' first tibble is all observations where some features are missing for some
#' observations, and the second tibble is the labels for each observation.
#'
#' @param path the path for the folder containing all csv files belonging to the
#' butterfly dataset
#'
#' @return a list of tibbles
load_butterfly <- function(path = '/data/Jaryd/missing_data_imputation/data/butterfly/1st_experiment/'){
  files <- paste0(path, list.files(path) %>% grep('\\.csv', x = ., value = T))
  dfs <- lapply(files, read_csv, col_types = cols(.default = col_double()), col_names = F)
  names(dfs) <- basename(files)

  # Adding an ID column for my own ease of joining
  dfs[['X1_ori.csv']] <- dfs[['X1_ori.csv']] %>% mutate(ID = 1:length(.[[1]]))
  dfs[['X2_ori.csv']] <- dfs[['X2_ori.csv']] %>%
    rename_all(.funs = ~ str_replace(., 'X', 'Y')) %>%
    mutate(ID = c(1:832, (length(dfs[['X1_ori.csv']][[1]])+1):(length(dfs[['X1_ori.csv']][[1]])+300)))
  dfs[['all_labels.csv']] <- dfs[['all_labels.csv']] %>%
    rename() %>%
    mutate(ID = 1:length(.[[1]]))

  list(dat =
         full_join(dfs[['X1_ori.csv']], dfs[['X2_ori.csv']], by = 'ID') %>%
         select(-ID),
       labels = dfs[['all_labels.csv']] %>% select(-ID))
}



