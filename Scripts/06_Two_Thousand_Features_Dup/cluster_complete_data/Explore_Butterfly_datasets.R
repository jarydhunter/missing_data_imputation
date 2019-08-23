# Explore the butterfly datasets, prior to running ben's clustering scripts on
# them, since I'll need to change ben's scripts to run on the new datasets.

library(tidyverse)
library(e1071)

paths <- c('./data/butterfly/1st_experiment/X1_ori.csv',
           './data/butterfly/1st_experiment/X2_ori.csv',
           './data/butterfly/simulation_script/butterfly_1.csv',
           './data/butterfly/simulation_script/butterfly_2.csv')
names(paths) <- basename(paths)

dfs <- lapply(paths, read_csv, col_names = F)

# There are no labels or headers in the data, they are just csvs of
# floating point values. I'm assuming that the convention is followed that rows
# are observations and columns are different features.

# This function will take a dataframe, and will return a list of image of a
# random ten feature's histograms, a histogram of the means, variances, and skews
# of the features.
get_summary_stats <- function(df){
  stopifnot(require(tidyverse) & require(e1071))
  res <- list()

  gathered <- df %>% gather()

  random_feats <- sample(names(df), 10)
  res[[1]] <- ggplot(gathered %>% filter(key %in% random_feats), aes(x = value)) +
    geom_histogram() + facet_wrap(~key) + ggtitle('random ten features histograms')

  res[[2]] <- gathered %>% group_by(key) %>%
    summarize(abs_mean = mean(value)) %>%
    ggplot(., aes(x = abs_mean)) + geom_histogram()

  res[[3]] <- gathered %>% group_by(key) %>% summarize(var = var(value)) %>%
    ggplot(., aes(x = var)) + geom_histogram()

  res[[4]] <- gathered %>% group_by(key) %>% summarize(skew = skewness(value)) %>%
    ggplot(., aes(x = skew)) + geom_histogram()

  names(res) <- c("ten_hist", "means", "variances", "skews")

  res
}

b1_stats <- get_summary_stats(dfs$butterfly_1.csv)
b2_stats <- get_summary_stats(dfs$butterfly_2.csv)
X1_stats <- get_summary_stats(dfs$X1_ori.csv)
X2_stats <- get_summary_stats(dfs$X2_ori.csv)

b_labels <- read_csv('./data/butterfly/1st_experiment/butterfly_labels.csv', col_names = F)
all_labels <- read_csv('./data/butterfly/1st_experiment/all_labels.csv', col_names = F)


