library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(dplyr)
library(fuzzyjoin)
library(stringdist)
load("data/biomarker-clean.RData")

## task 3 
library(tidyverse)


# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

var_names

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}


# read in data
biomarker_clean <- read_csv('data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  # reorder columns
  select(group, ados, everything())

# export as r binary
save(list = 'biomarker_clean', 
     file = 'data/biomarker-clean.RData')

# partition into training and test set before the variable selection process
set.seed(101422)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.8)

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- training(biomarker_split) %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  select(protein, p.adj)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- training(biomarker_split) %>%
  select(-c(group, ados))

response <- training(biomarker_split) %>% pull(group) %>% factor()



# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  select(protein, MeanDecreaseGini)

## LOGISTIC REGRESSION
#######################

# # select subset of interest
# proteins_sstar <- intersect(proteins_s1$protein, proteins_s2$protein)
# 


# Fuzzy intersection with a specified similarity threshold
fuzzy_intersection <- stringdist_inner_join(proteins_s1, proteins_s2,
                                            by = "protein",
                                            max_dist = 0.46,  # maximum allowed distance
                                            method = "jw") # Jaro-Winkler similarity

# Combine the matched protein names from both columns and keep unique values
proteins_sstar <- unique(c(fuzzy_intersection$protein.x, fuzzy_intersection$protein.y))

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_training <- training(biomarker_split) %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_testing <- testing(biomarker_split) %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = biomarker_sstar_training, 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

biomarker_sstar_testing %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

# before using fuzzy intersection(in class analysis"DERM" "IgD" "TSP4" "FSTL1")
# 1 sensitivity binary         0.75 
# 2 specificity binary         0.8  
# 3 accuracy    binary         0.774
# 4 roc_auc     binary         0.871

# after using fuzzy intersection(no significant improvement"DERM" "IgD" "TSP4" "FSTL1" "MAPK14")
# 1 sensitivity binary         0.75 
# 2 specificity binary         0.8  
# 3 accuracy    binary         0.774
# 4 roc_auc     binary         0.888


## task 4

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- training(biomarker_split) %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  select(protein, p.adj)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- training(biomarker_split) %>%
  select(-c(group, ados))

response <- training(biomarker_split) %>% pull(group) %>% factor()



# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  select(protein, MeanDecreaseGini)

## LOGISTIC REGRESSION
#######################

# # select subset of interest
# proteins_sstar <- intersect(proteins_s1$protein, proteins_s2$protein)

# use a fuzzy intersection by considering the adj. p-value from t-tests and Mean Decrease Gini
# from the RF together
top_proteins_s1 <- proteins_s1 %>%
  slice_min(p.adj, n = 10)

top_proteins_s2 <- proteins_s2 %>%
  slice_max(MeanDecreaseGini, n = 10)

# Get the intersection of proteins in proteins_s1 and proteins_s2
intersection_proteins <- intersect(proteins_s1$protein, proteins_s2$protein)

# Combine top proteins and intersection
proteins_top <- bind_rows(top_proteins_s1, top_proteins_s2) %>%
  filter(protein %in% c(top_proteins_s1$protein, top_proteins_s2$protein)) %>%
  distinct(protein, .keep_all = TRUE) %>%
  pull(protein)

proteins_sstar <- intersect(intersection_proteins, proteins_top)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_training <- training(biomarker_split) %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_testing <- testing(biomarker_split) %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = biomarker_sstar_training, 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

biomarker_sstar_testing %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

## task 4 

# choose 20 instead of 10 proteins from each method and 
# intersect the intersection of these two sets with the 
# union of the top 10 proteins to get the final panel of 
# proteins (5)

# 1 sensitivity binary         0.812
# 2 specificity binary         0.867
# 3 accuracy    binary         0.839
# 4 roc_auc     binary         0.871

