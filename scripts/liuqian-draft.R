library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(dplyr)
load("data/biomarker-clean.RData")

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

ttests_out <- biomarker_clean %>%
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
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()



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

# use a fuzzy intersection by considering the adj. p-value from t-tests and Mean Decrease Gini 
# from the RF together 
top_proteins_s1 <- proteins_s1 %>% 
  slice_min(p.adj, n = 3)
  
top_proteins_s2 <- proteins_s2 %>% 
  slice_max(MeanDecreaseGini, n = 3)

# Get the intersection of proteins in proteins_s1 and proteins_s2
intersection_proteins <- intersect(proteins_s1$protein, proteins_s2$protein)

# Combine top proteins and intersection, removing duplicates
proteins_sstar <- bind_rows(top_proteins_s1, top_proteins_s2) %>%
  filter(protein %in% intersection_proteins | protein %in% c(top_proteins_s1$protein, top_proteins_s2$protein)) %>%
  distinct(protein, .keep_all = TRUE) %>% 
  pull(protein)

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

# before using fuzzy intersection 
# 1 sensitivity binary         0.812
# 2 specificity binary         0.733
# 3 accuracy    binary         0.774
# 4 roc_auc     binary         0.883

# after using fuzzy intersection(better)
# 1 sensitivity binary         0.875
# 2 specificity binary         0.867
# 3 accuracy    binary         0.871
# 4 roc_auc     binary         0.908