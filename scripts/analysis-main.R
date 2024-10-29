

########## TASK 1


# examine distribution of non transformed data

setwd('/Users/davidpan/Desktop/PSTAT197/module1-group-1')
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()



# read in data
biomarker_no_transform <- read_csv('data/biomarker-raw.csv', 
                                   skip = 2,
                                   col_select = -2L,
                                   col_names = c('group', 
                                                 'empty',
                                                 pull(var_names, abbreviation),
                                                 'ados'),
                                   na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # reorder columns
  select(group, ados, everything())



# pick a random protein to see its distribution
ggplot(biomarker_no_transform, aes(x=SMAD2, color=group, fill=group)) +
  geom_histogram(alpha=0.5)
# compare with transformed data
ggplot(biomarker_clean, aes(x=SMAD2, color=group, fill=group)) +
  geom_histogram(alpha=0.5)



# pick another protein
ggplot(biomarker_no_transform, aes(x=BTC, color=group, fill=group)) +
  gem_histogram(alpha=0.5)
# compare with transformed data
ggplot(biomarker_clean, aes(x=BTC, color=group, fill=group)) +
  geom_histogram(alpha=0.5)





########## TASK 2


# outliers will have values of 3 or -3
# count number of outliers for each subject
outlier_count <- rowSums(biomarker_clean[, 3:1319] == 3 | biomarker_clean[, 3:1319] == -3)

outlier_by_id <- cbind(c(1:154), outlier_count) %>%
  as_tibble() %>% rename(id=V1) %>% arrange(desc(outlier_count))

# average number of outliers in all subjects
mean(outlier_count)

# top 10 subjects with most outliers
# these subjects have an outlier count much greater than the mean
# subjects with more than 50 outliers: 154, 108, 9, 121, 52, 77, 147
head(outlier_by_id, 10)

# examine whether group affects outliers
biomarker_no_transform <- biomarker_no_transform %>% mutate(id=row_number()) %>% relocate(id)
biomarker_outliers <- merge(biomarker_no_transform, outlier_by_id, by='id')[,c('id','group','ados','outlier_count')]

ASD_mean_out <- mean(biomarker_outliers[biomarker_no_transform$group=='ASD',]$outlier_count)
TD_mean_out <- mean(biomarker_outliers[biomarker_no_transform$group=='TD',]$outlier_count)

ASD_median_out <- median(biomarker_outliers[biomarker_no_transform$group=='ASD',]$outlier_count)
TD_median_out <- median(biomarker_outliers[biomarker_no_transform$group=='TD',]$outlier_count)

# the TD group has a higher mean outlier count
# the difference does not seem to be substantial

# look at the distribution of outliers between groups
# they have similar distributions, with the TD group having slightly more outliers that skew the mean
ggplot(biomarker_outliers, aes(x=outlier_count, color=group, fill=group)) + 
  geom_histogram(binwidth=10, alpha=0.5) + theme_minimal() + facet_grid(~group) + 
  labs(title='Distribition of Outlier Count by Group', x='Outliers', y='Frequency')


# in the ASD group, does ados score have any correlation with outlier count?
# visually, there is no relationship
ggplot(biomarker_outliers[biomarker_outliers$group=='ASD',], aes(x=ados, y=outlier_count)) + 
  geom_point()


library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(dplyr)
load("data/biomarker-clean.RData")

########## TASK 3

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

# use a fuzzy intersection by considering the adj. p-value from t-tests and Mean Decrease Gini
# from the RF together
top_proteins_s1 <- proteins_s1 %>%
  slice_min(p.adj, n = 3)

top_proteins_s2 <- proteins_s2 %>%
  slice_max(MeanDecreaseGini, n = 3)

# Get the intersection of proteins in proteins_s1 and proteins_s2
intersection_proteins <- intersect(proteins_s1$protein, proteins_s2$protein)

# Combine top proteins and intersection
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

# before using fuzzy intersection(in class analysis)
# 1 sensitivity binary         0.75 
# 2 specificity binary         0.8  
# 3 accuracy    binary         0.774
# 4 roc_auc     binary         0.871

# after using fuzzy intersection(worse)
# 1 sensitivity binary         0.562
# 2 specificity binary         0.867
# 3 accuracy    binary         0.710
# 4 roc_auc     binary         0.746


########## TASK 4

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

# Combine top proteins from both methods
proteins_top <- bind_rows(top_proteins_s1, top_proteins_s2) %>%
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

## task 4 results

# choose 20 instead of 10 proteins from each method and 
# intersect the intersection of these two sets with the 
# union of the top 10 proteins to get the final panel of 
# proteins (5)

# 1 sensitivity binary         0.812
# 2 specificity binary         0.867
# 3 accuracy    binary         0.839
# 4 roc_auc     binary         0.871
