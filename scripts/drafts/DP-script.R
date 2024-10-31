library(tidyverse)


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


