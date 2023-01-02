#### Read and Normalize data ####
# yeast
yeast_design <- model.matrix(~0+factor(rep(1:2, each = 3)))
colnames(yeast_design) <- paste0('ng', c(50, 100))

yeast_prnn <- load_data('yeast') %>%
   # calculate_mean_sd_trends(yeast_design) %>%
   # mutate(
   #   across(matches(get_conditions(yeast_design)), ~ (. - mean(mean))/sd(mean))
   # ) %>%
  calculate_mean_sd_trends(yeast_design)

# UPS
ups_design <- model.matrix(~0+factor(rep(1:3, each = 4)))
colnames(ups_design) <- paste0('fmol', c(25, 50, 100))

ups_prnn <- load_data('ups') %>%
   # calculate_mean_sd_trends(ups_design) %>%
   # mutate(
   #   across(matches(get_conditions(ups_design)), ~ (. - mean(mean))/sd(mean))
   # ) %>%
  calculate_mean_sd_trends(ups_design)

# If needed, download Ramus or Human DS
get_data()

# Ramus
ramus_design <- model.matrix(~0+factor(rep(1:9, each = 3)))
colnames(ramus_design) <- paste0('condi', 1:9)

ramus_prnn <- load_data('ramus') %>%
   # calculate_mean_sd_trends(ramus_design) %>%
   # mutate(
   #   across(matches(get_conditions(ramus_design)), ~ (. - mean(mean))/sd(mean))
   # ) %>%
  calculate_mean_sd_trends(ramus_design)

# Human
human_design <- model.matrix(~0+factor(rep(1:3, each = 23)))
colnames(human_design) <- paste0('spike_prop_', c(25, 12, 6))

human_prnn <- load_data('human') %>%
  #trend_imputation(human_design) %>%
   # calculate_mean_sd_trends(human_design) %>%
   # mutate(
   #   across(matches(get_conditions(human_design)), ~ (. - mean(mean))/sd(mean))
   # ) %>%
  calculate_mean_sd_trends(human_design)

# Auxiliary variables
full_page <- 170
half_page <- 85
uni       <- 'mm'
dpi       <- 300
workers   <- floor(parallel::detectCores()/2)

ramus_str_replace <- c(
  'condi1' = 'fmol0.05',
  'condi2' = 'fmol0.125',
  'condi3' = 'fmol0.25',
  'condi4' = 'fmol0.5',
  'condi5' = 'fmol2.5',
  'condi6' = 'fmol5',
  'condi7' = 'fmol12.5',
  'condi8' = 'fmol25',
  'condi9' = 'fmol50'
)

human_str_replace <- c(
  'spike_prop_6' = '1:6',
  'spike_prop_12' = '1:12',
  'spike_prop_25' = '1:25'
)
