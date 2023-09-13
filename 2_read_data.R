#### Setup design matrix, Read-, Normalize- data, and calc M-V trend ####
# yeast
yeast_design <- model.matrix(~0+factor(rep(1:2, each = 3)))
colnames(yeast_design) <- paste0('ng', c(50, 100))

yeast_prnn <- load_data('yeast') %>%
  calculate_mean_sd_trends(yeast_design)

# UPS
ups_design <- model.matrix(~0+factor(rep(1:3, each = 4)))
colnames(ups_design) <- paste0('fmol', c(25, 50, 100))

ups_prnn <- load_data('ups') %>%
  calculate_mean_sd_trends(ups_design)

# If needed, download Ramus, Bruderer, or Human DS
get_data()

# Ramus
ramus_design <- model.matrix(~0+factor(rep(1:9, each = 3)))
colnames(ramus_design) <- paste0('condi', 1:9)

ramus_prnn <- load_data('ramus') %>%
  calculate_mean_sd_trends(ramus_design)

# Human
human_design <- model.matrix(~0+factor(rep(1:3, each = 23)))
colnames(human_design) <- paste0('spike_prop_', c(25, 12, 6))

human_prnn <- load_data('human') %>%
  calculate_mean_sd_trends(human_design)

# Bruderer
bruder_design <- ramus_design[1:24,1:8]

bruder_prnn <- load_data('bruderer') %>%
  calculate_mean_sd_trends(bruder_design)

# Navarro
navarro_design <- yeast_design
colnames(navarro_design) <- paste0('condi', 1:2)

navarro_prnn <- load_data('navarro') %>%
  calculate_mean_sd_trends(navarro_design)

#### Define auxiliary variables ####
full_page <- 170
half_page <- 85
dpi       <- 300
workers   <- floor(parallel::detectCores()/2)

ramus_str_replace <- c(
  'condi1' = '0.05',
  'condi2' = '0.125',
  'condi3' = '0.25',
  'condi4' = '0.5',
  'condi5' = '2.5',
  'condi6' = '5',
  'condi7' = '12.5',
  'condi8' = '25',
  'condi9' = '50'
)

ramus_plot_order <- ramus_str_replace %>%
  combn(2) %>%
  t() %>%
  apply(1, str_flatten, ' vs ')

human_str_replace <- c(
  'spike_prop_6' = '1:6',
  'spike_prop_12' = '1:12',
  'spike_prop_25' = '1:25'
)

human_plot_order <- human_str_replace %>%
  rev() %>%
  combn(2) %>%
  t() %>%
  apply(1, str_flatten, ' vs ')

facet_order <- c(
  paste0(c('False', 'True'), ' Positive Rate'), 'Precision'
)
