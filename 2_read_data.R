#### Read and Normalize data ####
# yeast
yeast_prnn <- yeast %>%
  drop_na() %>%
  psrn(load_info = F, id_col = 'identifier')

# UPS
ups_prnn <- ups %>%
  drop_na() %>%
  psrn(load_info = F, id_col = 'identifier')

# If needed, download Ramus or Human DS
get_data()

# Ramus
ramus_prnn <- readr::read_csv('ramus_clean.csv') %>%
  drop_na() %>%
  rename_with(~paste0('condi', rep(1:9, each = 3), '_', rep(1:3, lenght.out = 9*3)), where(is.numeric)) %>%
  psrn(load_info = F, id_col = 'identifier')

# Human
human_prnn <- readxl::read_excel('diaWorkflowResults_allDilutions.xlsx', na = "NA") %>%
  select(-2:-24) %>%
  rename_with(
    ~str_replace_all(.,
      c(
        '\\.\\.\\.1' = 'identifier',
        '1-' = 'spike_prop_'
      )
    ),
    everything()
  ) %>%
  rename_with(
    ~str_remove(., '[0-9]*$') %>%
      paste0(., 1:23),
    where(is.numeric)
  ) %>%
  drop_na()

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
