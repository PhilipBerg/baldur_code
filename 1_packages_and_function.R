if (!("pacman" %in% .packages(all.available = T))) {
  install.packages("pacman")
  library("pacman")
} else if(!("pacman" %in% (.packages()))) {
  library("pacman")
}
if (!("devtools" %in% .packages(all.available = T))) {
  install.packages("devtools")
}

if (!("BiocManager" %in% .packages(all.available = T))) {
  install.packages("BiocManager")
  BiocManager::install(version = "3.12")
}
if (!("baldur" %in% .packages(all.available = T))) {
  devtools::install_github("PhilipBerg/baldur")
}

p_load(fitdistrplus, vroom, magrittr, tibble, plyr, dplyr,
       purrr, tidyr, stringr, readxl, multidplyr,
       ggplot2, future, furrr, janitor, baldur, rstan)
if (!("limma" %in% .packages(all.available = T))) {
  BiocManager::install("limma")
  library("limma")
}else{
  library("limma")
}

calc_mcfadden <- function(model){
  1 - model$deviance / model$null.deviance
}


#### Data related functions ####
get_human_data <- function() {
  download.file('https://zenodo.org/record/6379087/files/Source%20Data.zip?download=1', 'human_data.zip', method = 'curl')
  unzip('human_data.zip')
  file.rename('Source Data/MinimumDataset/diaWorkflowResults_allDilutions.xlsx', 'diaWorkflowResults_allDilutions.xlsx')
  file.remove('human_data.zip')
  unlink('__MACOSX', TRUE)
  unlink('Source Data', TRUE)
}
get_ramus_data <- function() {
  readr::read_csv('https://figshare.com/ndownloader/files/35592290?private_link=28e837bfe865e8f13479', show_col_types = FALSE) %>% 
    janitor::clean_names() %>%
    mutate(
      across(where(is.numeric), ~na_if(.x, 0))
    ) %>% 
    readr::write_csv('ramus_clean.csv')
}
get_data <- function() {
  if(!file.exists('diaWorkflowResults_allDilutions.xlsx')) {
    get_human_data()
  }
  if(!file.exists('ramus_clean.csv')) {
    get_ramus_data()
  }
}

#### Performance ####
calc_tpr_fpr <- function(alpha, hits, p, n, p_val_col){
  p_val <- pull(hits, p_val_col)
  tp <- hits$tp
  tn <- hits$tn
  tibble(
    TP = sum((p_val<=alpha) & tp),
    TN = sum((p_val>=alpha) & tn),
    FP = sum((p_val<=alpha) & tn),
    FN = sum((p_val>=alpha) & tp)
  ) %>%
    mutate(
      TPR = TP/p,
      FPR = FP/n,
      NPV = TN/(TN+FN),
      precision = TP/(TP + FP),
      MCC = sqrt(precision)*sqrt(TPR)*sqrt(1-FPR)*sqrt(NPV) - sqrt(1-precision)*sqrt(1-TPR)*sqrt(FPR)*sqrt(1-NPV)
    )
}

create_roc <- function(p_val_col, data, regex, cl = NULL){
  data <- data %>%
    split.data.frame(.$comparison)
  out <- map(data,
             pull, p_val_col
  ) %>%
    map(unique) %>%
    map(c, 1.1) %>%
    map(~c(-0.1, .)) %>%
    map(enframe, value = 'alpha', name = NULL)
  data <- data %>%
    map(
      mutate,
      tp = str_detect(identifier, regex),
      tn = !tp
    )
  p <- data %>%
    map(
      use_series, tp
    ) %>%
    map(sum)
  n <- data %>%
    map(
      use_series, tn
    ) %>%
    map(sum)
  pmap(list(out, data, p, n), create_roc_helper, cl, p_val_col) %>%
    bind_rows()
}

create_roc_helper <- function(out, data, p, n, cl, p_val_col) {
  if(!is.null(cl)){
    multidplyr::cluster_library(cl,
                                c("dplyr",
                                  "stringr",
                                  "tidyr",
                                  "purrr",
                                  "tibble"
                                )
    )
    multidplyr::cluster_copy(cl,
                             c(
                               "calc_tpr_fpr",
                               "data",
                               "regex",
                               "p",
                               "n"
                             )
    )
    out <- out %>%
      multidplyr::partition(cl)
  }
  out <- out %>%
    mutate(
      results = map(alpha, calc_tpr_fpr, data, p, n, p_val_col),
      comparison = data$comparison[[1]]
    )
  if(!is.null(cl)){
    out <- out %>%
      collect()
  }
  out %>%
    unnest(results)
}

integrater <- function(x, y) {
  x_na <- is.na(x)
  y_na <- is.na(y)
  x <- x[!x_na&!y_na]
  y <- y[!x_na&!y_na]
  x_order <- order(x)
  x <- x[x_order]
  y <- y[x_order]
  dx <- diff(x)
  end <- length(y)
  my <- (y[1:(end - 1)] + y[2:end]) / 2
  sum(dx * my)
}

#### Testing ####
limma_wrapper <- function(data, design, contrast, weights = NULL){
  data <- data %>%
    dplyr::select(-1) %>%
    as.data.frame() %>%
    magrittr::set_rownames(
      data %>%
        magrittr::use_series(identifier)
    )
  hit <- limma::lmFit(data, design, weights = weights) %>%
    limma::contrasts.fit(contrast) %>%
    limma::eBayes(robust = T, trend = is.null(weights))
  colnames(contrast) %>%
    stats::setNames(colnames(contrast)) %>%
    purrr::map(
      ~limma::topTable(hit, .x, number = Inf, adjust.method = "fdr") %>%
        tibble::rownames_to_column('identifier') %>%
        tibble::as_tibble()
    ) %>%
    purrr::imap(dplyr::mutate) %>%
    purrr::map(dplyr::rename, comparison = dplyr::last_col()) %>%
    dplyr::bind_rows() %>%
    dplyr::select(
      all_of('identifier'), lfc = logFC, mean = AveExpr,
      p_val = P.Value, comparison
    ) %>%
    dplyr::mutate(comparison = stringr::str_replace(comparison, "-", " vs "))
}
calc_weights <- function(df, gm){
  id <- df$identifier
  df %>%
    mutate(
      dplyr::across(where(is.numeric),
                    ~stats::predict.glm(
                      gm,
                      newdata = data.frame(mean = ., c = c),
                      type = "response"
                    ),
      )
    ) %>%
    select(-c,-1,-identifier) %>%
    as.matrix() %>%
    set_rownames(id)
}
ttest_wrapper <- function(contrast, data) {
  comp <- str_split(contrast, '-')[[1]]
  data %>%
    group_by(identifier) %>%
    mutate(
      comparison = str_flatten(comp, ' vs '),
      p_val = t.test(across(contains(comp[1])), across(contains(comp[2])), var.equal = T)$p.value
    ) %>%
    select(identifier, comparison, p_val)
}
baldur_wrapper <- function(data, design, contrast, gam_model, bald_model, workers) {
  uncertainty <- data %>%
    estimate_uncertainty('identifier', design, gam_model)
  data %>%
    estimate_gamma_priors(design, gam_model) %>%
    sample_posterior(
      id_col_name = 'identifier',
      design_matrix = design,
      contrast_matrix = contrast,
      uncertainty_matrix = uncertainty,
      bayesian_model = bald_model,
      clusters = workers
    )
}
#### Figure saving wrapper ####
ggsave_wrapper <- function(file_name, width, height = NA){
  dir.create('figures', F)
  ggsave(paste0('figures/', file_name, '.tiff'), width = width, height = height, units = 'mm', dpi = 320, compression = 'lzw')
  ggsave(paste0('figures/', file_name, '.pdf'), width = width, height = height, units = 'mm', dpi = 320)
}