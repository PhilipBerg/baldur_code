#### Models ####
#sd_semi_single <- rstan::stan_model('stan_c_sd_single_alpha.stan', auto_write = T)
#sd_semi_model_mixed  <- rstan::stan_model('stan_sd_semi_fixed.stan',     auto_write = T)
sd_semi_model_single <- rstan::stan_model('trans_ybar.stan',     auto_write = T)
#sd_semi_trans_mixed <- rstan::stan_model('trans_ybar_mixed.stan',     auto_write = T)
#### Model parameters ####
p_beta <- c(1, 1)*.5
bound  <- .25
prop   <- .01

#### Yeast ####
yeast_prnn <- yeast_prnn %>%
  arrange(sd)

tictoc::tic('sampling')
samp <- fit_lgmr(yeast_prnn, sd_semi_model_single, T, p_beta = p_beta, bound = bound, prop = prop)
tictoc::toc()

p_yeast <- samp %>%
  rstan::summary(pars = 'p') %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

reg <- samp %>%
  rstan::summary(pars = c('I_U', 'I_L', 'S_U', 'S_L')) %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

alpha_yeast <- samp %>%
  rstan::summary(pars = c('alpha')) %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

yeast_prnn <- yeast_prnn %>%
  mutate(
    p = p_yeast,
    alpha = alpha_yeast,
    beta = alpha / reg_fun(p, reg, mean)
  )

yeast_unc <- yeast_prnn %>%
  select(ng50_1:ng100_3, p, mean) %>%
  mutate(
    across(ng50_1:ng100_3, ~ reg_fun(p, reg, .))
  ) %>%
  select(-p, -mean) %>%
  as.matrix() %>%
  set_rownames(yeast_prnn$identifier)


rm(p_yeast, alpha_yeast, samp)

#### UPS ####
ups_prnn <- ups_prnn %>%
  calculate_mean_sd_trends(ups_design) %>%
  arrange(sd)

tictoc::tic('sampling ups')
samp_ups <- fit_lgmr(ups_prnn, sd_semi_model_single, TRUE, p_beta = p_beta, bound = bound, prop = prop)
tictoc::toc()

p_ups <- samp_ups %>%
  rstan::summary(pars = 'p') %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

reg_ups <- samp_ups %>%
  rstan::summary(pars = c('I_U', 'I_L', 'S_U', 'S_L')) %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

alpha_ups <- samp_ups %>%
  rstan::summary(pars = c('alpha')) %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

ups_prnn <- ups_prnn %>%
  mutate(
    p = p_ups,
    alpha = alpha_ups,
    beta = alpha / reg_fun(p, reg_ups, mean)
  )

ups_unc <- ups_prnn %>%
  select(fmol25_1:fmol100_4, p, mean) %>%
  mutate(
    across(fmol25_1:fmol100_4, ~ reg_fun(p, reg_ups, .))
  ) %>%
  select(-p, -mean) %>%
  as.matrix() %>%
  set_rownames(ups_prnn$identifier)

rm(p_ups, alpha_ups, samp_ups)

#### Ramus ####
ramus_prnn <- ramus_prnn %>%
  calculate_mean_sd_trends(ramus_design) %>%
  arrange(sd)

tictoc::tic('ramus')
samp_ramus <- fit_lgmr(ramus_prnn, sd_semi_model_single, 1L, p_beta = p_beta, bound = bound, prop = prop)
tictoc::toc()


p_ramus <- samp_ramus %>%
  rstan::summary(pars = 'p') %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

reg_ramus <- samp_ramus %>%
  rstan::summary(pars = c('I_U', 'I_L', 'S_U', 'S_L')) %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

alpha_ramus <- samp_ramus %>%
  rstan::summary(pars = 'alpha') %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

ramus_prnn <- ramus_prnn %>%
  mutate(
    p = p_ramus,
    alpha = alpha_ramus,
    beta = alpha / reg_fun(p_ramus, reg_ramus, mean)
  )

ramus_unc <- ramus_prnn %>%
  select(condi1_1:condi9_3, p) %>%
  mutate(
    across(condi1_1:condi9_3, ~ reg_fun(p, reg_ramus, .))
  ) %>%
  select(-p) %>%
  as.matrix() %>%
  set_rownames(ramus_prnn$identifier)

rm(p_ramus, alpha_ramus, samp_ramus)

#### Human ####
human_prnn <- human_prnn %>%
  calculate_mean_sd_trends(human_design) %>%
  arrange(sd)

tictoc::tic('human')
samp_human <- fit_lgmr(human_prnn, sd_semi_model_single, 0L, p_beta = p_beta, bound = bound, prop = prop)
tictoc::toc()

p_human <- samp_human %>%
  rstan::summary(pars = 'p') %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

reg_human <- samp_human %>%
  rstan::summary(pars = c('I_U', 'I_L', 'S_U', 'S_L')) %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

alpha_human <- samp_human %>%
  rstan::summary(pars = 'alpha') %>%
  magrittr::use_series(summary) %>%
  magrittr::extract(,1)

human_prnn <- human_prnn %>%
  mutate(
    p = p_human,
    alpha = alpha_human,
    beta = alpha / reg_fun(p, reg_human, mean, FALSE)
  )

human_unc <- human_prnn %>%
  select(spike_prop_25_1:spike_prop_6_23, p, mean) %>%
  mutate(
    across(spike_prop_25_1:spike_prop_6_23, ~ reg_fun(p, reg_human, ., FALSE))
  ) %>%
  select(-p, -mean) %>%
  as.matrix() %>%
  set_rownames(human_prnn$identifier)

rm(p_human, alpha_human, samp_human)
