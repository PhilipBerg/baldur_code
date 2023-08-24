#### Yeast ####
tictoc::tic('sampling')
samp <- fit_lgmr(yeast_prnn, 'identifier', simplify = TRUE, cores = 5)
tictoc::toc()

yeast_pars <- coef(samp, TRUE, c('coef', 'aux'))

#### UPS ####
ups_prnn <- ups_prnn %>%
  calculate_mean_sd_trends(ups_design)

tictoc::tic('sampling ups')
samp_ups <- baldur::fit_lgmr(ups_prnn, 'identifier', simplify = TRUE, cores = 5)
tictoc::toc()

ups_pars <- coef(samp_ups, TRUE, c('coef', 'aux'))

#### Ramus ####
tictoc::tic('ramus')
samp_ramus <- fit_lgmr(ramus_prnn, 'identifier', simplify = TRUE, cores = 5)
tictoc::toc()

ramus_pars <- coef(samp_ramus, TRUE, c('coef', 'aux'))

#### Human ####
tictoc::tic('human')
samp_human <- fit_lgmr(human_prnn, 'identifier', simplify = TRUE, cores = 5)
tictoc::toc()

human_pars <- coef(samp_human, TRUE, c('coef', 'aux'))

#### Bruderer ####
tictoc::tic('Bruderer')
samp_bruder <- fit_lgmr(bruder_prnn, 'identifier', simplify = TRUE, cores = 5)
tictoc::toc()

human_bruder <- coef(samp_bruder, TRUE, c('coef', 'aux'))
