#### Plot trends ####
p_y <- {
  plot_gamma(yeast_prnn) +
    ggtitle(NULL)
} %>%
  set_shape()
p_u <- {
  plot_gamma(ups_prnn) +
    ggtitle(NULL) +
    theme(
      axis.title.y = element_blank()
    )
} %>%
  set_shape()
p_r <- {
  plot_gamma(ramus_prnn) +
    ggtitle(NULL)+
    theme(
      axis.title.y = element_blank()
    )
} %>%
  set_shape()
p_h <- {
  plot_gamma(human_prnn) +
    ggtitle(NULL) +
    theme(
      axis.title.y = element_blank()
    )
} %>%
  set_shape()
p_b <- {
  plot_gamma(bruder_prnn) +
    ggtitle(NULL) +
    theme(
      axis.title.y = element_blank()
    )
} %>%
  set_shape()

cowplot::plot_grid(p_y, p_u, p_r, p_h, p_b, nrow = 1, labels = 'AUTO', align = 'hv')
ggsave_wrapper('single_trends', width = full_page, height = full_page/3)

rng <- 10
n    <- 100
f_y <- plot_regression_field(samp, rng = rng, n = n) +
  theme(
    legend.background = element_blank(),
    legend.key.size = unit(8, 'pt'),
    legend.text = element_text(size = unit(9, 'pt')),
    legend.position = c(.8, .8),
    legend.spacing.y = unit(0, 'pt'),
    legend.spacing.x = unit(0, 'pt')
  )

f_u <- plot_regression_field(samp_ups, rng = rng, n = n) +
    theme(
      legend.background = element_blank(),
      axis.title.y = element_blank(),
      legend.key.size = unit(8, 'pt'),
      legend.text = element_text(size = unit(9, 'pt')),
      legend.position = c(.8, .8),
      legend.spacing.y = unit(0, 'pt'),
      legend.spacing.x = unit(0, 'pt')
    )

f_r <- plot_regression_field(samp_ramus, rng = rng, n = n) +
  ggtitle(NULL) +
  theme(
    legend.background = element_blank(),
    axis.title.y = element_blank(),
    legend.key.size = unit(8, 'pt'),
    legend.text = element_text(size = unit(9, 'pt')),
    legend.position = c(.8, .8),
    legend.spacing.y = unit(0, 'pt'),
    legend.spacing.x = unit(0, 'pt')
  )

f_h <- plot_regression_field(samp_human, rng = rng, n = n) +
  theme(
    legend.background = element_blank(),
    legend.position = c(.8, .8),
    legend.key.size = unit(8, 'pt'),
    legend.text = element_text(size = unit(9, 'pt')),
    axis.title.y = element_blank(),
    legend.spacing.y = unit(0, 'pt'),
    legend.spacing.x = unit(0, 'pt'),
    plot.margin = margin(l = .5)
  )

f_b <- plot_regression_field(samp_bruder, rng = rng, n = n) +
  theme(
    legend.background = element_blank(),
    legend.position = c(.8, .8),
    legend.key.size = unit(8, 'pt'),
    legend.text = element_text(size = unit(9, 'pt')),
    axis.title.y = element_blank(),
    legend.spacing.y = unit(0, 'pt'),
    legend.spacing.x = unit(0, 'pt'),
    plot.margin = margin(l = .5)
  )

paste0('f_', c('y', 'u', 'r', 'h', 'b')) %>%
  map(get) %>%
  map(
    ~ {.x$layers[[1]]$aes_params$size <- .05; .x}
  ) %>%
  cowplot::plot_grid(plotlist = ., labels = 'AUTO', align = 'hv', nrow = 1)
ggsave_wrapper('field_mixed', width = full_page, height = full_page*1/3)

#### Plot performance ####
color_scheeme <- set_names(viridisLite::turbo(6, end = .9),
                           c("LGMR-Baldur EB", "LGMR-Baldur WI", 'GR-Baldur EB', 'GR-Baldur WI',
                                                              'Limma-Trend', 't-test'))
# yeast
yeast_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  geom_text(data = yeast_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.7, .3),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm')
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  )
ggsave_wrapper('yeast_roc', half_page, half_page)

yeast_roc %>%
  filter(between(alpha, 0, .2)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,.6,.2)) +
  scale_x_continuous(breaks = c(.01, .05, .1, .15, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.795, .16),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.spacing = unit(0, 'mm')
  ) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  )
ggsave_wrapper('yeast_mcc', width = half_page, height = half_page)

yeast_roc %>%
  filter(between(alpha, 0, .2)) %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  mutate(
    intercept = if_else(name == 'precision', sum(str_detect(yeast_prnn$identifier, 'YEAST'))/nrow(yeast_prnn), 0),
    slope = if_else(name == 'precision', 0, 1),
    name = str_replace(name, 'pre', 'Pre'),
    name = str_replace(name, 'FPR', 'False Positive Rate'),
    name = str_replace(name, 'TPR', 'True Positive Rate'),
    name = factor(name, levels = facet_order, ordered = T)
  ) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, value, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = c(.01, .05, .1, .15, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.76, .25),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.25, "cm"),
    legend.background = element_blank()
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(name~.)
ggsave_wrapper('yeast_decomposed', width = half_page, height = full_page)

yeast_roc %>%
  group_by(method) %>%
  filter(alpha <= .05) %>%
  slice_max(alpha) %>%
  pivot_longer(c(TP, FP)) %>%
  mutate(
    method = factor(method, levels = rev(unique(yeast_auroc$method))[c(1:2, 4:3, 6:5)])
  ) %>%
  ggplot(aes(method, value, fill = name, label = value)) +
  geom_col(position = position_dodge2(1, padding = .001, reverse = T)) +
  geom_text(aes(method, value+25, color = name, label = value), show.legend = F, position = position_dodge2(1, reverse = T), size = 3, hjust = 0) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1250), breaks = seq(0, 1100, by = 300)) +
  labs(
    y = 'Count', x = 'Method', fill = 'Category'
  ) +
  scale_color_viridis_d(option = 'H', direction = -1) +
  scale_fill_viridis_d(option = 'H', direction = -1) +
  theme(
    legend.position = c(0.8,.875),
    axis.text.x = element_text(size = 7),
    legend.background = element_blank(),
    plot.margin = margin()
  )
ggsave_wrapper('bar_alpha_05_yeast', width = half_page, height = half_page)

# UPS
ur <- ups_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_path(linewidth = 1/6) +
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  geom_text(data = ups_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 2.5,
            show.legend = F
  ) +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       override.aes = list(linewidth = 1)
                     )
  ) +
  facet_wrap(
    factor(
      comparison,
      levels = c(
        paste0('fmol', c(25, 50, 25), ' vs ', 'fmol', c(100, 100, 50)),
        "fmol25 vs fmol50 and fmol100", "fmol100 vs fmol25 and fmol50",
        "fmol50 vs fmol25 and fmol100"
      )
    )~., ncol = 3, dir = 'h') +
  theme(
    legend.position = "none",
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.box.margin = margin(),
    legend.margin = margin(),
    legend.background = element_rect(fill = 'white'),
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm'),
    legend.text = element_text(size = 7),
    axis.title = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  )

ups_roc %>%
  filter(between(alpha, 0, .3)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       #title.hjust = .5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = c(.01, .05, .1, .15, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.92, .2),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm')
  ) +
  facet_wrap(.~factor(
    comparison,
    levels = c(
      paste0('fmol', c(25, 50, 25), ' vs ', 'fmol', c(100, 100, 50)),
      "fmol25 vs fmol50 and fmol100", "fmol100 vs fmol25 and fmol50",
      "fmol50 vs fmol25 and fmol100"
    )
  ), nrow = 3, dir = 'v'
  ) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  ) +
  coord_cartesian(xlim = c(0, .2))
ggsave_wrapper('ups_mcc', width = full_page, height = full_page)


ups_roc %>%
  filter(between(alpha, 0, .3) & str_detect(comparison, 'and')) %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  group_by(comparison, method) %>%
  mutate(
    name = str_replace_all(name,
                           c(
                             'FPR' = 'False Positive Rate',
                             'TPR' = 'True Positive Rate',
                             'precision' = 'Precision'
                           )
    )
  ) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, value, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                         title.position = "top",
                         ncol = 1,
                         override.aes = list(linewidth = 1.5)
                     )
  ) +
  scale_y_continuous(position="right", labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = c(.01, .05, .1, .15, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.9, .925),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm'),
    legend.background = element_blank()
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(
    factor(name, levels = facet_order)~factor(
    comparison,
    levels = c(
      paste0('fmol', c(25, 50, 25), ' vs ', 'fmol', c(100, 100, 50)),
      "fmol25 vs fmol50 and fmol100", "fmol100 vs fmol25 and fmol50",
      "fmol50 vs fmol25 and fmol100"
    )
  ), scales = 'free', switch = 'y') +
  coord_cartesian(xlim = c(0, .2))
ggsave_wrapper('ups_decomposed_complex', width = full_page, height = full_page)

ups_roc %>%
  filter(between(alpha, 0, .3) & str_detect(comparison, 'and', T)) %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  group_by(comparison, method) %>%
  mutate(
    name = str_replace_all(name,
                           c(
                             'FPR' = 'False Positive Rate',
                             'TPR' = 'True Positive Rate',
                             'precision' = 'Precision'
                           )
    )
  ) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, value, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       title.position = "top",
                       ncol = 1,
                       override.aes = list(linewidth = 1.5)
                     )
  ) +
  scale_y_continuous(position="right", labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = c(.01, .05, .1, .15, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.9, .925),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm'),
    legend.background = element_blank()
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(
    factor(name, levels = facet_order)~factor(
      comparison,
      levels = c(
        paste0('fmol', c(25, 50, 25), ' vs ', 'fmol', c(100, 100, 50)),
        "fmol25 vs fmol50 and fmol100", "fmol100 vs fmol25 and fmol50",
        "fmol50 vs fmol25 and fmol100"
      )
    ), scales = 'free', switch = 'y') +
  coord_cartesian(xlim = c(0, .2))
ggsave_wrapper('ups_decomposed', width = full_page, height = full_page)

ups_roc %>%
  group_by(method, comparison) %>%
  filter(alpha <= .05) %>%
  slice_max(alpha) %>%
  arrange(TP) %>%
  pivot_longer(c(TP, FP)) %>%
  mutate(
    method = factor(method, levels = rev(unique(ups_auroc$method))[c(1:2, 4:3, 5:6)])
  ) %>%
  ggplot(aes(method, value, fill = name, label = value)) +
  geom_col(position = position_dodge2(1, padding = .01, reverse = T)) +
  geom_text(aes(method, value+5, color = name, label = value),
            show.legend = F,
            position = position_dodge2(1, padding = .01, reverse = T),
            size = 2.9, hjust = 0
  ) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 475)) +
  labs(
    y = 'Count', x = 'Method', fill = 'Category'
  ) +
  scale_color_viridis_d(option = 'H', direction = -1) +
  scale_fill_viridis_d(option = 'H', direction = -1,
                       guide = guide_legend(
                         title.position = "top",
                       )) +
  theme(
    legend.position = c(0.91, .2),
    axis.text.x = element_text(size = 7),
    legend.background = element_blank(),
    legend.direction = 'vertical',
    strip.text = element_text(size = 5, face = 'bold')
  ) +
  facet_wrap(factor(
    comparison,
    levels = c(
      paste0('fmol', c(25, 50, 25), ' vs ', 'fmol', c(100, 100, 50)),
      "fmol25 vs fmol50 and fmol100", "fmol100 vs fmol25 and fmol50",
      "fmol50 vs fmol25 and fmol100"
    )
  )~., scales = 'free_y')
ggsave_wrapper('bar_alpha_05_ups', width = full_page, height = half_page)

# Ramus
ramus_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  mutate(
    shape = case_when(
      (alpha - .05)^2 == min((alpha - .05)^2) ~ '0.05',
      (alpha - .01)^2 == min((alpha - .01)^2) ~ '0.01',
      (alpha - .001)^2 == min((alpha - .001)^2) ~ '0.001',
      T ~ 'F'
    ),
    comparison = factor(comparison, ramus_plot_order)
  ) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_path(linewidth = 1/6) +
  scale_y_continuous(breaks = seq(0,1,.5), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.5), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       nrow = 2,
                       title.hjust = .5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  facet_wrap(.~comparison) +
  theme(
    legend.position = 'bottom',
    plot.margin = margin(),
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    legend.margin=margin(0,0,.01,0),
    strip.text = element_text(size = 8)
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  )
ggsave_wrapper('ramus_roc', width = full_page, height = full_page)

ramus_roc %>%
  filter(between(alpha, 0, .2)) %>%
  mutate(
    comparison = factor(comparison, ramus_plot_order)
  ) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       nrow = 1,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1,.25), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = c(0, .05, .125, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = 'top',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 8)
  ) +
  facet_wrap(.~comparison) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  )
ggsave_wrapper('ramus_mcc', width = full_page, height = full_page)

ramus_roc %>%
  mutate(
    comparison = factor(comparison, ramus_plot_order)
  ) %>%
  arrange(alpha) %>%
  filter(between(alpha, 0, .2)) %>%
  ggplot(aes(alpha, TPR, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  scale_y_continuous(breaks = seq(0,1,.5), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = c(0, .05, .125, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       nrow = 1,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  facet_wrap(.~comparison) +
  theme(
    legend.position = 'top',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    legend.margin=margin(0,0,.01,0),
    strip.text = element_text(size = 8)
  ) +
  labs(
    x = expression(alpha),
    y = 'True Positive Rate'
  )
ggsave_wrapper('ramus_tpr', width = full_page, height = full_page)

ramus_roc %>%
  arrange(alpha) %>%
  filter(between(alpha, 0, .2)) %>%
  mutate(
    comparison = factor(comparison, ramus_plot_order)
  ) %>%
  ggplot(aes(alpha, FPR, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  scale_y_continuous(breaks = seq(0,1,.5), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = c(0, .05, .125, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       nrow = 1,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  facet_wrap(.~comparison) +
  theme(
    legend.position = 'top',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    legend.margin=margin(0,0,.01,0),
    strip.text = element_text(size = 8)
  ) +
  labs(
    x = expression(alpha),
    y = 'False Positive Rate'
  )
ggsave_wrapper('ramus_fpr', width = full_page, height = full_page)

ramus_roc %>%
  group_by(method, comparison) %>%
  filter(alpha <= .05) %>%
  slice_max(alpha) %>%
  mutate(
    comparison = factor(comparison, ramus_plot_order)
  ) %>%
  arrange(FP) %>%
  pivot_longer(c(TP, FP)) %>%
  mutate(
    name = str_replace(name, 'P', 'Positive'),
    name = str_replace_all(name, c('T' = 'True ', 'F' = 'False ')),
    method = factor(method, levels = rev(unique(ramus_auroc$method))[c(2:1, 4:3, 6:5)])
  ) %>%
  ggplot() +
  geom_col(aes(method, value, fill = name), position = position_dodge2(1, reverse = T)) +
  geom_text(aes(method, value+25, color = name, label = value), position = position_dodge2(1, reverse = T), hjust = 0, size = 1.8, show.legend = F) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 3100, by = 1500), limits = c(0, 3500)) +
  labs(
    y = 'Count', x = 'Method', fill = 'Category'
  ) +
  theme(
    legend.position = 'bottom',
    axis.text = element_text(size = 7.5),
    strip.text = element_text(size = 8),
    legend.box.margin = margin(),
    legend.margin = margin()
  ) +
  coord_flip() +
  scale_color_viridis_d(option = 'H', direction = -1) +
  scale_fill_viridis_d(option = 'H', direction = -1) +
  facet_wrap(comparison~.)
ggsave_wrapper('bar_alpha_05_ramus', width = full_page, height = full_page)

# Human
hr <- human_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  filter(between(alpha, 0, 1)) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_path(linewidth = 1/6) +
  scale_y_continuous(breaks = seq(0, 1, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0, 1, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       override.aes = list(linewidth = 1)
                     )
  ) +
  geom_text(data = human_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 2.5,
            show.legend = F
  ) +
  facet_wrap(factor(
    comparison,
    levels = paste0('1:', c(25, 12, 25), ' vs ', '1:', c(6, 6, 12))
  )~., nrow = 1) +
  theme(
    legend.position = c(.9, .36),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.box.margin = margin(),
    legend.margin = margin(),
    legend.background = element_rect(fill = 'white'),
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm'),
    axis.title = element_blank()
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  )

cowplot::plot_grid(ur, hr, labels = 'AUTO', nrow = 2, rel_heights = c(.6, .3), rel_widths = c(1, 1)) + #, hjust = c(-1.9, -1.8), rel_widths = c(1, 1)) +
  theme(
    plot.margin = margin(0, 0, 15, 15)
  ) +
  cowplot::draw_label("False Positive Rate", x = 0.5, y = 0,   angle = 0,  vjust = 1, hjust = 0.45, size = 11) +
  cowplot::draw_label("True Positive Rate",  x = 0,   y = 0.5, angle = 90, vjust = -.35, hjust = .52, size = 11)
ggsave_wrapper('ups_human_roc', width = full_page, height = full_page)

human_roc %>%
  filter(between(alpha, 0, .2)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.hjust = .5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = c(.01, .05, .1, .15, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.92, .16),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7)
  ) +
  facet_grid(.~comparison) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  )
ggsave_wrapper('human_mcc', width = full_page, height = half_page)

human_roc %>%
  filter(between(alpha, 0, .2)) %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  mutate(
    intercept = if_else(name == 'precision', sum(str_detect(human_prnn$identifier, 'ECOLI'))/nrow(human_prnn), 0),
    slope = if_else(name == 'precision', 0, 1),
    name = str_replace(name, 'pre', 'Pre'),
    name = str_replace(name, 'FPR', 'False Positive Rate'),
    name = str_replace(name, 'TPR', 'True Positive Rate')
  ) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, value, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       nrow = 2,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(position="right") +
  scale_x_continuous(breaks = c(.01, .05, .1, .15, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = 'top',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.margin = margin()
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  )  +
  facet_grid(factor(name, levels = facet_order)~comparison, scales = 'free', switch = 'y')
ggsave_wrapper('human_decomp', width = full_page, height = full_page)

human_roc %>%
  group_by(method, comparison) %>%
  filter(alpha <= .05) %>%
  slice_max(alpha) %>%
  pivot_longer(c(TP, FP)) %>%
  mutate(
    name = str_replace(name, 'P', 'Positive'),
    name = str_replace_all(name, c('T' = 'True ', 'F' = 'False ')),
    method = factor(method, levels = rev(unique(human_auroc$method))[c(4:3, 1:2, 5:6)])
  ) %>%
  ggplot(aes(method, value, fill = name, label = value)) +
  geom_col(position = position_dodge2(1, reverse = T)) +
  geom_text(aes(method, value+25, color = name, label = value), position = position_dodge2(1, reverse = T), hjust = 0, size = 1.9, show.legend = F) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 3100, by = 200), limits = c(0, 500)) +
  labs(
    y = 'Count', x = 'Method', fill = 'Category'
  ) +
  theme(
    legend.position = 'bottom',
    axis.text = element_text(size = 7.5),
    strip.text = element_text(size = 8),
    legend.box.margin = margin(),
    legend.margin = margin()
  ) +
  coord_flip() +
  scale_color_viridis_d(option = 'H', direction = -1) +
  scale_fill_viridis_d(option = 'H', direction = -1) +
  facet_wrap(comparison~.)
ggsave_wrapper('bar_alpha_05_human', width = full_page, height = half_page)

# Bruderer
bruder_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  geom_text(data = bruder_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.7, .3),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm')
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  )
ggsave_wrapper('bruderer_roc', half_page, half_page)

bruder_roc %>%
  filter(between(alpha, 0, .2)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_segment(aes(x = .05 , y = 0, xend = .05, yend = .62), linetype = 'dotted', color = 'red') +
  # geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 2,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,.6,.2)) +
  scale_x_continuous(breaks = c(.01, .05, .1, .15, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.6, .94),
    plot.margin = margin(),
    legend.background = element_blank(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.25, "cm"),
    legend.spacing = unit(0, 'mm')
  ) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  )
ggsave_wrapper('bruderer_mcc', width = half_page, height = half_page)

bruder_roc %>%
  filter(between(alpha, 0, .2)) %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  mutate(
    intercept = if_else(name == 'precision', sum(str_detect(yeast_prnn$identifier, 'YEAST'))/nrow(yeast_prnn), 0),
    slope = if_else(name == 'precision', 0, 1),
    name = str_replace(name, 'pre', 'Pre'),
    name = str_replace(name, 'FPR', 'False Positive Rate'),
    name = str_replace(name, 'TPR', 'True Positive Rate'),
    name = factor(name, levels = facet_order, ordered = T)
  ) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, value, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = c(.01, .05, .1, .15, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.9, .8),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.25, "cm"),
    legend.background = element_blank()
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(.~name)
ggsave_wrapper('bruderer_decomposed', width = full_page, height = half_page)

bruder_roc %>%
  group_by(method) %>%
  filter(alpha <= .05) %>%
  slice_max(alpha) %>%
  pivot_longer(c(TP, FP)) %>%
  mutate(
    method = factor(method, levels = rev(unique(bruder_auroc$method))[c(1:2, 4:3, 6:5)])
  ) %>%
  ggplot(aes(method, value, fill = name, label = value)) +
  geom_col(position = position_dodge2(1, padding = .001, reverse = T)) +
  geom_text(aes(method, value+25, color = name, label = value), show.legend = F, position = position_dodge2(1, reverse = T), size = 3, hjust = 0) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7500), breaks = seq(0, 7500, by = 1000)) +
  labs(
    y = 'Count', x = 'Method', fill = 'Category'
  ) +
  scale_color_viridis_d(option = 'H', direction = -1) +
  scale_fill_viridis_d(option = 'H', direction = -1) +
  theme(
    legend.position = c(0.8,.875),
    axis.text.x = element_text(size = 7),
    legend.background = element_blank(),
    plot.margin = margin()
  )
ggsave_wrapper('bar_alpha_05_bruderer', width = half_page, height = half_page)


#### Power ####
load('human_3_columns.RData')
results_concat <- results

load('human_6_columns.RData')
results_concat <- results %>%
  bind_rows(results_concat)

load('human_9_columns.RData')
results_concat <- results %>%
  bind_rows(results_concat)

load('human_12_columns.RData')
results_concat <- results %>%
  bind_rows(results_concat)

load('human_15_columns.RData')
results_concat <- results %>%
  bind_rows(results_concat)

load('human_18_columns.RData')
results_concat <- results %>%
  bind_rows(results_concat)

load('human_21_columns.RData')
results_concat <- results %>%
  bind_rows(results_concat)

results_concat <- human_roc %>%
  mutate(
    columns = 23
  ) %>%
  bind_rows(results_concat)

results_concat %>%
  group_by(comparison, method, columns, rep) %>%
  filter(alpha <= .05) %>%
  slice_max(alpha) %>%
  pivot_longer(c(TP, FP)) %>%
  group_by(comparison, method, columns, name) %>%
  summarise(
    mean = mean(value),
    se   = sd(value)/sqrt(n()),
    ub   = mean + se,
    lb   = mean - se
  ) %>%
  mutate(
    name = str_replace_all(name, c('TP' = 'True Positives', 'FP' = 'False Positives')),
  ) %>%
  ggplot(aes(factor(columns), mean, fill = method, ymin = lb, ymax = ub)) +
  geom_col(position = position_dodge(1)) +
  geom_errorbar(position = position_dodge(1), color = 'grey50') +
  scale_fill_manual('Method',
                    values = color_scheeme,
                    guide = guide_legend(
                      nrow = 2,
                      title.position = "top",
                    )
  ) +
  facet_grid(factor(comparison, levels = unique(comparison)[c(2, 1, 3)])~name, drop = F, as.table = F) +
  theme_classic() +
  theme(
    legend.position = 'bottom',
    legend.margin = margin()
  ) +
  labs(
    y = 'Counts',
    x = 'Replicates'
  )
ggsave_wrapper('tpfp_subsets', width = full_page, height = full_page*2/3)

#### Running time ####
color_data <- set_names(
  viridisLite::turbo(5, end = .9), c('Yeast', 'UPS', 'Ramus', 'Human', 'Bruder')
)
tdd <- time_data_dec_model %>%
  ggplot(aes(workers, time/60)) +
  stat_summary(aes(color = dataset), fun = median, geom = "line", linewidth = .1) +
  geom_boxplot(aes(fill = dataset, group = interaction(workers, dataset)), position = position_dodge(width = 0),
               outlier.shape = '.', size = .1) +
  # geom_text(
  #   inherit.aes = F,
  #   data = summarise(time_data_dec_model, time = compose(round, median)(time/60), .by = c(dataset, workers)),
  #   aes(workers + .75, time + 1, label = time, color = dataset), size = 1.5
  # ) +
  scale_x_continuous(breaks = c(seq(2, 10, 2), 16, 32, 64)) +
  theme_classic() +
  theme(
    legend.position = c(.5, .8),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    panel.grid.major.y = element_line(linewidth = .5)
  ) +
  labs(
    x = "Parallel Workers", y = 'Time (Min)', color = 'Dataset', fill = 'Dataset'
  ) +
  scale_color_manual(values = color_data) +
  scale_fill_manual(values = color_data)

lt <- lgmr_time %>%
  ggplot(aes(p, time)) +
  geom_smooth(linewidth = .3, se = F, method = lm, formula = y ~ x) +
  geom_boxplot(aes(fill = dataset, group = interaction(p, dataset)), outlier.shape = '.', size = .1) +
  # geom_text(
  #   inherit.aes = F,
  #   data = summarise(lgmr_time, time = compose(round, median)(time), .by =  c(p, dataset)),
  #   aes(p + 325, time - .25, label = time, color = dataset), size = 3
  # ) +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.grid.major.y = element_line(linewidth = .5)
  ) +
  scale_x_continuous(breaks = unique(lgmr_time$p)) +
  scale_color_manual("Dataset", values = color_data) +
  scale_fill_manual("Dataset", values = color_data) +
  labs(
    y = 'Time (min)', x = "Number of Peptides"
  )
cowplot::plot_grid(tdd, lt, nrow = 2, labels = "AUTO")
ggsave_wrapper('running_time', width = full_page, height = full_page*1/2)

#### Tables ####
mget(rev(ls(pattern = '_pars'))) %>%
  map(
    ~ c(
      .x$aux[1], .x$coef, .x$aux[2]
    )
  ) %>%
  map(round, 3) %>%
  imap(
    ~ c(.y, .x)
  ) %>%
  map(
    str_remove, '_pars'
  ) %>%
  map(
    str_to_title
  ) %>%
  map(
    str_replace, 'Ups', 'UPS'
  ) %>%
  map_chr(str_flatten, '\t&\t') %>%
  c(., '') %>%
  cat(sep = '\\\\\n')

mget(ls(pattern = '_sin_gam'))  %>%
  map(MASS::gamma.shape) %>%
  map_dbl(use_series, alpha) %>%
  map(round, 2)

mget(ls(pattern = '_sin_gam')) %>%
  map(coef) %>%
  map(round, 2) %>%
  map_chr(str_flatten, ' & ')

mget(ls(pattern = '_sin_gam')) %>%
  map_dbl(calc_nrmse) %>%
  map(round, 3)
