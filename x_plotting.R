#### Plotting trends ####
p_y <- plot_gamma(yeast_prnn) +
  ggtitle(NULL) +
  theme(
    axis.title.x = element_blank()
  )
p_u <- plot_gamma(ups_prnn) +
  ggtitle(NULL) +
  theme(
    axis.title = element_blank()
  )
p_r <- plot_gamma(ramus_prnn) +
  ggtitle(NULL)
p_h <- plot_gamma(human_prnn) +
  ggtitle(NULL) +
  theme(
    axis.title.y = element_blank()
  )

cowplot::plot_grid(p_y, p_u, p_r, p_h, nrow = 2, labels = 'AUTO', align = 'hv')
ggsave('single_trends.pdf', width = full_page, height = full_page, units = uni, dpi = dpi)

c_y <- plot_gamma_partition(yeast_prnn, yeast_design) +
  ggtitle(NULL) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )
c_u <- plot_gamma_partition(ups_prnn, ups_design) +
  ggtitle(NULL) +
  theme(
    axis.title = element_blank()
  )
c_r <- plot_gamma_partition(ramus_prnn, ramus_design) +
  ggtitle(NULL) +
  theme(
    legend.position = "none"
  )
c_h <- plot_gamma_partition(human_prnn, human_design, formula = sd ~ mean*c) +
  ggtitle(NULL) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank()
  )

cowplot::plot_grid(c_y, c_u, c_r, c_h, nrow = 2, labels = 'AUTO', align = 'hv')
ggsave('partitioned_trends.pdf', width = full_page, height = full_page, units = uni, dpi = dpi)


#### Plot performance ####
color_scheeme <- set_names(viridisLite::turbo(4, end = .9), c("Mix-Baldur", 'Single-Baldur', 'Limma-Trend', 't-test'))
load('performance.RData')
# yeast
yr <- yeast_roc %>%
  filter(between(alpha, 0, 1)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(size = 1/4) +
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
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.1)) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = "none", #c(.775, .4),
    plot.margin = margin(l = 13.5),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm'),
    axis.title.y = element_blank()
  ) +
  labs(
    x = 'False Positive Rate',
    y = ''
  )

yeast_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dashed') +
  geom_path(size = 1/4) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,.6,.1)) +
  theme(
    legend.position = c(.775, .85),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm")
  ) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  )
ggsave_wrapper('yeast_mcc2', width = single, height = single)

yeast_roc %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  mutate(
    intercept = if_else(name == 'precision', 280/1608, 0),
    slope = if_else(name == 'precision', 0, 1),
    name = str_replace(name, 'pre', 'Pre'),
    name = str_replace(name, 'FPR', 'False Positive Rate'),
    name = str_replace(name, 'TPR', 'True Positive Rate')
  ) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, value, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dashed', color = 'red') +
  geom_abline(aes(intercept = intercept, slope = slope), linetype = 'dashed') +
  geom_path(size = 1/4) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.1)) +
  theme(
    legend.position = c(.897, .17),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm")
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(.~name)
ggsave_wrapper('yeast_decomposed2', width = double, height = single)

# UPS
ur <- ups_roc %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(size = 1/4) +
  #lims(x = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,.1), limits = c(0,1)) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  geom_text(data = ups_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  facet_wrap(comparison~., ncol = 1) +
  theme(
    legend.position = c(.8, .85),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.box.margin = margin(),
    legend.margin = margin(),
    legend.background = element_rect(fill = 'white'),
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm'),
    legend.text = element_text(size = 7),
    axis.title.y=element_text(hjust = .25)
  ) +
  labs(
    x = '',
    y = 'True Positive Rate'
  )
cowplot::plot_grid(ur, yr, labels = 'AUTO', rel_heights = c(3, 1.02), hjust = c(-1.5, -1.5), vjust = c(1.1, 0), ncol = 1)
ggsave_wrapper('ups_yeast_roc2', width = single*.95, height = single*2)

ups_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dashed') +
  geom_path(size = 1/4) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.1)) +
  theme(
    legend.position = c(.915, .855),
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
ggsave_wrapper('ups_mcc2', width = double, height = single)


ups_roc %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  mutate(
    intercept = if_else(name == 'precision', 345/10123, 0),
    slope = if_else(name == 'precision', 0, 1),
    name = str_replace(name, 'pre', 'Pre'),
    name = str_replace(name, 'FPR', 'False Positive Rate'),
    name = str_replace(name, 'TPR', 'True Positive Rate')
  ) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, value, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dashed', color = 'red') +
  geom_abline(aes(intercept = intercept, slope = slope), linetype = 'dashed') +
  geom_path(size = 1/4) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  lims(x = c(0,1.02)) +
  scale_y_continuous(breaks = seq(0,1,.1), limits = c(0,1.02)) +
  theme(
    legend.position = c(.89, .75),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm")
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(comparison~name)
ggsave_wrapper('ups_decomposed2', width = double, height = double)

# Ramus
ramus_roc %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(size = 1/4) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = seq(0,1,.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0,1)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 3,
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  geom_text(data = ramus_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  facet_wrap(.~comparison) +
  theme(
    legend.position = 'bottom',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    legend.margin=margin(0,0,.01,0)
  ) +
  labs(
    x = 'False Positive Rate',
    y = 'True Positive Rate'
  )
ggsave_wrapper('ramus_roc', width = double, height = double)

ramus_auroc %>%
  mutate(
    method = factor(method, names(color_scheeme))
  ) %>%
  ggplot(aes(comparison, auROC, fill = method)) +
  geom_col(position = 'dodge') +
  theme_classic() +
  scale_y_continuous(breaks = seq(0,1,.2), labels = seq(0,1,.2), limits = c(0,1)) +
  scale_fill_manual('Method',
                    values = color_scheeme,
                    guide = guide_legend(
                      ncol = 3,
                      title.hjust = .5,
                      override.aes = list(size = 1)
                    )
  ) +
  facet_wrap(.~comparison, scales = 'free_x') +
  theme(
    legend.position = 'bottom',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.margin=margin(0,0,0.01,0),
    panel.grid.major.y =  element_line('grey', .1)
  )
ggsave_wrapper('ramus_auroc', width = double, height = double)

ramus_roc %>%
  arrange(alpha) %>%
  filter(between(alpha, 0, 1)) %>%
  ggplot(aes(alpha, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(size = 1/4) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = seq(0,1,.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0,1)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 3,
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  geom_text(data = ramus_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  facet_wrap(.~comparison) +
  theme(
    legend.position = 'bottom',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    legend.margin=margin(0,0,.01,0)
  ) +
  labs(
    x = expression(alpha),
    y = 'True Positive Rate'
  )

ramus_roc %>%
  arrange(alpha) %>%
  filter(between(alpha, 0, 1)) %>%
  ggplot(aes(alpha, FPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(size = 1/4) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = seq(0,1,.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0,1)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 3,
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  facet_wrap(.~comparison) +
  theme(
    legend.position = 'bottom',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    legend.margin=margin(0,0,.01,0)
  ) +
  labs(
    x = expression(alpha),
    y = 'True Positive Rate'
  )

# Human
human_roc %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(size = 1/4) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = seq(0,1,.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0,1)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 3,
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  geom_text(data = human_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  facet_wrap(.~comparison) +
  theme(
    legend.position = 'bottom',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    legend.margin=margin(0,0,.01,0)
  ) +
  labs(
    x = 'False Positive Rate',
    y = 'True Positive Rate'
  )

human_roc %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  mutate(
    intercept = if_else(name == 'precision', 345/10123, 0),
    slope = if_else(name == 'precision', 0, 1),
    name = str_replace(name, 'pre', 'Pre'),
    name = str_replace(name, 'FPR', 'False Positive Rate'),
    name = str_replace(name, 'TPR', 'True Positive Rate')
  ) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, value, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dashed', color = 'red') +
  geom_abline(aes(intercept = intercept, slope = slope), linetype = 'dashed') +
  geom_path(size = 1/4) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.hjust = .5,
                       override.aes = list(size = 1)
                     )
  ) +
  lims(x = c(0,1.02)) +
  scale_y_continuous(breaks = seq(0,1,.1), limits = c(0,1.02)) +
  theme(
    legend.position = c(.89, .75),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm")
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(comparison~name)
