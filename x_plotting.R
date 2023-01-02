#### Plot trends ####
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
ggsave_wrapper('single_trends', width = full_page, height = full_page)

c_y <- plot_regression(yeast_prnn, reg, mixed = FALSE) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )
c_u <- plot_regression(ups_prnn, reg_ups, mixed = FALSE) +
  theme(
    axis.title = element_blank()
  )
c_r <- plot_regression(ramus_prnn, reg_ramus, mixed = FALSE) +
  ggtitle(NULL) +
  theme(
    legend.position = "none"
  )
c_h <- plot_regression(human_prnn, reg_human, mixed = FALSE) +
  ggtitle(NULL) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank()
  )

cowplot::plot_grid(c_y, c_u, c_r, c_h, nrow = 2, labels = 'AUTO', align = 'hv')
ggsave_wrapper('partitioned_trends_mixed', width = full_page, height = full_page)


#### Plot performance ####
color_scheeme <- set_names(viridisLite::turbo(4, end = .9), c("LGMR-Baldur", 'GR-Baldur',
                                                              'Limma-Trend', 't-test'))
# yeast
yeast_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  mutate(
    shape = case_when(
      (alpha - .05)^2 == min((alpha - .05)^2) ~ '0.05',
      (alpha - .01)^2 == min((alpha - .01)^2) ~ '0.01',
      (alpha - .001)^2 == min((alpha - .001)^2) ~ '0.001',
      T ~ 'F'
    )
  ) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
  theme_classic() +
  geom_text(data = yeast_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  geom_point(
    data = \(x) subset(x, shape != 'F'),
    aes(FPR, TPR, color = method, shape = shape)
  ) +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       #title.hjust = .5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.775, .4),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
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
  filter(between(alpha, 0, 1)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       #title.hjust = .5,
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,.6,.2)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.775, .85),
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
  filter(between(alpha, 0, 1)) %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  mutate(
    intercept = if_else(name == 'precision', sum(str_detect(yeast_prnn$identifier, 'YEAST'))/nrow(yeast_prnn), 0),
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
  geom_path(linewidth = 1/4) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       #title.hjust = .5,
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
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
ggsave_wrapper('yeast_decomposed', width = full_page, height = half_page)

yeast_roc %>%
  group_by(method) %>%
  filter(min((alpha - .05)^2) == (alpha - .05)^2) %>%
  pivot_longer(c(TP, FP)) %>%
  mutate(
    name = str_replace_all(name, c('T' = 'True ', 'F' = 'False ', 'P' = 'Positive'))
  ) %>%
  ggplot(aes(method, value, fill = name, label = value)) +
  geom_col(position = position_dodge(1)) +
  geom_text(position = position_dodge(1), vjust = 1) +
  theme_classic() +
  scale_y_continuous(expand = c(0,1)) +
  labs(
    y = 'Count', x = 'Method', fill = 'Category'
  ) +
  theme(
    legend.position = c(0.2,.8),
    axis.text.x = element_text(size = 7)
  )
ggsave_wrapper('bar_alpha_05_yeast', width = half_page, height = half_page)

# UPS
ur <- ups_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  mutate(
    shape = case_when(
      (alpha - .05)^2 == min((alpha - .05)^2) ~ '0.05',
      (alpha - .01)^2 == min((alpha - .01)^2) ~ '0.01',
      (alpha - .001)^2 == min((alpha - .001)^2) ~ '0.001',
      T ~ 'F'
    )
  ) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  geom_text(data = ups_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  geom_point(
    data = \(x) subset(x, shape != 'F'),
    aes(FPR, TPR, color = method, shape = shape)
  ) +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       #title.hjust = .5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  facet_wrap(
    factor(
      comparison,
      levels = paste0('fmol', c(25, 50, 25), ' vs ', 'fmol', c(100, 100, 50))
    )~., ncol = 1) +
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
    axis.title = element_blank()
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  )
#cowplot::plot_grid(ur, yr, labels = 'AUTO', rel_heights = c(3, 1.02), hjust = c(-1.5, -1.5), vjust = c(1.1, 0), ncol = 1)
#ggsave_wrapper('ups_yeast_roc', width = half_page, height = half_page*4)

ups_roc %>%
  filter(between(alpha, 0, 1)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
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
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.915, .855),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm')
  ) +
  facet_grid(.~comparison) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  )
ggsave_wrapper('ups_mcc', width = full_page, height = half_page)


ups_roc %>%
  filter(between(alpha, 0, 1)) %>%
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
  geom_path(linewidth = 1/4) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       #title.hjust = .5,
                       override.aes = list(linewidth = 1.5)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.89, .75),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm')
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(comparison~name)
ggsave_wrapper('ups_decomposed', width = full_page, height = full_page)

ups_roc %>%
  group_by(method, comparison) %>%
  filter(min((alpha - .05)^2) == (alpha - .05)^2) %>%
  pivot_longer(c(TP, FP)) %>%
  ggplot(aes(method, value, fill = name, label = value)) +
  geom_col(position = position_dodge(1)) +
  geom_text(position = position_dodge(1), vjust = .5) +
  theme_classic() +
  scale_y_continuous(expand = c(0,20)) +
  labs(
    y = 'Count', x = 'Method', fill = 'Category'
  ) +
  theme(
    legend.position = c(0.1,.9),
    axis.text.x = element_text(size = 7)
  ) +
  facet_grid(comparison~.)
ggsave_wrapper('bar_alpha_05_ups', width = half_page, height = full_page)

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
    )
  ) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
  geom_point(
    data = \(x) subset(x, shape != 'F'),
    aes(FPR, TPR, color = method, shape = shape)
  ) +
  scale_y_continuous(breaks = seq(0,1,.5), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.5), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       nrow = 1,
                       title.hjust = .5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  # geom_text(data = ramus_auroc,
  #           aes(FPR, TPR, label = round(auROC, 3)),
  #           size = 3,
  #           show.legend = F
  # ) +
  facet_wrap(.~comparison) +
  theme(
    legend.position = 'bottom',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    legend.margin=margin(0,0,.01,0),
    strip.text = element_text(size = 7)
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  )
ggsave_wrapper('ramus_roc', width = full_page, height = full_page)

ramus_roc %>%
  filter(between(alpha, 0, 1)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
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
  scale_y_continuous(breaks = seq(0,1,.25), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.25), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.915, .855),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7)
  ) +
  facet_wrap(.~comparison) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  )
ggsave_wrapper('ramus_mcc', width = full_page, height = full_page)

ramus_auroc %>%
  mutate(
    method = factor(method, names(color_scheeme))
  ) %>%
  ggplot(aes(comparison, auROC, fill = method, label = round(auROC, 2))) +
  geom_col(aes(comparison, auROC, fill = method), position = 'dodge') +
  geom_text(aes(comparison, auROC - .06, label = round(auROC, 2)), color = 'white', position = position_dodge(width = .9), size = 2) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_fill_manual('Method',
                    values = color_scheeme,
                    guide = guide_legend(
                      nrow = 1,
                      title.hjust = .5,
                      override.aes = list(linewidth = 1)
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
    panel.grid.major.y =  element_line('grey', .1),
    strip.text = element_text(size = 7)
  )
ggsave_wrapper('ramus_auroc', width = full_page, height = full_page)

ramus_roc %>%
  filter(between(alpha, 0, 1)) %>%
  arrange(alpha) %>%
  filter(between(alpha, 0, 1)) %>%
  ggplot(aes(alpha, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
  scale_y_continuous(breaks = seq(0,1,.25), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.25), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 3,
                       title.hjust = .5,
                       override.aes = list(linewidth = 1)
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
ggsave_wrapper('ramus_tpr', width = full_page, height = full_page)

ramus_roc %>%
  arrange(alpha) %>%
  filter(between(alpha, 0, 1)) %>%
  ggplot(aes(alpha, FPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
  scale_y_continuous(breaks = seq(0,1,.25), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.25), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 3,
                       title.hjust = .5,
                       override.aes = list(linewidth = 1)
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
    y = 'False Positive Rate'
  )
ggsave_wrapper('ramus_fpr', width = full_page, height = full_page)

ramus_roc %>%
  group_by(method, comparison) %>%
  filter(min((alpha - .05)^2) == (alpha - .05)^2) %>%
  pivot_longer(c(TP, FP)) %>%
  ggplot(aes(method, value, fill = name, label = value)) +
  geom_col(position = position_dodge(1)) +
  geom_text(position = position_dodge(1), vjust = .5) +
  theme_classic() +
  scale_y_continuous(expand = c(0,20)) +
  labs(
    y = 'Count', x = 'Method', fill = 'Category'
  ) +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(size = 7)
  ) +
  facet_wrap(comparison~.)
ggsave_wrapper('bar_alpha_05_ramus', width = full_page, height = full_page)

# Human
hr <- human_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  mutate(
    shape = case_when(
      (alpha - .05)^2 == min((alpha - .05)^2) ~ '0.05',
      (alpha - .01)^2 == min((alpha - .01)^2) ~ '0.01',
      (alpha - .001)^2 == min((alpha - .001)^2) ~ '0.001',
      T ~ 'F'
    )
  ) %>%
  filter(between(alpha, 0, 1)) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
  geom_point(
    data = \(x) subset(x, shape != 'F'),
    aes(FPR, TPR, color = method, shape = shape)
  ) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_scheeme,
                     guide = guide_legend(
                       ncol = 3,
                       title.hjust = .5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  geom_text(data = human_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  facet_wrap(factor(
    comparison,
    levels = paste0('1:', c(25, 12, 25), ' vs ', '1:', c(6, 6, 12))
  )~., nrow = 3) +
  theme(
    legend.position = 'none',
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.text = element_text(size = 7),
    legend.margin=margin(0,0,.01,0),
    axis.title = element_blank()
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  )

cowplot::plot_grid(ur, hr, labels = 'AUTO', ncol = 2) + #, hjust = c(-1.9, -1.8), rel_widths = c(1, 1)) +
  theme(
    plot.margin = margin(0, 0, 15, 15)
  ) +
  cowplot::draw_label("False Positive Rate", x = 0.5, y = 0,   angle = 0,  vjust = 1, hjust = 0.45) +
  cowplot::draw_label("True Positive Rate",  x = 0,   y = 0.5, angle = 90, vjust = -.35, hjust = .52)
ggsave_wrapper('ups_human_roc', width = full_page, height = full_page)

human_roc %>%
  filter(between(alpha, 0, 1)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dashed') +
  geom_path(linewidth = 1/4) +
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
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
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
ggsave_wrapper('human_mcc', width = full_page, height = half_page)

human_roc %>%
  filter(between(alpha, 0, 1)) %>%
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
  geom_path(linewidth = 1/4) +
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
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
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
ggsave_wrapper('human_decomp', width = full_page, height = half_page)

human_roc %>%
  group_by(method, comparison) %>%
  filter(min((alpha - .05)^2) == (alpha - .05)^2) %>%
  pivot_longer(c(TP, FP)) %>%
  ggplot(aes(method, value, fill = name, label = value)) +
  geom_col(position = position_dodge(1)) +
  geom_text(position = position_dodge(1), vjust = .5) +
  theme_classic() +
  scale_y_continuous(expand = c(0,20)) +
  labs(
    y = 'Count', x = 'Method', fill = 'Category'
  ) +
  theme(
    legend.position = c(0.1,.9),
    axis.text.x = element_text(size = 7)
  ) +
  facet_grid(comparison~.)
ggsave_wrapper('bar_alpha_05_human', width = half_page, height = full_page)
