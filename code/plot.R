library(patchwork)
library(tidyverse)

grayscale <- FALSE 

sims_plot <- 
  sims[!sims$type %in% c('sample_size') & !stringr::str_detect(sims$type, "(cover)|(len)"),]

sims_plot$type <- factor(
  sims_plot$type,
  levels = c(
    "cohort_reg_U_bias",
    "cohort_reg_noU_bias",
    "did_reg_bias",
    "logit_reg_bias",
    "rrv_om_bias",
    "rrv_ipw_bias",
    "rrv_dr_bias"
  ),
  labels = c(
    "cohort,\n U measured",
    "cohort,\n U unmeasured",
    "cohort,\n DiD",
    "TND,\n logit",
    "TND,\n om",
    "TND,\n ipw",
    "TND,\n dr"
  )
)
  
true_ve <- data.frame(
  .id = c(
    "scenario 1: no unmeasured confounding",
    "scenario 2: equi-confounding",
    "scenario 3: direct effect of V on I = 1",
    "scenario 4: equi-confounding violated",
    "scenario 5: equi-selection violated",
    "scenario 6: equal effect of V on T",
    "scenario 7: unequal effect of V on T"
  ),
  value = c(rep(1 - exp(-1), 5), 1 - exp(-1.25), rep(1 - exp(-1), 1)),
  alt_value = c(NA, NA, 1 - exp(-0.9), NA, NA, 1 - exp(-1), NA)
)



# plot 1 ------------------------------------------------------------------

p <- ggplot(
  subset(sims_plot, 
         !stringr::str_detect(.id, "scenario 8") &
           type %in% c("cohort,\n U measured",
                       "cohort,\n U unmeasured",
                       "cohort,\n DiD",
                       "TND,\n logit")),
  aes(
    x = type,
    y = 1 - (value + exp(-1)), # convert back to raw scale and make it VE
    fill = type,
    color = type
  )
) +
  facet_wrap(~gsub(":", ":\n", .id)) +
  geom_hline(aes(yintercept = value),
             data = true_ve,
             linetype = 'dashed') +
  geom_hline(aes(yintercept = alt_value), 
             data = true_ve,
             linetype = 'dashed', 
             color = "grey") +
  geom_boxplot(
    aes(color = type),
    width = 0.4,
    # fill = "white",
    alpha = 0.4,
    outlier.shape = NA
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  #scale_color_manual(values = brewer_pal(palette = "Blues")(5)[1:4]) +
  #scale_fill_manual(values = c("#e6194B", "#3cb44b", "#4363d8")) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  labs(
    x = NULL,
    y = "Estimated VE"
  ) +
  ylim(c(0.35, 0.85)) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
    )

ggsave(
  filename = "results/sims1.pdf",
  plot = p,
  device = "pdf",
  width = 6.5,
  height = 5
)


# plot 2 ------------------------------------------------------------------

true_ve8 <- data.frame(
  .id = c(
    "scenario 8: effect heterogeneity",
    "scenario 8: propensity score model misspecified",
    "scenario 8: outcome model misspecified",
    "scenario 8: both misspecified"
  ),
  truth = c(1 - 0.4229406, 1 - 0.37351, 1 - 0.4664839, 1 - 0.4019365)
)

levels <- c(
  "scenario 8: effect heterogeneity",
  "scenario 8: propensity score model misspecified",
  "scenario 8: outcome model misspecified",
  "scenario 8: both misspecified"
)

sims_plot <- 
  left_join(sims_plot, true_ve8, by = ".id") |>
  mutate(
    id = factor(.id, 
                levels = levels, 
                labels = gsub(":", ":\n", levels))
  )

p1 <- ggplot(filter(sims_plot, .id %in% c("scenario 8: effect heterogeneity")),
             aes(
               x = type,
               y = value + truth,
               fill = type,
               color = type
             )
) +
  facet_wrap(~id, ncol = 3) +
  geom_hline(aes(yintercept = truth),
             linetype = 'dashed') +
  geom_boxplot(
    width = 0.4,
    # fill = "white",
    alpha = 0.4,
    outlier.shape = NA
  ) +
  scale_fill_manual(values = scales::brewer_pal(palette = "Set1")(8)[c(1:7)]) +
  scale_color_manual(values = scales::brewer_pal(palette = "Set1")(8)[c(1:7)]) +
  #scale_color_manual(values = brewer_pal(palette = "Blues")(5)[1:4]) +
  #scale_fill_manual(values = c("#e6194B", "#3cb44b", "#4363d8")) +
  labs(
    x = NULL,
    y = "Estimated VE"
  ) +
  ylim(c(0.4, 0.8)) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

p2 <- ggplot(filter(sims_plot, 
                    .id %in% c("scenario 8: propensity score model misspecified",
                               "scenario 8: outcome model misspecified",
                               "scenario 8: both misspecified") &
                      type %in% c("TND,\n logit",
                                  "TND,\n om",
                                  "TND,\n ipw",
                                  "TND,\n dr")),
             aes(
               x = type,
               y = value,
               fill = type,
               color = type
             )) +
  facet_wrap(~id, ncol = 3) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed') +
  geom_boxplot(
    width = 0.4,
    # fill = "white",
    alpha = 0.4,
    outlier.shape = NA
  ) +
  scale_fill_manual(values = scales::brewer_pal(palette = "Set1")(8)[c(4:7)]) +
  scale_color_manual(values = scales::brewer_pal(palette = "Set1")(8)[c(4:7)]) +
  #scale_color_manual(values = brewer_pal(palette = "Blues")(5)[1:4]) +
  #scale_fill_manual(values = c("#e6194B", "#3cb44b", "#4363d8")) +
  labs(
    x = NULL,
    y = "Bias"
  ) +
  ylim(c(-0.3, 0.2)) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = "results/sims2.pdf",
  plot = p1 / p2,
  device = "pdf",
  width = 6.5,
  height = 5
)
