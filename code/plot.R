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
    "cohort, U measured",
    "cohort, U unmeasured",
    "DiD",
    "TND, logit",
    "TND, om",
    "TND, ipw",
    "TND, dr"
  )
)
  
true_rr <- data.frame(
  .id = c(
    "scenario 1: no unmeasured confounding",
    "scenario 2: equi-confounding",
    "scenario 3: direct effect of V on I = 1",
    "scenario 4: equi-confounding violated",
    "scenario 5: equi-selection violated",
    "scenario 6: equal effect of V on T",
    "scenario 7: unequal effect of V on T",
    "scenario 8: I = 1 and I = 2 not mutually exclusive",
    "scenario 9: effect heterogeneity"
  ),
  value = c(rep(1 - exp(-1), 5), 1 - exp(-1.25), rep(1 - exp(-1), 2), 0.5539383),
  alt_value = c(NA, NA, 1 - exp(-0.9), NA, NA, 1 - exp(-1), NA, NA, NA)
)



# plot 1 ------------------------------------------------------------------

p <- ggplot(
  subset(sims_plot, 
         !.id %in% c("scenario 8: I = 1 and I = 2 not mutually exclusive",
                     "scenario 9: effect heterogeneity") &
           type %in% c("cohort, U measured",
                       "cohort, U unmeasured",
                       "DiD",
                       "TND, logit")),
  aes(
    x = gsub(", logit", "", type, fixed = TRUE),
    y = 1 - (RR + exp(-1)),
    fill = gsub(", logit", "", type, fixed = TRUE),
    color = gsub(", logit", "", type, fixed = TRUE)
  )
) +
  facet_wrap(~gsub(":", ":\n", .id)) +
  geom_hline(aes(yintercept = value),
             data = subset(
               true_rr,
               !.id %in% c(
                 "scenario 8: I = 1 and I = 2 not mutually exclusive",
                 "scenario 9: effect heterogeneity"
               )
             ),
             linetype = 'dashed') +
  geom_hline(aes(yintercept = alt_value), 
             data = subset(
               true_rr,
               !.id %in% c(
                 "scenario 8: I = 1 and I = 2 not mutually exclusive",
                 "scenario 9: effect heterogeneity"
               )
             ),
             linetype = 'dashed', 
             color = "grey") +
  # geom_hline(aes(yintercept = 1 - est), data = de_est, linetype = 'dashed', color = "grey") +
  #geom_violin() +
  geom_boxplot(
    aes(color = gsub(", logit", "", type, fixed = TRUE)),
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


p <- ggplot(
  subset(sims_plot, 
         .id %in% c("scenario 9: effect heterogeneity")),
  aes(
    x = type,
    y = 1 - (RR + exp(-1)),
    fill = type,
    color = type
  )
) +
  facet_wrap(~gsub(":", ":\n", .id)) +
  geom_hline(aes(yintercept = value),
             data = subset(
               true_rr,
               .id %in% c(
                 "scenario 9: effect heterogeneity"
               )
             ),
             linetype = 'dashed') +
  geom_hline(aes(yintercept = alt_value), 
             data = subset(
               true_rr,
               .id %in% c(
                 "scenario 9: effect heterogeneity"
               )
             ),
             linetype = 'dashed', 
             color = "grey") +
  # geom_hline(aes(yintercept = 1 - est), data = de_est, linetype = 'dashed', color = "grey") +
  #geom_violin() +
  geom_boxplot(
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
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = "results/sims2.pdf",
  plot = p,
  device = "pdf",
  width = 6.5,
  height = 5
)
