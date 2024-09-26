library(dplyr)
library(tidyr)
library(stringr)
library(kableExtra)

sim_results <-
  sims[!sims$type %in% c('sample_size'), .(est = mean(RR), sd = sd(RR)),
       by = list(.id, type)]

sim_results |>
  mutate(
    estimator = type,
    estimator = str_remove(estimator, "_[a-z]+$"),
    type = str_remove(str_extract(type, "_[a-z]+$"), "_"),
    estimator = factor(
      estimator,
      levels = c(
        "cohort_reg_U",
        "cohort_reg_noU",
        "did_reg",
        "logit_reg",
        "rrv_om",
        "rrv_ipw",
        "rrv_dr"
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
  ) |>
  pivot_longer(c(est, sd)) |>
  filter(name == "est" | (name == "sd" & type == "bias")) |>
  mutate(
    type = replace(type, name == "sd", "sd"),
    type = factor(
      type,
      levels = c(
        "bias",
        "sd",
        "cover",
        "len"
      ),
      labels = c(
        "Bias",
        "SE",
        "Coverage",
        "CIL"
      )
    )
    ) |>
  select(-name) |>
  pivot_wider(names_from = estimator, values_from = value) |>
  select(-.id) |>
  filter(type != "CIL") |>
  rename(Statistic = type) |>
  kable(
    format = "latex",
    digits = 3,
    booktabs = TRUE,
    align = c("l", rep("d", 7)),
    col.names = c(
      "Statistic",
      "{\\TableHead{TND,\\\\logit}}",
      "{\\TableHead{TND,\\\\om}}",
      "{\\TableHead{TND,\\\\ipw}}",
      "{\\TableHead{TND,\\\\dr}}",
      "{DiD}",
      "{\\TableHeadd{cohort,\\\\U}}",
      "{\\TableHeadd{cohort,\\\\no U}}"
    ),
    linesep = ""
  ) |>
  kable_styling() |>
  group_rows(
    group_label = "scenario 1: no unmeasured confounding",
    start_row = 1,
    end_row = 3,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "scenario 2: equi-confounding",
    start_row = 4,
    end_row = 6,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "scenario 3: direct effect of V on I = 1",
    start_row = 7,
    end_row = 9,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "scenario 4: equi-confounding violated",
    start_row = 10,
    end_row = 12,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "scenario 5: equi-selection violated",
    start_row = 13,
    end_row = 15,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "scenario 6: equal effect of V on T",
    start_row = 16,
    end_row = 18,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |> 
  group_rows(
    group_label = "scenario 7: unequal effect of V on T",
    start_row = 19,
    end_row = 21,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "scenario 8: I = 1 and I = 2 not mutually exclusive",
    start_row = 22,
    end_row = 24,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "scenario 9: effect heterogeneity",
    start_row = 25,
    end_row = 27,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |> 
  save_kable("results/sims.tex")


