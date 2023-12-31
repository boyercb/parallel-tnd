library(data.table)
library(ggplot2)
library(progressr)

source("dgp.R")

NSIMS <- 100
N <- 100000

types <- c("TND, om (marg)",
           "TND, om",
           "cohort, U unmeasured",
           "cohort, U measured", 
           "cohort, U measured (OR)")


# scenario 1: no unmeasured confounding -----------------------------------

dgp1 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X),
    fI1 = function(V, X, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(V, X, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT = function(I, V, X, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 1: no unmeasured confounding", class = "sticky")
  
  sim1 <- dgp1[, .(
    VE = unlist(runsim(.SD, p)),
    type = types
  ), by = sim]
})


# scenario 2: equi-confounding --------------------------------------------

dgp2 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(V, X, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(V, X, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT = function(I, V, X, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 2: equi-confounding", class = "sticky")
  
  sim2 <- dgp2[, .(
    VE = unlist(runsim(.SD, p)),
    type = types
  ), by = sim]
})


# scenario 3: exclusion restriction violated ------------------------------

dgp3 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(V, X, U) plogis(-2.75 - 0.2 * V + 0.135 * X + 0.4 * U),
    fI2 = function(V, X, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT = function(I, V, X, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 3: direct effect of V on I1", class = "sticky")
  
  sim3 <- dgp3[, .(
    VE = unlist(runsim(.SD, p)),
    type = types
  ), by = sim]
})


# scenario 4: equi-confounding violated -----------------------------------

dgp4 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(V, X, U) plogis(-2.75 + 0.135 * X + 0.1 * U),
    fI2 = function(V, X, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT = function(I, V, X, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 4: equi-confounding violated", class = "sticky")
  
  sim4 <- dgp4[, .(
    VE = unlist(runsim(.SD, p)),
    type = types
  ), by = sim]
})


# scenario 5: nonexclusive ------------------------------------------------

dgp5 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(V, X, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(V, X, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT = function(I, V, X, U) plogis(0.5 + 0.5 * X + 0.5 * U),
    exclusive = FALSE
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 5: nonexclusive", class = "sticky")
  
  sim5 <- dgp5[, .(
    VE = unlist(runsim(.SD, p)),
    type = types
  ), by = sim]
})


# scenario 6: equi-selection violated -------------------------------------

dgp6 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(V, X, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(V, X, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT = function(I, V, X, U) plogis(0.5 + 0.5 * X + 0.5 * U - I(I == 2) * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 6: equi-selection violated", class = "sticky")
  
  sim6 <- dgp6[, .(
    VE = unlist(runsim(.SD, p)),
    type = types
  ), by = sim]
})


# scenario 7: direct effect of V on T -------------------------------------

dgp7 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(V, X, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(V, X, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT = function(I, V, X, U) plogis(0.5 - V + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 7: direct effect of V on T", class = "sticky")
  
  sim7 <- dgp7[, .(
    VE = unlist(runsim(.SD, p, effect_modification = TRUE)),
    type = types
  ), by = sim]
})


# scenario 8: effect heterogeneity ----------------------------------------

dgp8 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(V, X, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(V, X, U) plogis(-2.75  - 0.25 * V - 0.25 * V * X + 0.125 * X + 0.4 * U),
    fT = function(I, V, X, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 8: effect heterogeneity", class = "sticky")
  
  sim8 <- dgp8[, .(
    VE = unlist(runsim(.SD, p, effect_modification = TRUE)),
    type = types
  ), by = sim]
})


# plot --------------------------------------------------------------------

sims <- list(
  "scenario 1: no unmeasured confounding" = sim1, 
  "scenario 2: equi-confounding" = sim2, 
  "scenario 3: direct effect of V on I1" = sim3, 
  "scenario 4: equi-confounding violated" = sim4,
  "scenario 5: I1 and I2 not mutually exclusive" = sim5,
  "scenario 6: equi-selection violated" = sim6,
  "scenario 7: direct effect of V on T" = sim7,
  "scenario 8: effect heterogeneity" = sim8 
)

sims <- rbindlist(sims, idcol = TRUE)

sim_results <-
  sims[, .(est = mean(VE)),
  by = list(.id, type)]

truth_est <- sim_results[type == 'cohort, U measured']
or_est <- sim_results[type == 'cohort, U measured (OR)']
or_est <- or_est[5, ]
sims <- sims[sims$type != 'cohort, U measured (OR)', ]

ggplot(
  sims,
  aes(x = type, y = 1 - VE, fill = type)
) +
  facet_wrap(~.id) +
  geom_hline(aes(yintercept = 1 - est), data = truth_est, linetype = 'dashed') +
  geom_hline(aes(yintercept = 1 - est), data = or_est, linetype = 'dashed', color = "grey") +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  # scale_color_manual(values = brewer_pal(palette = "Blues")(5)[1:4]) +
  # scale_fill_manual(values = brewer_pal(palette = "Blues")(5)[1:4]) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  labs(
    x = NULL,
    y = "Estimated VE"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

