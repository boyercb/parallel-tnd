library(data.table)
library(progressr)

source("dgp.R")

NSIMS <- 100
N <- 100000


# scenario 1: no unmeasured confounding -----------------------------------

dgp1 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X),
    fI1 = function(X, V, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(X, V, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT1 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U),
    fT2 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 1: no unmeasured confounding", class = "sticky")
  
  sim1 <- dgp1[, .(
    VE = unlist(runsim(.SD, p)),
    type = c("TND, om",
             "TND, ipw",
             "cohort, U unmeasured",
             "cohort, U measured")
  ), by = sim]
})


# scenario 2: equi-confounding --------------------------------------------

dgp2 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(X, V, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(X, V, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT1 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U),
    fT2 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 2: equi-confounding", class = "sticky")
  
  sim2 <- dgp2[, .(
    VE = unlist(runsim(.SD, p)),
    type = c("TND, om",
             "TND, ipw",
             "cohort, U unmeasured",
             "cohort, U measured")
  ), by = sim]
})


# scenario 3: exclusion restriction violated ------------------------------

dgp3 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(X, V, U) plogis(-2.75 - 0.2 * V + 0.135 * X + 0.4 * U),
    fI2 = function(X, V, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT1 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U),
    fT2 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 3: exclusion restriction violated", class = "sticky")
  
  sim3 <- dgp3[, .(
    VE = unlist(runsim(.SD, p)),
    type = c("TND, om",
             "TND, ipw",
             "cohort, U unmeasured",
             "cohort, U measured")
  ), by = sim]
})


# scenario 4: equi-confounding violated -----------------------------------

dgp4 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(X, V, U) plogis(-2.75 + 0.135 * X + 0.1 * U),
    fI2 = function(X, V, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT1 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U),
    fT2 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 4: equi-confounding violated", class = "sticky")
  
  sim4 <- dgp4[, .(
    VE = unlist(runsim(.SD, p)),
    type = c("TND, om",
             "TND, ipw",
             "cohort, U unmeasured",
             "cohort, U measured")
  ), by = sim]
})


# scenario 5: nonexclusive ------------------------------------------------

dgp5 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(X, V, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(X, V, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT1 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U),
    fT2 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U),
    exclusive = FALSE
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 5: nonexclusive", class = "sticky")
  
  sim5 <- dgp5[, .(
    VE = unlist(runsim(.SD, p)),
    type = c("TND, om",
             "TND, ipw",
             "cohort, U unmeasured",
             "cohort, U measured")
  ), by = sim]
})


# scenario 6: equi-selection violated -------------------------------------

dgp6 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(X, V, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(X, V, U) plogis(-2.75 - V + 0.125 * X + 0.4 * U),
    fT1 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.1 * U),
    fT2 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 6: equi-selection violated", class = "sticky")
  
  sim6 <- dgp6[, .(
    VE = unlist(runsim(.SD, p)),
    type = c("TND, om",
             "TND, ipw",
             "cohort, U unmeasured",
             "cohort, U measured")
  ), by = sim]
})



# scenario 7: effect heterogeneity ----------------------------------------


dgp7 <- 
  gendata(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-1 + 0.125 * X + 0.3 * U),
    fI1 = function(X, V, U) plogis(-2.75 + 0.135 * X + 0.4 * U),
    fI2 = function(X, V, U) plogis(-2.75  - 0.25 * V - 0.25 * V * X + 0.125 * X + 0.4 * U),
    fT1 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U),
    fT2 = function(X, V, U) plogis(0.5 + 0.5 * X + 0.5 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 7: effect heterogeneity", class = "sticky")
  
  sim7 <- dgp7[, .(
    VE = unlist(runsim(.SD, p, effect_modification = TRUE)),
    type = c("TND, om",
             "TND, ipw",
             "cohort, U unmeasured",
             "cohort, U measured")
  ), by = sim]
})


# plot --------------------------------------------------------------------

sims <- list(
  "scenario 1: no unmeasured confounding" = sim1, 
  "scenario 2: equi-confounding" = sim2, 
  "scenario 3: exclusion restriction violated" = sim3, 
  "scenario 4: equi-confounding violated" = sim4,
  "scenario 5: I1 and I2 not mutually exclusive" = sim5,
  "scenario 6: equi-selection violated" = sim6,
  "scenario 7: effect heterogeneity" = sim7 
)

sims <- rbindlist(sims, idcol = TRUE)

sim_results <-
  sims[, .(est = mean(VE)),
  by = list(.id, type)]
sim_results
truth_est <- sim_results[type == 'cohort, U measured']

ggplot(
  sims,
  aes(x = type, y = 1 - VE, fill = type)
) +
  facet_wrap(~.id) +
  geom_hline(aes(yintercept = 1 - est), data = truth_est, linetype = 'dashed') +
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

