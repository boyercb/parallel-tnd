library(data.table)
library(ggplot2)
library(progressr)
library(survival)
library(lmtest)
library(sandwich)
library(readr)

source("code/datagen.R")
source("code/estimators.R")
source("code/sim.R")

types <- c(
  "logit_reg_bias",
  "logit_reg_cover",
  "logit_reg_len",
  "rrv_om_bias",
  "rrv_om_cover",
  "rrv_om_len",
  "rrv_ipw_bias",
  "rrv_ipw_cover",
  "rrv_ipw_len",
  "rrv_dr_bias",
  "rrv_dr_cover",
  "rrv_dr_len",
  "did_reg_bias",
  "did_reg_cover",
  "did_reg_len",
  "cohort_reg_U_bias",
  "cohort_reg_U_cover",
  "cohort_reg_U_len",
  "cohort_reg_noU_bias",
  "cohort_reg_noU_cover",
  "cohort_reg_noU_len",
  "sample_size"
)

NSIMS <- 1000
N <- 15000

set.seed(9453602)

# scenario 1: no unmeasured confounding -----------------------------------

dgp1 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) + 0.25 * X + 0.25 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 1: no unmeasured confounding", class = "sticky")
  
  sim1 <- dgp1[, .(
    RR = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})


# scenario 2: equi-confounding --------------------------------------------

dgp2 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) + 0.25 * X + 0.25 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 2: equi-confounding", class = "sticky")
  
  sim2 <- dgp2[, .(
    RR = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})


# scenario 3: exclusion restriction violated ------------------------------

dgp3 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.1 * V - 0.5 * X + U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) + 0.25 * X + 0.25 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 3: exclusion restriction violated", class = "sticky")
  
  sim3 <- dgp3[, .(
    RR = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})


# scenario 4: equi-confounding violated -----------------------------------

dgp4 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + 0.25 * U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) + 0.25 * X + 0.25 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 4: equi-confounding violated", class = "sticky")
  
  sim4 <- dgp4[, .(
    RR = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})


# scenario 5: equi-selection violated -------------------------------------

dgp5 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) 
      exp(-1.1 + 0.5 * I(I == 2) + 0.25 * X + (0.25 - 2 * I(I == 2)) * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 5: equi-selection violated", class = "sticky")
  
  sim5 <- dgp5[, .(
    RR = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})


# scenario 6: equal effect of V on T --------------------------------------

dgp6 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) - 0.25 * V + 0.25 * X + 0.25 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 6: equal effect of V on T", class = "sticky")
  
  sim6 <- dgp6[, .(
    RR = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})


# scenario 6: equal effect of V on T --------------------------------------

dgp6 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) - 0.25 * V + 0.25 * X + 0.25 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 6: equal effect of V on T", class = "sticky")
  
  sim6 <- dgp6[, .(
    RR = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})

# scenario 7 unequal effect of V on T ------------------------------------

dgp7 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) 
      exp(-1.1 + 0.5 * I(I == 2) - 0.25 * V * I(I == 1) + 0.25 * X + 0.25 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 7: unequal effect of V on T", class = "sticky")
  
  sim7 <- dgp7[, .(
    RR = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})


# scenario 8: nonexclusive ------------------------------------------------

dgp8 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) + 0.25 * X + 0.25 * U),
    exclusive = FALSE
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 8: nonexclusive", class = "sticky")
  
  sim8 <- dgp8[, .(
    RR = unlist(runsim(.SD, p,  exp(-1))), #exclusive = FALSE)),
    type = types
  ), by = sim]
})


# scenario 9: effect heterogeneity ----------------------------------------

dgp9 <-
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + U),
    fI2 = function(V, X, U) exp(-2.4 - 0.25 * V  - 1.5 * V * X - 0.625 * X + U),
    fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) + 0.25 * X + 0.25 * U)
  )

with_progress({
  p <- progressor(steps = NSIMS)
  p("scenario 9: effect heterogeneity", class = "sticky")

  sim9 <- dgp9[, .(
    RR = unlist(runsim(.SD, p,  0.5539383)),
    type = types
  ), by = sim]
})


sims <- list(
  "scenario 1: no unmeasured confounding" = sim1,
  "scenario 2: equi-confounding" = sim2,
  "scenario 3: direct effect of V on I = 1" = sim3,
  "scenario 4: equi-confounding violated" = sim4,
  "scenario 5: equi-selection violated" = sim5,
  "scenario 6: equal effect of V on T" = sim6,
  "scenario 7: unequal effect of V on T" = sim7,
  "scenario 8: I = 1 and I = 2 not mutually exclusive" = sim8,
  "scenario 9: effect heterogeneity" = sim9 
)

sims <- rbindlist(sims, idcol = TRUE)

write_rds(sims, "data/sims.rds")


