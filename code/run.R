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
  message("scenario 1: no unmeasured confounding")
  p <- progressr::progressor(steps = NSIMS)
  
  sim1 <- dgp1[, .(
    value = unlist(runsim(.SD, p,  exp(-1))),
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
  message("scenario 2: equi-confounding")
  p <- progressor(steps = NSIMS)
  
  sim2 <- dgp2[, .(
    value = unlist(runsim(.SD, p,  exp(-1))),
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
  message("scenario 3: exclusion restriction violated")
  p <- progressor(steps = NSIMS)
  
  sim3 <- dgp3[, .(
    value = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})


# scenario 4: equi-confounding violated -----------------------------------

dgp4 <- 
  datagen(
    N = N,
    sims = NSIMS,
    fV = function(X, U, Z) plogis(-0.9 - X + 2 * U),
    fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + 0.25 * U),
    fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
    fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) + 0.25 * X + 0.25 * U)
  )

with_progress({
  message("scenario 4: equi-confounding violated")
  p <- progressor(steps = NSIMS)
  
  sim4 <- dgp4[, .(
    value = unlist(runsim(.SD, p,  exp(-1))),
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
  message("scenario 5: equi-selection violated")
  p <- progressor(steps = NSIMS)
  
  sim5 <- dgp5[, .(
    value = unlist(runsim(.SD, p,  exp(-1))),
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
  message("scenario 6: equal effect of V on T")
  p <- progressor(steps = NSIMS)
  
  sim6 <- dgp6[, .(
    value = unlist(runsim(.SD, p,  exp(-1))),
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
  message("scenario 7: unequal effect of V on T")
  p <- progressor(steps = NSIMS)
  
  sim7 <- dgp7[, .(
    value = unlist(runsim(.SD, p,  exp(-1))),
    type = types
  ), by = sim]
})

# scenario 8: effect heterogeneity ----------------------------------------

set.seed(2873645)

NSIMS <- 2000
N <- 15000


# part 1: all models correct ----------------------------------------------

dgp8_correct <-
  datagen2(
    N = N,
    sims = NSIMS,
    fV = function(X, W, U) exp(-1.5 - X + 1.5 * U),
    fI1 = function(V, X, W, U) exp(-2.1 + U),
    fI2 = function(V, X, W, U) exp(-2.4 - 0.25 * V - 1.5 * V * X + 0.6 * X + U),
    fT = function(I, V, X, W, U) exp(-1.1 + 0.25 * I(I == 2) + 0.25 * U),
  )

with_progress({
  message("scenario 8, part 1: both correctly specified")
  p <- progressor(steps = NSIMS)
  
  sim8_correct <- dgp8_correct[, .(
    value = unlist(runsim(.SD, p,  0.4229406, family = "poisson")),
    type = types
  ), by = sim]
})

# log scale
with_progress({
  message("scenario 8, part 1: both correctly specified (log)")
  p <- progressor(steps = NSIMS)
  
  sim8_log_correct <- dgp8_correct[, .(
    value = unlist(runsim(.SD, p,  0.4229406, family = "poisson", log_bias = TRUE)),
    type = types
  ), by = sim]
})


# part 2: propensity score model misspecified -----------------------------

dgp8_miss_V <-
  datagen2(
    N = N,
    sims = NSIMS,
    fV = function(X, W, U) exp(-1.75 - I(X > 0.5) + 1.5 * U),
    fI1 = function(V, X, W, U) exp(-2.1 + U),
    fI2 = function(V, X, W, U) exp(-2.4 - 0.25 * V - 1.5 * V * X + 0.6 * X + U),
    fT = function(I, V, X, W, U) exp(-1.1 + 0.25 * I(I == 2) + 0.25 * U),
  )

with_progress({
  message("scenario 8, part 2: propensity score model misspecified")
  p <- progressor(steps = NSIMS)
  
  sim8_miss_V <- dgp8_miss_V[, .(
    value = unlist(runsim(.SD, p,  0.4436199, family = "poisson", log_bias = TRUE)),
    type = types
  ), by = sim]
})


# part 3: outcome model misspecified --------------------------------------

dgp8_miss_I2 <-
  datagen2(
    N = N,
    sims = NSIMS,
    fV = function(X, W, U) exp(-1.5 - X + 1.5 * U),
    fI1 = function(V, X, W, U) exp(-2.1 + U),
    fI2 = function(V, X, W, U) exp(-2.4 - 0.25 * V - 1.5 * V * X + 2.4 * (X - 0.5)^2 + U),
    fT = function(I, V, X, W, U) exp(-1.1 + 0.25 * I(I == 2) + 0.25 * U),
  )

with_progress({
  message("scenario 8, part 3: outcome model misspecified")
  p <- progressor(steps = NSIMS)
  
  sim8_miss_I2 <- dgp8_miss_I2[, .(
    value = unlist(runsim(.SD, p,  0.4664839, family = "poisson", log_bias = TRUE)),
    type = types
  ), by = sim]
})


# part 4: both misspecified -----------------------------------------------

dgp8_miss_all <-
  datagen2(
    N = N,
    sims = NSIMS,
    fV = function(X, W, U) exp(-1.5 - I(X > 0.5) + 1.5 * U),
    fI1 = function(V, X, W, U) exp(-2.1 + U),
    fI2 = function(V, X, W, U) exp(-2.4 - 0.25 * V - 1.5 * V * X + 2.4 * (X - 0.5)^2 + U),
    fT = function(I, V, X, W, U) exp(-1.1 + 0.25 * I(I == 2) + 0.25 * U),
  )

with_progress({
  message("scenario 8, part 4: both misspecified")
  p <- progressor(steps = NSIMS)
  
  sim8_miss_all <- dgp8_miss_all[, .(
    value = unlist(runsim(.SD, p,  0.4834371, family = "poisson", log_bias = TRUE)),
    type = types
  ), by = sim]
})


# scenario 9: nonexclusive ------------------------------------------------

# dgp9 <- 
#   datagen(
#     N = N,
#     sims = NSIMS,
#     fV = function(X, U) plogis(-0.9 - X + 2 * U),
#     fI1 = function(V, X, U) exp(-2.1 - 0.5 * X + U),
#     fI2 = function(V, X, U) exp(-2.4 - V - 0.625 * X + U),
#     fT = function(I, V, X, U) exp(-1.1 + 0.5 * I(I == 2) + 0.25 * X + 0.25 * U),
#     exclusive = FALSE
#   )
# 
# with_progress({
#   p <- progressor(steps = NSIMS)
#   p("scenario 8: nonexclusive", class = "sticky")
#   
#   sim9 <- dgp9[, .(
#     RR = unlist(runsim(.SD, p,  exp(-1))), #exclusive = FALSE)),
#     type = types
#   ), by = sim]
# })


sims <- list(
  "scenario 1: no unmeasured confounding" = sim1,
  "scenario 2: equi-confounding" = sim2,
  "scenario 3: direct effect of V on I = 1" = sim3,
  "scenario 4: equi-confounding violated" = sim4,
  "scenario 5: equi-selection violated" = sim5,
  "scenario 6: equal effect of V on T" = sim6,
  "scenario 7: unequal effect of V on T" = sim7,
  "scenario 8: effect heterogeneity" = sim8_correct,
  "scenario 8: both correctly specified" = sim8_log_correct,
  "scenario 8: propensity score model misspecified" = sim8_miss_V,
  "scenario 8: outcome model misspecified" = sim8_miss_I2,
  "scenario 8: both misspecified" = sim8_miss_all
)

sims <- rbindlist(sims, idcol = TRUE)

write_rds(sims, "data/sims.rds")


