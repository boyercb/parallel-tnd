# data generation ---------------------------------------------------------

# function to generate data from test-negative design under specified process
gendata <- function(
    N = 1000,                                  # population size
    sims = 1,                                  # number of simulations
    fV = function(X, U) plogis(X),             # function for generating V
    fI1 = function(X, V, U) plogis(X + U),     # function for generating I1
    fI2 = function(X, V, U) plogis(X - V + U), # function for generating I2
    fT1 = function(X, V, U) plogis(X + U),     # function for generating T among I1 == 1
    fT2 = function(X, V, U) plogis(X + U),     # function for generating T among I2 == 1
    exclusive = TRUE,                          # are I1 and I2 mutually exclusive?
    sens = 1,                                  # test sensitivity
    spec = 1                                   # test specificity
) {
  n <- N * sims
  
  U <- runif(n, 0, 5) 
  X <- runif(n, 0, 5)
  V <- rbinom(n, 1, fV(X, U))
  
  I1 <- rbinom(n, 1, fI1(X, V, U))
  
  if (exclusive) {
    I2 <- (1 - I1) * rbinom(n, 1, fI2(X, V, U) / (1 - fI1(X, V, U)))
    
  } else {
    I2 <- rbinom(n, 1, fI2(X, V, U))
    
  }
  
  T <- I(I1 == 1) * rbinom(n, 1, fT1(X, V, U)) +
    I(I2 == 1) * rbinom(n, 1, fT2(X, V, U))
  
  S <- T
  
  Istar <- I(I2 == 1) * S * rbinom(n, 1, spec) +
    I(I1 == 1) * S * (1 - rbinom(n, 1, sens))
  
  dt <- data.table(
    sim = rep(1:sims, each = N),
    U = U,
    X = X,
    V = V,
    I = I1 + 2 * I2,
    T = T,
    S = S,
    Istar = Istar
  )
  setkey(dt, sim)
  return(dt)
}


# simulation function -----------------------------------------------------

runsim <- function(dt, p, effect_modification = FALSE) {
  p()
  dt <- as.data.frame(dt)
  if (effect_modification) {
    f <- ~ V + X + V:X
  } else {
    f <- ~ V + X 
  }
  
  d <- dt[dt$S == 1, ]
  d1 <- d[d$V == 1, ]
  d0 <- transform(d[d$V == 1, ], V = 0)
  dt1 <- dt[dt$V == 1, ]
  dt0 <- transform(dt[dt$V == 1, ], V = 0)
  # dt1 <- transform(dt[dt$V == 1], V = 1)
  # dt0 <- transform(dt[dt$V == 1], V = 0)

  # fit conditional odds ratio estimator
  fit_om <-
    glm(
      formula = update.formula(f, Istar ~ .),
      family = binomial(link = "logit"),
      data = d
    )
  
  p1 <- predict(fit_om, newdata = d1, type = "response")
  p0 <- predict(fit_om, newdata = d0, type = "response")
  o1 <- exp(predict(fit_om, newdata = d1))
  o0 <- exp(predict(fit_om, newdata = d0))
  
  # fit ipw estimator
  fit_ipw <- 
    glm(
      formula = V ~ X,
      family = binomial(link = "logit"),
      data = d[d$Istar == 0, ]
    )

  pV <- predict(fit_ipw, newdata = d, type = "response")
  oV <- exp(predict(fit_ipw, newdata = d))
  # fit_ipw2 <- 
  #   glm(
  #     formula = V ~ X + W + X:W,
  #     family = binomial(link = "logit"),
  #     data = d
  #   )
  # 
  # Vor <- 1 / exp(predict(fit_ipw2, newdata = d))
  # 
  # fit cohort
  fit_cohort <-
    glm(
      formula = update.formula(f, I(I == 2) ~ .),
      family = binomial(link = "logit"),
      data = dt
    )
  
  r1_cohort <- predict(fit_cohort, newdata = dt1, type = 'response')
  r0_cohort <- predict(fit_cohort, newdata = dt0, type = 'response')
  
  # fit truth
  fit_truth <-
    glm(
      formula = update.formula(f, I(I == 2) ~ . + U),
      family = binomial(link = "logit"),
      data = dt
    )

  r1_truth <- predict(fit_truth, newdata = dt1, type = 'response')
  r0_truth <- predict(fit_truth, newdata = dt0, type = 'response')
  # print(head(p1))
  # print(head(p0))
  # a <- o1 / o0
  # b <- p0
  # print(summary(a))

  c(
    mean(d1$Istar) / mean(o0 / o1 * p1),

    #mean((a * b) / (a * b + 1 - b)),
    # mean(d1$Istar) / mean((o0 / o1 * p1) / (o0 / o1 * p1 + 1 - p1)),
    # mean(d1$Istar) / mean((1 - p1) / (1 - p0) * p0),
    mean(d$V / pV * d$Istar) /
      mean((1 - d$V) / (1 - pV) * d$Istar),
    # exp(coef(fit_om)[2]),
    # mean(d1$Istar) / (mean((1 - d$V) * oV * d$Istar) / mean((1 - d$V) * oV)),
    mean(r1_cohort) / mean(r0_cohort),
    mean(r1_truth) / mean(r0_truth)
  )
  
}
