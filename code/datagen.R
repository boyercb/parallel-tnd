# function to generate data from test-negative design 
datagen <- function(
    N = 1000,                                  # population size
    sims = 1,                                  # number of simulations
    fV = function(X, U) plogis(X),             # function for generating V
    fI1 = function(V, X, U) plogis(X + U),     # function for generating I1
    fI2 = function(V, X, U) plogis(X - V + U), # function for generating I2
    fT  = function(I, V, X, U) plogis(X + U),  # function for generating T 
    exclusive = TRUE,                          # are I1 and I2 mutually exclusive?
    sens = 1,                                  # test sensitivity
    spec = 1                                   # test specificity
) {
  n <- N * sims
  
  U <- runif(n, 0, 1) 
  X <- runif(n, 0, 1)
  V <- rbinom(n, 1, fV(X, U))
  I1 <- rbinom(n, 1, fI1(V, X, U))
  
  if (exclusive) {
    I2 <- rbinom(n, 1, (1 - I1) * fI2(V, X, U) / (1 - fI1(V, X, U)))
    I <- I1 + 2 * I2 - I1 * I2
    
  } else {
    I2 <- rbinom(n, 1, fI2(V, X, U))
    I <- I1 + 2 * I2 - I1 * I2
  }
  
  T <- rbinom(n, 1, I(I > 0) * fT(I, V, X, U)) 
  
  S <- T
  
  Istar <- I(I2 == 1) * S * rbinom(n, 1, spec) +
    I(I1 == 1) * S * (1 - rbinom(n, 1, sens))
  
  dt <- data.table(
    sim = rep(1:sims, each = N),
    U = U,
    X = X,
    V = V,
    I1 = I1,
    I2 = I2,
    I = I,
    T = T,
    S = S,
    Istar = Istar
  )
  setkey(dt, sim)
  return(dt)
}

expit <- function(x) {exp(x) / (1 + exp(x))}
logit <- function(x) {log(x / (1 + x))}
