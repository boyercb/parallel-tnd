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
    spec = 1,                                  # test specificity
    ncov = 1                                   # number of covariates in X
) {
  n <- N * sims
  
  U <- runif(n, 0, 1) 
  X <- runif(n * ncov, 0, 1) 
  #Z <- rnorm(n, 0.5 * U, 1)
    
  if (ncov > 1) {
    X <- matrix(X, ncol = ncov)
    colnames(X) <- paste0("X", 1:ncov)
  }
  
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
  
  Y <- S * I
  
  Istar <- I(I2 == 1) * S * rbinom(n, 1, spec) +
    I(I1 == 1) * S * (1 - rbinom(n, 1, sens))
  
  dt <- data.table(
    sim = rep(1:sims, each = N),
    U = U,
    X,
    #Z = Z,
    V = V,
    I1 = I1,
    I2 = I2,
    I = I,
    T = T,
    S = S,
    Y = Y,
    Istar = Istar
  )
  setkey(dt, sim)
  return(dt)
}


datagen2 <- function(
    N = 1000,                                     # population size
    sims = 1,                                     # number of simulations
    fV = function(X, W, U) plogis(X),             # function for generating V
    fI1 = function(V, X, W, U) plogis(X + U),     # function for generating I1
    fI2 = function(V, X, W, U) plogis(X - V + U), # function for generating I2
    fT  = function(I, V, X, W, U) plogis(X + U),  # function for generating T 
    exclusive = TRUE,                             # are I1 and I2 mutually exclusive?
    sens = 1,                                     # test sensitivity
    spec = 1                                      # test specificity
) {
  n <- N * sims
  
  U <- runif(n, 0, 1) 
  X <- runif(n, 0, 1)
  W <- runif(n, 0, 1)
  V <- rbinom(n, 1, fV(X, W, U))
  
  if (exclusive) {
    p1 <- fI1(V, X, W, U)
    p2 <- fI2(V, X, W, U)
    p0 <- 1 - p1 - p2
    mat <- cbind(p0, p1, p2)
    colnames(mat) <- NULL
    
    
    I <- as.vector(Hmisc::rMultinom(mat, 1)) - 1
    I2 <- as.numeric(I == 2)
    I1 <- as.numeric(I == 1)
    
  } else {
    I1 <- rbinom(n, 1, fI1(V, X, W, U))
    I2 <- rbinom(n, 1, fI2(V, X, W, U))
    I <- I1 + 2 * I2 - I1 * I2
  }
  
  T <- rbinom(n, 1, I(I > 0) * fT(I, V, X, W, U)) 
  
  S <- T
  
  Y <- S * I
  
  Istar <- I(I2 == 1) * S * rbinom(n, 1, spec) +
    I(I1 == 1) * S * (1 - rbinom(n, 1, sens))
  
  dt <- data.table(
    sim = rep(1:sims, each = N),
    U = U,
    X = X,
    W = W,
    V = V,
    I1 = I1,
    I2 = I2,
    I = I,
    T = T,
    S = S,
    Y = Y,
    Istar = Istar
  )
  setkey(dt, sim)
  return(dt)
}

expit <- function(x) {exp(x) / (1 + exp(x))}
logit <- function(x) {log(x / (1 + x))}
