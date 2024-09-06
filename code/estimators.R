
# TND estimators  ---------------------------------------------------------

LogitReg <- function(I2, V, X) {
  glm_model <- glm(I2 ~ V + X, family = binomial)
  ESTIMATE <- coef(glm_model)["V"]
  CI <- coefci(glm_model, vcov = vcovHC, type = "HC0")["V", ]
  RR <- exp(ESTIMATE)
  RR_CI <- exp(CI)
  return(list(ESTIMATE = ESTIMATE,
              CI = CI,
              RR = RR,
              RR_CI = RR_CI))
}

PoissonReg <- function(I2, V, X) {
  glm_model <- glm(I2 ~ V + X, family = poisson)
  ESTIMATE <- coef(glm_model)["V"]
  CI <- coefci(glm_model, vcov = vcovHC, type = "HC0")["V", ]
  RR <- exp(ESTIMATE)
  RR_CI <- exp(CI)
  return(list(ESTIMATE = ESTIMATE,
              CI = CI,
              RR = RR,
              RR_CI = RR_CI))
}

## estimate RR in the vaccinated using outcome regression
RRV_OM <- function(I2, V, X) {
  om_model <- glm(I2 ~ V * X, family = binomial)
  newdata_V1 <- data.frame(X = X, V = 1)
  newdata_V0 <- data.frame(X = X, V = 0)
  
  mu0 <- predict(om_model, newdata = newdata_V0, type = "response")
  mu1 <- predict(om_model, newdata = newdata_V1, type = "response")
  
  rrv_om <- sum(I2 * V) / sum(V * mu0 * (1 - mu1) / (1 - mu0))
  ## variables to facilitate variance calculation
  om_coef <- coef(om_model)
  params <- c(log(rrv_om), om_coef)
  
  om_matrix <- model.matrix(om_model)
  om_matrix_v1 <- cbind(1, 1, X, X)
  om_matrix_v0 <- cbind(1, 0, X, X * 0)
  
  # par <- params
  U_om <- function(par) {
    log_rrv <- par[1]
    reg_coef <- par[-1]
    
    mu1est <- expit(c(om_matrix_v1 %*% reg_coef))
    mu0est <- expit(c(om_matrix_v0 %*% reg_coef))
    
    U_rrv <- V * (I2 - exp(log_rrv) * mu0est *
                    (1 - mu1est) / (1 - mu0est))
    U_glm <- om_matrix * (I2 - expit(c(om_matrix %*% reg_coef)))
    
    return(cbind(U_rrv, U_glm))
  }
  
  GMMF <- function(mrf, par) {
    g0 <- mrf(par = par)
    g <- apply(g0, 2, mean)
    gmmf <- sum(g ^ 2)
    
    return(gmmf)
  }
  
  G <- function(bfun, par){
    G <- numDeriv::jacobian(func = G1, bfun = bfun, x = par)
    return(G)
  }
  
  G1 <- function(bfun, par) {
    G1 <- apply(bfun(par), 2, mean, na.rm=T)
    return(G1)
  }
  
  ## sandwich variance
  var_gmmf <- function(bfun, par){
    bG <- solve(G(bfun, par))
    bg <- bfun(par)
    spsz <- dim(bg)[1]
    Omega <- t(bg) %*% bg / spsz
    Sigma <- bG %*% Omega %*% t(bG) / spsz
    return(Sigma)
  }
  
  VAR_all <- var_gmmf(U_om, params)
  SE_all <- sqrt(diag(VAR_all))
  
  RR_CI <- exp(log(rrv_om) + qnorm(c(0.025, 0.975)) * SE_all[1])
  return(list(RR = rrv_om, RR_CI = RR_CI, VAR = VAR_all))
}

## estimate RR in the vaccinated using IPW
RRV_IPW <- function(I2, V, X) {
  ps_model <- glm(V ~ I2 * X, family = binomial)
  newdata_I0 <- data.frame(X = X, I2 = 0)
  
  pi0 <- predict(ps_model, newdata = newdata_I0, type = "response")
  
  rrv_ipw <- sum(I2 * V) / sum((1 - V) * I2 * pi0 / (1 - pi0))
  
  ## variables to facilitate variance calculation
  ps_coef <- coef(ps_model)
  params <- c(log(rrv_ipw), ps_coef)
  
  ps_matrix <- model.matrix(ps_model)
  ps_matrix_i0 <- cbind(1, 0, X, X * 0)
  
  # par <- params
  U_ipw <- function(par) {
    log_rrv <- par[1]
    reg_coef <- par[-1]
    
    pi0est <- expit(c(ps_matrix_i0 %*% reg_coef))
    
    U_rrv <- V * I2 - exp(log_rrv) * (1 - V) * I2 * pi0est / (1 - pi0est)
    U_glm <- ps_matrix * (V - expit(c(ps_matrix %*% reg_coef)))
    
    return(cbind(U_rrv, U_glm))
  }
  
  GMMF <- function(mrf, par) {
    g0 <- mrf(par = par)
    g <- apply(g0, 2, mean)
    gmmf <- sum(g ^ 2)
    
    return(gmmf)
  }
  
  G <- function(bfun, par){
    G <- numDeriv::jacobian(func = G1, bfun = bfun, x = par)
    return(G)
  }
  
  G1 <- function(bfun, par) {
    G1 <- apply(bfun(par), 2, mean, na.rm=T)
    return(G1)
  }
  
  ## sandwich variance
  var_gmmf <- function(bfun, par){
    bG <- solve(G(bfun, par))
    bg <- bfun(par)
    spsz <- dim(bg)[1]
    Omega <- t(bg) %*% bg / spsz
    Sigma <- bG %*% Omega %*% t(bG) / spsz
    return(Sigma)
  }
  
  VAR_all <- var_gmmf(U_ipw, params)
  SE_all <- sqrt(diag(VAR_all))
  
  RR_CI <- exp(log(rrv_ipw) + qnorm(c(0.025, 0.975)) * SE_all[1])
  return(list(RR = rrv_ipw, RR_CI = RR_CI, VAR = VAR_all))
}



## estimate RR in the vaccinated using the DR estimator
RRV_DR <- function(I2, V, X) {
  ps_model <- glm(V ~ X, family = binomial)
  om_model <- glm(I2 ~ V * X, family = binomial)
  
  newdata_V1 <- data.frame(X = X, V = 1)
  newdata_V0 <- data.frame(X = X, V = 0)
  newdata_x <- data.frame(X = X)
  
  mu0 <- predict(om_model, newdata = newdata_V0, type = "response")
  mu1 <- predict(om_model, newdata = newdata_V1, type = "response")
  pix <- predict(ps_model, type = "response")
  
  rrv_dr <- sum(I2 * V) / sum((1 - V) * (I2 - mu0) * (pix / (1 - pix)) * 
                                ((1 - mu1) / (1 - mu0) ^ 2) +
                                V * (1 - I2) * mu0 / (1 - mu0))
  
  ## variables to facilitate variance calculation
  ps_coef <- coef(ps_model); n_ps_coef <- length(ps_coef)
  om_coef <- coef(om_model); n_om_coef <- length(om_coef)
  params <- c(log(rrv_dr), ps_coef, om_coef)
  
  om_matrix <- model.matrix(om_model)
  om_matrix_v1 <- cbind(1, 1, X, X)
  om_matrix_v0 <- cbind(1, 0, X, X * 0)
  
  ps_matrix <- model.matrix(ps_model)
  
  # par <- params
  U_dr <- function(par) {
    log_rrv <- par[1]
    
    ps_coef <- par[2:(1 + n_ps_coef)]
    pixest <- expit(c(ps_matrix %*% ps_coef))
    
    om_coef <- par[(2 + n_ps_coef):(1 + n_ps_coef + n_om_coef)]
    mu1est <- expit(c(om_matrix_v1 %*% om_coef))
    mu0est <- expit(c(om_matrix_v0 %*% om_coef))
    
    
    U_rrv <- V * I2 - exp(log_rrv) * 
      ((1 - V) * (I2 - mu0est) * pixest * (1 - mu1est) / ((1 - pixest) * (1 - mu0est) ^ 2) +
         V * (1 - I2) * mu0est / (1 - mu0est))
    
    U_om <- om_matrix * (I2 - expit(c(om_matrix %*% om_coef)))
    U_ps <- ps_matrix * (V - expit(c(ps_matrix %*% ps_coef)))
    
    return(cbind(U_rrv, U_ps, U_om))
  }
  
  GMMF <- function(mrf, par) {
    g0 <- mrf(par = par)
    g <- apply(g0, 2, mean)
    gmmf <- sum(g ^ 2)
    
    return(gmmf)
  }
  
  G <- function(bfun, par){
    G <- numDeriv::jacobian(func = G1, bfun = bfun, x = par)
    return(G)
  }
  
  G1 <- function(bfun, par) {
    G1 <- apply(bfun(par), 2, mean, na.rm=T)
    return(G1)
  }
  
  ## sandwich variance
  var_gmmf <- function(bfun, par){
    bG <- solve(G(bfun, par))
    bg <- bfun(par)
    spsz <- dim(bg)[1]
    Omega <- t(bg) %*% bg / spsz
    Sigma <- bG %*% Omega %*% t(bG) / spsz
    return(Sigma)
  }
  
  VAR_all <- var_gmmf(U_dr, params)
  SE_all <- sqrt(diag(VAR_all))
  
  RR_CI <- exp(log(rrv_dr) + qnorm(c(0.025, 0.975)) * SE_all[1])
  return(list(RR = rrv_dr, RR_CI = RR_CI, VAR = VAR_all))
}


create_folds <- function(n, nfolds = 10) {
  n_resample <- c(sample(n, replace = F), rep(NA, nfolds - n %% nfolds))
  folds_matrix <- matrix(n_resample, nrow = 10)
  
  folds <- vector("list", nfolds)
  
  for (i in 1:nfolds) {
    ff <- folds_matrix[i, ]
    folds[[i]] <- ff[!is.na(ff)]
  }
  return(folds)
}


## estimate RR in the vaccinated using the DR estimator
## GAM with cross-fitting
RRV_DR_NP <- function(I2, V, X, nfolds = 10) {
  
  ## 10-folds cross-fitting
  nn <- length(I2)
  
  Pi_cf <- Psi_cf <- 1:nfolds * NA
  EFUNC <- NULL ## estimating functions
  
  cf_folds <- create_folds(nn, nfolds) 
  len_folds <- unlist(lapply(cf_folds, length))
  
  for (k in 1:nfolds) {
    ## validation data
    I2_val <- I2[cf_folds[[k]]]
    V_val <- V[cf_folds[[k]]]
    X_val <- X[cf_folds[[k]]]
    
    ## training data
    I2_tr <- I2[-cf_folds[[k]]]
    V_tr <- V[-cf_folds[[k]]]
    X_tr <- X[-cf_folds[[k]]]
    
    mu1_data <- data.frame(I2 = I2_tr, X = X_tr)[V_tr == 1,]
    mu0_data <- data.frame(I2 = I2_tr, X = X_tr)[V_tr == 0,]
    pi_data <- data.frame(V = V_tr, X = X_tr)
    
    val_data <- data.frame(I2 = I2_val, X = X_val, V = V_val)
    ## fit the outcome regression models and propensity model
    ## using gam
    mu1_model <- mgcv::gam(I2 ~ s(X), data = mu1_data, family = binomial)
    mu1est <- predict(mu1_model, newdata = val_data, type = "response")
    
    mu0_model <- mgcv::gam(I2 ~ s(X), data = mu0_data, family = binomial)
    mu0est <- predict(mu0_model, newdata = val_data, type = "response")
    
    pi_model <- mgcv::gam(V ~ s(X), data = pi_data, family = binomial)
    piest <- predict(pi_model, newdata = val_data, type = "response")
    
    Pi_cf[k] <- mean((1 - V_val) * (I2_val - mu0est) * piest * (1 - mu1est) /
                       ((1 - piest) * (1 - mu0est) ^ 2) +
                       V_val * (1 - I2_val) * mu0est / (1 - mu0est))
    Psi_cf[k] <- mean(V_val * I2_val) / Pi_cf[k]
    EFUNC <- c(EFUNC,
               (1 - V_val) * (I2_val - mu0est) * piest * (1 - mu1est) /
                 ((1 - piest) * (1 - mu0est) ^ 2) +
                 V_val * (1 - I2_val) * mu0est / (1 - mu0est))
  }
  
  Psi_est <- sum(Psi_cf * len_folds) / sum(len_folds)
  Pi_est <- sum(Pi_cf * len_folds) / sum(len_folds)
  
  ## obtain influence functions
  I2_rs <- I2[unlist(cf_folds)]
  V_rs <- V[unlist(cf_folds)]
  X_rs <- X[unlist(cf_folds)]
  
  EIF <- V_rs * I2_rs / Pi_est - Psi_est / Pi_est * EFUNC
  se_rrv <- sd(EIF) / sqrt(nn)
  
  RR_CI <- exp(log(Psi_est) + qnorm(c(0.025, 0.975)) * se_rrv)
  return(list(RR = Psi_est, RR_CI = RR_CI, SE = se_rrv))
}


# Cohort estimators -------------------------------------------------------

Cohort_OM <- function(I2, T, V, X, U = NULL) {
  if (length(U) != 0) {
    om_model <- glm(I(I2 == 1 & T == 1) ~ V * X + U, family = poisson)
  } else {
    om_model <- glm(I(I2 == 1 & T == 1) ~ V * X, family = poisson)
  }
  
  newdata_V1 <- data.frame(X = X, V = 1)
  newdata_V0 <- data.frame(X = X, V = 0)
  
  mu0 <- predict(om_model, newdata = newdata_V0, type = "response")
  mu1 <- predict(om_model, newdata = newdata_V1, type = "response")
  
  rr_om <- sum(V * mu1) / sum(V * mu0)
  
  ## variables to facilitate variance calculation
  om_coef <- coef(om_model)
  params <- c(log(rr_om), om_coef)
  
  om_matrix <- model.matrix(om_model)
  
  if (length(U) != 0) {
    om_matrix_v1 <- cbind(1, 1, X, U, X)
    om_matrix_v0 <- cbind(1, 0, X, U, X * 0)
  } else {
    om_matrix_v1 <- cbind(1, 1, X, X)
    om_matrix_v0 <- cbind(1, 0, X, X * 0)
  }
  
  # par <- params
  U_om <- function(par) {
    log_rr <- par[1]
    reg_coef <- par[-1]
    
    mu1est <- exp(c(om_matrix_v1 %*% reg_coef))
    mu0est <- exp(c(om_matrix_v0 %*% reg_coef))
    
    U_rr <- V * (mu1est - exp(log_rr) * mu0est)
    U_glm <- om_matrix * (I2 * T - exp(c(om_matrix %*% reg_coef)))
    
    return(cbind(U_rr, U_glm))
  }
  
  GMMF <- function(mrf, par) {
    g0 <- mrf(par = par)
    g <- apply(g0, 2, mean)
    gmmf <- sum(g ^ 2)
    
    return(gmmf)
  }
  
  G <- function(bfun, par){
    G <- numDeriv::jacobian(func = G1, bfun = bfun, x = par)
    return(G)
  }
  
  G1 <- function(bfun, par) {
    G1 <- apply(bfun(par), 2, mean, na.rm=T)
    return(G1)
  }
  
  ## sandwich variance
  var_gmmf <- function(bfun, par){
    bG <- solve(G(bfun, par))
    bg <- bfun(par)
    spsz <- dim(bg)[1]
    Omega <- t(bg) %*% bg / spsz
    Sigma <- bG %*% Omega %*% t(bG) / spsz
    return(Sigma)
  }
  
  VAR_all <- var_gmmf(U_om, params)
  SE_all <- sqrt(diag(VAR_all))
  
  RR_CI <- exp(log(rr_om) + qnorm(c(0.025, 0.975)) * SE_all[1])
  return(list(RR = rr_om, RR_CI = RR_CI, VAR = VAR_all))
}


## estimate RR in the vaccinated using outcome regression
Cohort_DID_OM <- function(I2, I1, T, V, X, exclusive = TRUE) {
  om_model2 <- glm(I(I2==1 & T==1) ~ V * X, family = binomial)
  
  if (exclusive) {
    om_model1 <- glm(I(I1==1 & T==1) ~ V * X, family = binomial)
  } else {
    om_model1 <- glm(I(I1==1 & I2 == 0 & T==1) ~ V * X, family = binomial)
  }
  
  newdata_V1 <- data.frame(X = X, V = 1)
  newdata_V0 <- data.frame(X = X, V = 0)
  
  mu2_0 <- predict(om_model2, newdata = newdata_V0, type = "response")
  mu1_0 <- predict(om_model1, newdata = newdata_V0, type = "response")
  mu1_1 <- predict(om_model1, newdata = newdata_V1, type = "response")
  
  rrv_om <- sum(I2 * T * V) / sum(V * mu2_0 * (mu1_1) / (mu1_0))
  
  ## variables to facilitate variance calculation
  om_coef2 <- coef(om_model2)
  om_coef1 <- coef(om_model1)
  params <- c(log(rrv_om), om_coef2, om_coef1)
  
  om_matrix2 <- model.matrix(om_model2)
  om_matrix1 <- model.matrix(om_model1)
  
  om_matrix_v1 <- cbind(1, 1, X, X)
  om_matrix_v0 <- cbind(1, 0, X, X * 0)
  
  # par <- params
  U_om <- function(par) {
    log_rrv <- par[1]
    reg_coef2 <- par[-1][1:(length(par - 1) / 2)]
    reg_coef1 <- par[-1][(length(par - 1) / 2 + 1):length(par - 1)]
    
    mu2_0est <- expit(c(om_matrix_v0 %*% reg_coef2))
    mu1_1est <- expit(c(om_matrix_v1 %*% reg_coef1))
    mu1_0est <- expit(c(om_matrix_v0 %*% reg_coef1))
    
    U_rrv <- V * (I2 * T - exp(log_rrv) * mu2_0est *
                    (mu1_1est) / (mu1_0est))
    
    U_glm2 <- om_matrix2 * (I2 * T - expit(c(om_matrix2 %*% reg_coef2)))
    if (exclusive) {
      U_glm1 <- om_matrix1 * (I1 * T - expit(c(om_matrix1 %*% reg_coef1)))
    } else {
      U_glm1 <- om_matrix1 * (I1 * (1 - I2) * T - expit(c(om_matrix1 %*% reg_coef1)))
    }
    
    return(cbind(U_rrv, U_glm2, U_glm1))
  }
  
  GMMF <- function(mrf, par) {
    g0 <- mrf(par = par)
    g <- apply(g0, 2, mean)
    gmmf <- sum(g ^ 2)
    
    return(gmmf)
  }
  
  G <- function(bfun, par){
    G <- numDeriv::jacobian(func = G1, bfun = bfun, x = par)
    return(G)
  }
  
  G1 <- function(bfun, par) {
    G1 <- apply(bfun(par), 2, mean, na.rm=T)
    return(G1)
  }
  
  ## sandwich variance
  var_gmmf <- function(bfun, par){
    bG <- solve(G(bfun, par))
    bg <- bfun(par)
    spsz <- dim(bg)[1]
    Omega <- t(bg) %*% bg / spsz
    Sigma <- bG %*% Omega %*% t(bG) / spsz
    return(Sigma)
  }
  
  VAR_all <- var_gmmf(U_om, params)
  SE_all <- sqrt(diag(VAR_all))
  
  RR_CI <- exp(log(rrv_om) + qnorm(c(0.025, 0.975)) * SE_all[1])
  return(list(RR = rrv_om, RR_CI = RR_CI, VAR = VAR_all))
}

