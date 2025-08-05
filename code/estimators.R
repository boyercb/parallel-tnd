
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
  
  #rrv_om <- sum(I2 * V) / sum(V * mu0 * (1 - mu1) / (1 - mu0))
  rrv_om <- sum(I2 * V) / sum(V * (1 - I2) * mu0 / (1 - mu0))
  
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
    
    U_rrv <- V * (I2 - exp(log_rrv) * (1 - I2) * mu0est / (1 - mu0est))
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
RRV_IPW <- function(I2, V, X, family = "binomial") {
  
  if (family == "binomial") {
    ps_model <- glm(V ~ I2 * X, family = binomial(link = "logit"))
  } else if (family == "poisson") {
    ps_model <- glm(V ~ I2 * X, family = poisson(link = "log"))
  }
  
  newdata_I0 <- data.frame(X = X, I2 = 0)
  
  pi0 <- predict(ps_model, newdata = newdata_I0, type = "response")
  
  rrv_ipw <- sum(I2 * V) / sum((1 - V) * I2 * pi0 / (1 - pi0))
  
  ## variables to facilitate variance calculation
  ps_coef <- coef(ps_model)
  params <- c(log(rrv_ipw), ps_coef)
  
  ps_matrix <- model.matrix(ps_model)
  #ps_matrix_i0 <- cbind(1, X)
  ps_matrix_i0 <- cbind(1, 0, X, X * 0)
  
  # par <- params
  U_ipw <- function(par) {
    log_rrv <- par[1]
    reg_coef <- par[-1]
    
    if (family == "binomial") {
      pi0est <- expit(c(ps_matrix_i0 %*% reg_coef))# * (I2 == 0)
    } else if (family == "poisson") {
      pi0est <- exp(c(ps_matrix_i0 %*% reg_coef))# * (I2 == 0)
    }
    
    U_rrv <- V * I2 - exp(log_rrv) * (1 - V) * I2 * pi0est / (1 - pi0est)
    
    if (family == "binomial") {
      U_glm <- ps_matrix * (V - expit(c(ps_matrix %*% reg_coef)))# * (I2 == 0)
    } else if (family == "poisson") {
      U_glm <- ps_matrix * (V - exp(c(ps_matrix %*% reg_coef)))# * (I2 == 0)
    }
    
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


RRV_DR <- function(I2, V, X) {
  
  ps_model <- glm(V ~ I2 * X, family = binomial)
  om_model <- glm(I2 ~ V * X, family = binomial)
  
  newdata_X <- data.frame(X = X)
  newdata_V0 <- data.frame(X = X, V = 0)
  newdata_I0 <- data.frame(X = X, I2 = 0)
  
  pi <- predict(ps_model, newdata = newdata_I0, type = "response")
  mu <- predict(om_model, newdata = newdata_V0, type = "response")
  
  eqn <- function(theta) {
    eqn <- 
      cbind(1, X) * (V - pi) * exp(-(theta[1] + theta[2] * X) * V * I2) * (I2 - mu)
      #(V - pi) * exp(-(theta) * V * I2) * (I2 - mu)
    
    apply(eqn, 2, function(x) mean(x))
    #mean(eqn)
  }
  est <- rootSolve::multiroot(eqn, c(0,0))
  psi <- est$root

  rrv_dr <- sum(I2 * V) / 
    sum(I2 * exp(-(psi[1] + psi[2] * X)) * V)
  
  # variables to facilitate variance calculation
  ps_coef <- coef(ps_model); n_ps_coef <- length(ps_coef)
  om_coef <- coef(om_model); n_om_coef <- length(om_coef)
  n_psi <- length(psi)
  params <- c(log(rrv_dr), psi, ps_coef, om_coef)

  om_matrix <- model.matrix(om_model)
  om_matrix_v0 <- cbind(1, 0, X, X * 0)

  ps_matrix <- model.matrix(ps_model)
  ps_matrix_i0 <- cbind(1, 0, X, X * 0)
  
  or_matrix <- cbind(1, X)

  # par <- params
  U_dr <- function(par) {
    log_rrv <- par[1]
    psi_coef <- par[2:(1 + n_psi)]
    orest <- #exp(psi_coef)
      exp(c(or_matrix %*% psi_coef))
    ps_coef <- par[(2 + n_psi):(1 + n_psi + n_ps_coef)]
    pixest <- expit(c(ps_matrix_i0 %*% ps_coef))

    om_coef <- par[(2 + n_psi + n_ps_coef):(1 + n_psi + n_ps_coef + n_om_coef)]
    mu0est <- expit(c(om_matrix_v0 %*% om_coef))

    U_rrv <- V * I2 - exp(log_rrv) * (V * I2 / orest)
    U_or <- #(V - pixest) * exp(-psi_coef * V * I2) * (I2 - mu0est)
      or_matrix * (V - pixest) * exp(-c(or_matrix %*% psi_coef) * V * I2) * (I2 - mu0est)
    U_om <- om_matrix * (I2 - expit(c(om_matrix %*% om_coef)))
    U_ps <- ps_matrix * (V - expit(c(ps_matrix %*% ps_coef)))

    return(cbind(U_rrv, U_or, U_ps, U_om))
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

## estimate RR in the vaccinated using an NCE
RRV_NCE <- function(I2, V, X, Z) {
  stage1_model <- glm(I(I2 == 0) ~ V + X + Z, family = poisson)
  mu0 <- predict(stage1_model, type = "link")
  stage2_model <- glm(I2 ~ V + X + mu0, family = binomial)

  rrv_nce <- exp(coef(stage2_model)[2])
  
  # ## variables to facilitate variance calculation
  # stage1_coef <- coef(stage1_model); n_stage1_coef <- length(stage1_coef)
  # stage2_coef <- coef(stage2_model); n_stage2_coef <- length(stage2_coef)
  # params <- c(stage2_coef, stage1_coef)
  # 
  # stage1_matrix <- model.matrix(stage1_model)
  # stage2_matrix <- model.matrix(stage2_model)
  # 
  # # par <- params
  # U_nce <- function(par) {
  #   stage2_coef <- par[1:n_stage2_coef]
  #   stage1_coef <- par[(n_stage2_coef + 1):(n_stage1_coef + n_stage2_coef)]
  #   stage2_matrix[, 4] <- stage1_coef[4]
  #   
  #   U_glm1 <- stage1_matrix * ((1 - I2) - exp(c(stage1_matrix %*% stage1_coef)))
  #   U_glm2 <- stage2_matrix * (I2 - expit(c(stage2_matrix %*% stage2_coef)))
  #   
  #   return(cbind(U_glm2, U_glm1))
  # }
  # 
  # GMMF <- function(mrf, par) {
  #   g0 <- mrf(par = par)
  #   g <- apply(g0, 2, mean)
  #   gmmf <- sum(g ^ 2)
  #   
  #   return(gmmf)
  # }
  # 
  # G <- function(bfun, par){
  #   G <- numDeriv::jacobian(func = G1, bfun = bfun, x = par)
  #   return(G)
  # }
  # 
  # G1 <- function(bfun, par) {
  #   G1 <- apply(bfun(par), 2, mean, na.rm=T)
  #   return(G1)
  # }
  # 
  # ## sandwich variance
  # var_gmmf <- function(bfun, par){
  #   bG <- solve(G(bfun, par))
  #   bg <- bfun(par)
  #   spsz <- dim(bg)[2]
  #   Omega <- t(bg) %*% bg / spsz
  #   Sigma <- bG %*% Omega %*% t(bG) / spsz
  #   return(Sigma)
  # }
  # 
  # VAR_all <- var_gmmf(U_nce, params)
  # SE_all <- sqrt(diag(VAR_all))
  # 
  # RR_CI <- exp(log(rrv_nce) + qnorm(c(0.025, 0.975)) * SE_all[2])
  return(list(RR = rrv_nce, RR_CI = NA, VAR = NA))#RR_CI, VAR = VAR_all))
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
Cohort_DID_OM <- function(I2, I1, T, Y, V, X, exclusive = TRUE) {
  om_model2 <- glm(I(I2==1 & T==1) ~ V + X, family = binomial)
  
  if (exclusive) {
    om_model1 <- glm(I(I1==1 & T==1) ~ V + X, family = binomial)
  } else {
    om_model1 <- glm(I(I1==1 & I2 == 0 & T==1) ~ V + X, family = binomial)
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
  
  om_matrix_v1 <- cbind(1, 1, X)
  om_matrix_v0 <- cbind(1, 0, X)
  
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

## estimate RR in the vaccinated using inverse-probability weighting
Cohort_DID_IPW <- function(I2, I1, T, Y, V, X, exclusive = TRUE) {
  ps_model <- glm(V ~ X + I(Y==1), family = binomial, subset = Y != 2)
  newdata_I1 <- data.frame(X = X, Y = 1)
  
  pi0 <- predict(ps_model, newdata = newdata_I1, type = "response")
  
  rrv_ipw <- sum(I(Y==2) * V) / sum((1 - V) * I(Y==2) * pi0 / (1 - pi0))
  
  ## variables to facilitate variance calculation
  ps_coef <- coef(ps_model)
  params <- c(log(rrv_ipw), ps_coef)
  
  ps_matrix <- cbind(1, X, as.numeric(Y==1))
  ps_matrix_i0 <- cbind(1, X, 1)
 
  # par <- params
  U_ipw <- function(par) {
    log_rrv <- par[1]
    reg_coef <- par[-1]
    
    pi0est <- expit(c(ps_matrix_i0 %*% reg_coef))
    
    U_rrv <- V * I2 - exp(log_rrv) * (1 - V) * I2 * pi0est / (1 - pi0est)
    U_glm <- (ps_matrix * (V - expit(c(ps_matrix %*% reg_coef)))) * (Y != 2)
    
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


Cohort_DID_DR1 <- function(I2, I1, T, Y, V, X, exclusive = TRUE, multinomial = FALSE, recalc = FALSE, link = "logit") {
  if (!multinomial) {
    om_model2 <- glm(I(Y==2) ~ X, family = binomial(link), subset = V == 0)
    om_model11 <- glm(I(Y==1) ~ X, family = binomial(link), subset = V == 0)
    om_model10 <- glm(I(Y==1) ~ X, family = binomial(link), subset = V == 0)
    
    ps_model <- glm(V ~ X, family = binomial, subset = Y==1)
    
    newdata_V1 <- data.frame(X = X, V = 1)
    newdata_V0 <- data.frame(X = X, V = 0)
    newdata_I1 <- data.frame(X = X, Y = 1)
    newdata_I0 <- data.frame(X = X, Y = 0)
    
    pi0 <- predict(ps_model, newdata = newdata_I1, type = "link")
    
    mu2_0 <- predict(om_model2, newdata = newdata_V0, type = "response")
    mu1_0 <- predict(om_model10, newdata = newdata_V0, type = "response")
    mu1_1 <- predict(om_model11, newdata = newdata_V1, type = "response")
  } else {
    om_model <- nnet::multinom(factor(Y)~X, subset = V==0)
    om_model1 <- nnet::multinom(factor(Y)~X, subset = V==1)
    
    ps_model <- glm(V ~ X, family = binomial, subset = Y==1)
    
    newdata_V1 <- data.frame(X = X, V = 1)
    newdata_V0 <- data.frame(X = X, V = 0)
    newdata_I1 <- data.frame(X = X, Y = 1)
    newdata_I0 <- data.frame(X = X, Y = 0)
    
    pi0 <- predict(ps_model, newdata = newdata_I1, type = "link")

    muY0 <- predict(om_model, newdata = newdata_V0, type = "probs")
    mu2_0 <- muY0[, 3]
    mu1_0 <- muY0[, 2]
    mu0_0 <- muY0[, 1]
    muY1 <- predict(om_model1, newdata = newdata_V1, type = "probs")
    mu1_1 <- muY1[, 2]
  }
  
  
  eqn <- function(theta) {
    eqn <- cbind(1, X) *
      ((1 - V) * ((1 + exp((theta[1] + theta[2] * X) * I(Y==0) + pi0 * I(Y!=0)))) - 1)
    
    apply(eqn, 2, function(x) mean(x))
  }
  est <- rootSolve::multiroot(eqn, c(0, 0))
  
  eta <- est$root

  eqn <- function(theta) {
    eqn <- 
      (V - expit(eta[1] + eta[2] * X)) * exp(-(0 * I(Y!=1) + (theta) * I(Y==1)) * V) *
      (I(Y!=1) * (0) + I(Y==1) * (1 - mu1_0))
    
    mean(eqn)
  }
  est <- rootSolve::multiroot(eqn, c(0))
  
  beta <- est$root
  beta_ipw <- mean(pi0 - eta[1] - eta[2] * X)

  if (recalc) {
    eqn <- function(theta) {
      eqn <- cbind(1, X) *
        ((1 - V) * ((1 + exp(theta[1] + theta[2] * X + (beta) * I(Y!=0)))) - 1)
      
      apply(eqn, 2, function(x) mean(x))
    }
    est <- rootSolve::multiroot(eqn, c(0, 0))
    
    eta <- est$root
  }
  
  den <- exp(0) * (1 - mu1_0 - mu2_0) + exp(beta) * mu1_0 + 
    exp(beta) * mu2_0
  
  rrv_dr <- sum(I(Y==2) * V) / 
    sum((1 - V) * exp(eta[1] + eta[2] * X + (beta) * I(Y!=0)) * 
          (I(Y==2) - exp((beta)) / den * mu2_0) + V * exp((beta)) / den * mu2_0)
  
  return(list(RR = rrv_dr, RR_CI = NA, VAR = NA))
}

Cohort_DID_DRX <- function(I2, I1, T, Y, V, X, exclusive = TRUE, multinomial = FALSE, recalc = FALSE, link = "logit") {
  if (!multinomial) {
    om_model2 <- glm(I(Y==2) ~ X, family = binomial(link), subset = V == 0)
    om_model11 <- glm(I(Y==1) ~ X, family = binomial(link), subset = V == 0)
    om_model10 <- glm(I(Y==1) ~ X, family = binomial(link), subset = V == 0)
    
    ps_model <- glm(V ~ X, family = binomial, subset = Y==1)
    
    newdata_V1 <- data.frame(X = X, V = 1)
    newdata_V0 <- data.frame(X = X, V = 0)
    newdata_I1 <- data.frame(X = X, Y = 1)
    newdata_I0 <- data.frame(X = X, Y = 0)
    
    pi0 <- predict(ps_model, newdata = newdata_I1, type = "link")
    
    mu2_0 <- predict(om_model2, newdata = newdata_V0, type = "response")
    mu1_0 <- predict(om_model10, newdata = newdata_V0, type = "response")
    mu1_1 <- predict(om_model11, newdata = newdata_V1, type = "response")
  } else {
    om_model <- nnet::multinom(factor(Y)~X, subset = V==0)
    om_model1 <- nnet::multinom(factor(Y)~X, subset = V==1)
    
    ps_model <- glm(V ~ X, family = binomial, subset = Y==1)
    
    newdata_V1 <- data.frame(X = X, V = 1)
    newdata_V0 <- data.frame(X = X, V = 0)
    newdata_I1 <- data.frame(X = X, Y = 1)
    newdata_I0 <- data.frame(X = X, Y = 0)
    
    pi0 <- predict(ps_model, newdata = newdata_I1, type = "link")
    
    muY0 <- predict(om_model, newdata = newdata_V0, type = "probs")
    mu2_0 <- muY0[, 3]
    mu1_0 <- muY0[, 2]
    mu0_0 <- muY0[, 1]
    muY1 <- predict(om_model1, newdata = newdata_V1, type = "probs")
    mu1_1 <- muY1[, 2]
  }
  
  
  eqn <- function(theta) {
    eqn <- cbind(1, X) *
      ((1 - V) * ((1 + exp((theta[1] + theta[2] * X) * I(Y==0) + pi0 * I(Y!=0)))) - 1)
    
    apply(eqn, 2, function(x) mean(x))
  }
  est <- rootSolve::multiroot(eqn, c(0, 0))
  
  eta <- est$root
  
  eqn <- function(theta) {
    eqn <- 
      (V - expit(eta[1] + eta[2] * X)) * exp(-(0 * I(Y!=1) + (theta * X) * I(Y==1)) * V) *
      (I(Y!=1) * (0) + I(Y==1) * (1 - mu1_0))
    
    mean(eqn)
  }
  est <- rootSolve::multiroot(eqn, c(0))
  
  beta <- est$root
  beta_ipw <- mean(pi0 - eta[1] - eta[2] * X)
  
  if (recalc) {
    eqn <- function(theta) {
      eqn <- cbind(1, X) *
        ((1 - V) * ((1 + exp(theta[1] + theta[2] * X + (beta * X) * I(Y!=0)))) - 1)
      
      apply(eqn, 2, function(x) mean(x))
    }
    est <- rootSolve::multiroot(eqn, c(0, 0))
    
    eta <- est$root
  }
  
  den <- exp(0) * (1 - mu1_0 - mu2_0) + exp(beta * X) * mu1_0 + 
    exp(beta * X) * mu2_0
  
  rrv_dr <- sum(I(Y==2) * V) / 
    sum((1 - V) * exp(eta[1] + eta[2] * X + (beta * X) * I(Y!=0)) * 
          (I(Y==2) - exp((beta * X)) / den * mu2_0) + V * exp((beta * X)) / den * mu2_0)
  
  return(list(RR = rrv_dr, RR_CI = NA, VAR = NA))
}


Cohort_DID_DR2 <- function(I2, I1, T, Y, V, X, exclusive = TRUE, multinomial = FALSE, recalc = FALSE, link = "logit") {
  if (!multinomial) {
    om_model2 <- glm(I(Y==2) ~ X, family = binomial(link), subset = V == 0)
    om_model11 <- glm(I(Y==1) ~ X, family = binomial(link), subset = V == 0)
    om_model10 <- glm(I(Y==1) ~ X, family = binomial(link), subset = V == 0)
    
    ps_model <- glm(V ~ X, family = binomial, subset = Y==1)
    
    newdata_V1 <- data.frame(X = X, V = 1)
    newdata_V0 <- data.frame(X = X, V = 0)
    newdata_I1 <- data.frame(X = X, Y = 1)
    newdata_I0 <- data.frame(X = X, Y = 0)
    
    pi0 <- predict(ps_model, newdata = newdata_I1, type = "link")
    
    mu2_0 <- predict(om_model2, newdata = newdata_V0, type = "response")
    mu1_0 <- predict(om_model10, newdata = newdata_V0, type = "response")
    mu1_1 <- predict(om_model11, newdata = newdata_V1, type = "response")
  } else {
    om_model <- nnet::multinom(factor(Y)~X, subset = V==0)
    om_model1 <- nnet::multinom(factor(Y)~X, subset = V==1)
    
    ps_model <- glm(V ~ X, family = binomial, subset = Y==1)
    
    newdata_V1 <- data.frame(X = X, V = 1)
    newdata_V0 <- data.frame(X = X, V = 0)
    newdata_I1 <- data.frame(X = X, Y = 1)
    newdata_I0 <- data.frame(X = X, Y = 0)
    
    pi0 <- predict(ps_model, newdata = newdata_I1, type = "link")
    
    muY0 <- predict(om_model, newdata = newdata_V0, type = "probs")
    mu2_0 <- muY0[, 3]
    mu1_0 <- muY0[, 2]
    mu0_0 <- muY0[, 1]
    muY1 <- predict(om_model1, newdata = newdata_V1, type = "probs")
    mu1_1 <- muY1[, 2]
  }
  
  
  eqn <- function(theta) {
    eqn <- cbind(1, X) *
      ((1 - V) * ((1 + exp((theta[1] + theta[2] * X) * I(Y==0) + pi0 * I(Y!=0)))) - 1)
    
    apply(eqn, 2, function(x) mean(x))
  }
  est <- rootSolve::multiroot(eqn, c(0, 0))
  
  eta <- est$root
  
  eqn <- function(theta) {
    eqn <- cbind(1, X) *
      (V - expit(eta[1] + eta[2] * X)) * exp(-(0 * I(Y!=1) + (theta[1] + theta[2] * X) * I(Y==1)) * V) *
      (I(Y!=1) * (0) + I(Y==1) * (1 - mu1_0))
    
    apply(eqn, 2, function(x) mean(x))
  }
  est <- rootSolve::multiroot(eqn, c(0, 0))
  
  beta <- est$root
  beta_ipw <- mean(pi0 - eta[1] - eta[2] * X)
  
  if (recalc) {
    eqn <- function(theta) {
      eqn <- cbind(1, X) *
        ((1 - V) * ((1 + exp(theta[1] + theta[2] * X + (beta[1] + beta[2] * X) * I(Y!=0)))) - 1)
      
      apply(eqn, 2, function(x) mean(x))
    }
    est <- rootSolve::multiroot(eqn, c(0, 0))
    
    eta <- est$root
  }
  
  den <- exp(0) * (1 - mu1_0 - mu2_0) + exp(beta[1] + beta[2] * X) * mu1_0 + 
    exp(beta[1] + beta[2] * X) * mu2_0
  
  rrv_dr <- sum(I(Y==2) * V) / 
    sum((1 - V) * exp(eta[1] + eta[2] * X + (beta[1] + beta[2] * X) * I(Y!=0)) * 
          (I(Y==2) - exp((beta[1] + beta[2] * X)) / den * mu2_0) + V * exp((beta[1] + beta[2] * X)) / den * mu2_0)
  
  return(list(RR = rrv_dr, RR_CI = NA, VAR = NA))
}

## estimate RR in the vaccinated using inverse-probability weighting
Cohort_DID_DR <- function(I2, I1, T, Y, V, X, exclusive = TRUE) {
  # om_model2 <- glm(I(Y==2) ~ V + X, family = binomial)
  # 
  # if (exclusive) {
  #   om_model1 <- glm(I(Y==1) ~ V + X, family = binomial)
  # } else {
  #   om_model1 <- glm(I(I1==1 & I2 == 0 & T==1) ~ V + X, family = binomial)
  # }
  om_model <- nnet::multinom(factor(Y)~X, subset = V==0)
  om_model1 <- nnet::multinom(factor(Y)~X, subset = V==1)
  # ps_model <- glm(V ~ X + I(Y==1), family = binomial, subset = Y!=2)
  ps_model <- glm(V ~ X, family = binomial, subset = Y==1)
  
  newdata_V1 <- data.frame(X = X, V = 1)
  newdata_V0 <- data.frame(X = X, V = 0)
  newdata_I1 <- data.frame(X = X, Y = 1)
  newdata_I0 <- data.frame(X = X, Y = 0)
  
  pi0 <- predict(ps_model, newdata = newdata_I1, type = "link")
  #nu <- predict(ps_model, newdata = newdata_I0, type = "response")
  
  muY0 <- predict(om_model, newdata = newdata_V0, type = "probs")
  mu2_0 <- muY0[, 3]
  mu1_0 <- muY0[, 2]
  mu0_0 <- muY0[, 1]
  muY1 <- predict(om_model1, newdata = newdata_V1, type = "probs")
  mu1_1 <- muY1[, 2]
  # mu2_0 <- predict(om_model2, newdata = newdata_V0, type = "response")
  # mu1_0 <- predict(om_model1, newdata = newdata_V0, type = "response")
  # mu1_1 <- predict(om_model1, newdata = newdata_V1, type = "response")
  
  eqn <- function(theta) {
    eqn <- cbind(1, X) *
      ((1 - V) * ((1 + exp((theta[1] + theta[2] * X) * I(Y==0) + pi0 * I(Y!=0)))) - 1)
    
    apply(eqn, 2, function(x) mean(x))
  }
  est <- rootSolve::multiroot(eqn, c(0, 0))
  
  eta <- est$root
  print(paste0("eta: ", round(eta, 2)))
  
  eqn <- function(theta) {
    eqn <- cbind(1, X) *
      (V - expit(eta[1] + eta[2] * X)) * exp(-(0 * I(Y!=1) + (theta[1] + theta[2] * X) * I(Y==1)) * V) *
         (I(Y!=1) * (0) + I(Y==1) * (1 - mu1_0))

    apply(eqn, 2, function(x) mean(x))
  }
  est <- rootSolve::multiroot(eqn, c(0, 0))
  
  beta <- est$root
  print(paste0("beta (dr): ", round(mean(beta[1] + beta[2] * X), 4)))
  beta_ipw <- mean(pi0 - eta[1] - eta[2] * X)
  print(paste0("beta (ipw): ", round(beta_ipw, 4)))
  
  eqn <- function(theta) {
    eqn <- cbind(1, X) *
      ((1 - V) * ((1 + exp(theta[1] + theta[2] * X + (beta[1] + beta[2] * X) * I(Y!=0)))) - 1)

    apply(eqn, 2, function(x) mean(x))
  }
  est <- rootSolve::multiroot(eqn, c(0, 0))

  eta <- est$root
  print(paste0("eta (2nd): ", round(eta, 2)))
  
  # ALT
  # print(paste0("eta (ipw): ", round(coef(ps_model)[c(1,2)], 4)))
  # print(paste0("beta (ipw): ", round(coef(ps_model)[3], 4)))
  # beta_ipw <- coef(ps_model)[3]
  
  
  # eqn <- function(theta) {
  #   eqn <- 
  #     I(Y!=2) * ((V - nu) * exp(-(0 * I(Y!=1) + theta * I(Y==1)) * V) * 
  #        (I(Y!=1) * (0) + I(Y==1) * (1 - mu1_0)))# + I(Y==2) * (1 - mu2_0)))
  #   
  #   sum(eqn)
  # }
  # est <- rootSolve::multiroot(eqn, c(0))
  # 
  # beta <- est$root
  # print(paste0("beta (dr): ", round(beta, 4)))
  # 
  # eqn <- function(theta) {
  #   eqn <- cbind(1, X) *
  #     ((1 - V) * ((1 + exp(theta[1] + theta[2] * X + beta * I(Y!=0)))) - 1)
  # 
  #   apply(eqn, 2, function(x) mean(x))
  # }
  # est <- rootSolve::multiroot(eqn, c(0, 0))
  # 
  # eta <- est$root
  # print(paste0("eta: ", round(eta, 2)))

  #beta <- beta_ipw
  #den <- mean((1 - V) * exp(beta)) / mean(1 - V)
  #den <- mean((1 - V) * exp(0 * I(Y==0) + beta * I(Y==1) + beta * I(Y==2))) / mean(1 - V)
  den <- exp(0) * (1 - mu1_0 - mu2_0) + exp(beta[1] + beta[2] * X) * mu1_0 + 
    exp(beta[1] + beta[2] * X) * mu2_0
  # den <- mean(exp(0 * I(Y[V==0]==0) + (beta[1] + beta[2] * X) * I(Y[V==0]==1) + (beta[1] + beta[2] * X) * I(Y[V==0]==2)))
  den_ipw <- mean((1 - V) * exp(0 * I(Y==0) + beta_ipw * I(Y==1) + beta_ipw * I(Y==2))) / mean(1 - V)
  
  # print(paste0("beta: ", round(beta, 2)))
  print(paste0("den (ipw): ", round(den_ipw, 4)))
  print(paste0("den (dr): ", round(mean(den), 4)))
  
  print(paste0("alpha (ipw): ", round(exp(beta_ipw) / den_ipw, 2)))
  print(paste0("alpha (dr): ", round(mean(exp(beta[1] + beta[2] * X) / den), 2)))
  print(paste0("alpha (om): ", round(mean(mu1_1/mu1_0), 2)))
  
  df0 <- data.frame(X = seq(0, 1, 0.01), V = 0)
  df1 <- data.frame(X = seq(0, 1, 0.01), V = 1)
  # mu2_0p <- predict(om_model2, newdata = df0, type = "response")
  # mu1_0p <- predict(om_model1, newdata = df0, type = "response")
  # mu1_1p <- predict(om_model1, newdata = df1, type = "response")
  muY0 <- predict(om_model, newdata = df0, type = "probs")
  mu2_0p <- muY0[, 3]
  mu1_0p <- muY0[, 2]
  mu0_0p <- muY0[, 1]
  muY1 <- predict(om_model1, newdata = df1, type = "probs")
  mu1_1p <- muY1[, 2]
  mu0_0p <- muY0[, 1]
  plot(seq(0, 1, 0.01), mu1_1p/mu1_0p * (exp(0) * (1 - mu1_0p - mu2_0p) + 
                                           exp(beta[1] + beta[2] * seq(0, 1, 0.01)) * mu1_0p +
                                           exp(beta[1] + beta[2] * seq(0, 1, 0.01)) * mu2_0p))
  plot(seq(0, 1, 0.01), exp(beta[1] + beta[2] * seq(0, 1, 0.01)))
  
  rrv_dr <- sum(I(Y==2) * V) / 
    sum((1 - V) * exp(eta[1] + eta[2] * X + (beta[1] + beta[2] * X) * I(Y!=0)) * 
          (I(Y==2) - exp((beta[1] + beta[2] * X)) / den * mu2_0) + V * exp((beta[1] + beta[2] * X)) / den * mu2_0)
  print(paste0("rrv (om): ", sum(I(Y==2) * V) / 
                 sum(V * exp((beta[1] + beta[2] * X)) / den * mu2_0)))
  return(list(RR = rrv_dr, RR_CI = NA, VAR = NA))
  
}
