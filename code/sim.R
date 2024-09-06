runsim <- function(dt, p, true_rrv, exclusive = TRUE) {
  p()
  
  I2 <- dt$I2
  I1 <- dt$I1
  T <- dt$T
  V <- dt$V
  X <- dt$X
  U <- dt$U
  
  cohort_reg_U_res <- Cohort_OM(I2, T, V, X, U)
  cohort_reg_U_bias <- cohort_reg_U_res$RR - true_rrv
  cohort_reg_U_cover <- cohort_reg_U_res$RR_CI[1] < true_rrv & 
    cohort_reg_U_res$RR_CI[2] > true_rrv
  cohort_reg_U_len <- abs(cohort_reg_U_res$RR_CI[2] - cohort_reg_U_res$RR_CI[1])
     
  
  cohort_reg_noU_res <- Cohort_OM(I2, T, V, X)
  cohort_reg_noU_bias <- cohort_reg_noU_res$RR - true_rrv
  cohort_reg_noU_cover <- cohort_reg_noU_res$RR_CI[1] < true_rrv & 
    cohort_reg_noU_res$RR_CI[2] > true_rrv
  cohort_reg_noU_len <- abs(cohort_reg_noU_res$RR_CI[2] - cohort_reg_noU_res$RR_CI[1])
  
  did_reg_res <- Cohort_DID_OM(I2, I1, T, V, X, exclusive)
  did_reg_bias <- did_reg_res$RR - true_rrv
  did_reg_cover <- did_reg_res$RR_CI[1] < true_rrv & 
    did_reg_res$RR_CI[2] > true_rrv
  did_reg_len <- abs(did_reg_res$RR_CI[2] - did_reg_res$RR_CI[1])
  
  tnd <- dt[dt$S == 1, ]
  
  I2 <- tnd$I2
  V <- tnd$V
  X <- tnd$X
  
  logit_reg_res <- LogitReg(I2, V, X)
  logit_reg_bias <- logit_reg_res$RR - true_rrv
  logit_reg_cover <- logit_reg_res$RR_CI[1] < true_rrv & 
    logit_reg_res$RR_CI[2] > true_rrv
  logit_reg_len <- abs(logit_reg_res$RR_CI[2] - logit_reg_res$RR_CI[1])
    
  rrv_om_res <- RRV_OM(I2, V, X)
  rrv_om_bias <- rrv_om_res$RR - true_rrv
  rrv_om_cover <- rrv_om_res$RR_CI[1] < true_rrv & 
    rrv_om_res$RR_CI[2] > true_rrv
  rrv_om_len <- abs(rrv_om_res$RR_CI[2] - rrv_om_res$RR_CI[1])
  
  rrv_ipw_res <- RRV_IPW(I2, V, X)
  rrv_ipw_bias <- rrv_ipw_res$RR - true_rrv
  rrv_ipw_cover <- rrv_ipw_res$RR_CI[1] < true_rrv & 
    rrv_ipw_res$RR_CI[2] > true_rrv
  rrv_ipw_len <- abs(rrv_ipw_res$RR_CI[2] - rrv_ipw_res$RR_CI[1])
  
  rrv_dr_res <- RRV_DR(I2, V, X)
  rrv_dr_bias <- rrv_dr_res$RR - true_rrv
  rrv_dr_cover <- rrv_dr_res$RR_CI[1] < true_rrv & 
    rrv_dr_res$RR_CI[2] > true_rrv
  rrv_dr_len <- abs(rrv_dr_res$RR_CI[2] - rrv_dr_res$RR_CI[1])
  
  return(c(#true_rrv = true_rrv,
    logit_reg_bias = as.numeric(logit_reg_bias),
    logit_reg_cover = as.numeric(logit_reg_cover),
    logit_reg_len = as.numeric(logit_reg_len),
    rrv_om_bias = as.numeric(rrv_om_bias),
    rrv_om_cover = as.numeric(rrv_om_cover),
    rrv_om_len = as.numeric(rrv_om_len),
    rrv_ipw_bias = as.numeric(rrv_ipw_bias),
    rrv_ipw_cover = as.numeric(rrv_ipw_cover),
    rrv_ipw_len = as.numeric(rrv_ipw_len),
    rrv_dr_bias = as.numeric(rrv_dr_bias),
    rrv_dr_cover = as.numeric(rrv_dr_cover),
    rrv_dr_len = as.numeric(rrv_dr_len),
    did_reg_bias = as.numeric(did_reg_bias),
    did_reg_cover = as.numeric(did_reg_cover),
    did_reg_len = as.numeric(did_reg_len),
    cohort_reg_U_bias = as.numeric(cohort_reg_U_bias),
    cohort_reg_U_cover = as.numeric(cohort_reg_U_cover),
    cohort_reg_U_len = as.numeric(cohort_reg_U_len),
    cohort_reg_noU_bias = as.numeric(cohort_reg_noU_bias),
    cohort_reg_noU_cover = as.numeric(cohort_reg_noU_cover),
    cohort_reg_noU_len = as.numeric(cohort_reg_noU_len),
    sample_size = as.numeric(nrow(tnd))))
}
