# 文件路径: R/main.R

#' Run TLMR Estimation
#'
#' This is the main function that runs the full TLMR procedure.
#'
#' @param exposure_data The exposure data. The data should be in list format, with variable names specified as follows: instrumental variables as Z, exposure as X, outcome as Y, and covariates as C.
#' @param outcome_data The outcome data. The data should be in list format, with variable names specified as follows: instrumental variables as Z, exposure as X, outcome as Y, and covariates as C.
#' @param p_threshold The setting of the p-value threshold used to test the significance of effect modification.The default is 0.005.
#' @param covariates_type 0 represents continuous type, 1 represents classified type. The default is 0.
#' @param summary_cov Covariates that need to be adjusted for calculating summary data, vector format(i.e. c(1,2)). The default is NULL.
#' @param n_divides When using IPW or AIPW estimation, this specifies the number of categories used to discretize continuous covariates. The default is 3.
#' @return A list with results.
#' @export
#'
run_TLMR <- function(exposure_data, outcome_data, p_threshold, covariates_type = 0, summary_cov =NULL
                     , n_divides =3) {

  # exposure_data <- data.table(exposure_data)
  # outcome_data <- data.table(outcome_data)

  Z_1 <- exposure_data$Z
  Z_0 <- outcome_data$Z

  X_1 <- exposure_data$X
  Y_0 <- outcome_data$Y

  C_1 <- exposure_data$C
  C_0<- outcome_data$C

  n_0 <- dim(Z_0)[1]
  n_1 <- dim(Z_1)[1]

  p_Z <- dim(Z_0)[2]
  p_C <- dim(C_0)[2]

  AS_0 <- C_0[,summary_cov]
  AS_1 <- C_1[,summary_cov]

  quantile_vec_all <- seq(0,1,1/n_divides)
  quantile_vec <- quantile_vec_all[-c(1,length(quantile_vec_all))]
  #------------------------------
  C_1 <- as.data.frame(C_1)
  C_1_cl <- C_1
  for (i in c(1:dim(C_1)[2])) {
    #print(i)
    if( length(unique(C_1[,i])) > 10 ){
      # 计算分位数
      C_qt_vector <- unique(quantile(C_1[, i], quantile_vec))

      # 检查去重后是否有足够的分位数
      if (length(C_qt_vector) >= n_divides - 1) {
        # 使用去重后的分位数作为 cut 的 breaks 参数
        C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf),
                          labels = 1:n_divides, include.lowest = TRUE)
        C_1_cl[, i] <- C_qt_label
      } else {
        C_qt_label <- cut(C_1[, i], breaks = c(-Inf, C_qt_vector, Inf),
                          labels = 1:(length(C_qt_vector)+1), include.lowest = TRUE)
        C_1_cl[, i] <- C_qt_label

        warning(paste("第", i, "列的分位数不足以进行划分,","现分为",(length(C_qt_vector)+1),"类"))
      }

    }
  }
  C_1_origin <- C_1
  C_1 <- C_1_cl
  C_1 <- as.matrix(C_1)
  C_1_cl <- as.matrix(C_1_cl)

  C_1 <- apply(C_1, 2, function(x) as.numeric(as.character(x)))
  Z_1 <- apply(Z_1, 2, function(x) as.numeric(as.character(x)))
  C_1_origin <- apply(C_1_origin, 2, function(x) as.numeric(as.character(x)))

  res_1 <- summary_compute(Z=Z_1,X=X_1,C=AS_1)

  beta_1 <- res_1$beta
  sd_1 <- res_1$sd
  p_1 <- res_1$p
  eaf_1 <- res_1$eaf
  R2_1 <- res_1$R2
  F_stats_1 <- res_1$F_stats

  data_1 <- list(Z=Z_1,C=C_1,X=X_1,beta=beta_1,sd=sd_1,p=p_1,eaf=eaf_1,R2=R2_1,F_stats=F_stats_1)

  #--------------------------------------------------
  res_0 <- summary_compute(Z=Z_0,X=Y_0,C=AS_0)
  beta_0 <- res_0$beta
  sd_0 <- res_0$sd
  p_0 <- res_0$p
  eaf_0 <- res_0$eaf
  R2_0 <- res_0$R2
  F_stats_0 <- res_0$F_stats
  #------------------
  C_0 <- as.data.frame(C_0)
  C_0_cl <- C_0
  for (i in c(1:dim(C_0)[2])) {
    if( length(unique(C_0[,i])) > 10 ){
      # 计算分位数
      C_qt_vector <- unique(quantile(C_0[, i], quantile_vec))

      # 检查去重后是否有足够的分位数
      if (length(C_qt_vector) >= n_divides - 1) {
        # 使用去重后的分位数作为 cut 的 breaks 参数
        C_qt_label <- cut(C_0[, i], breaks = c(-Inf, C_qt_vector, Inf),
                          labels = 1:n_divides, include.lowest = TRUE)
        C_0_cl[, i] <- C_qt_label
      } else {
        C_qt_label <- cut(C_0[, i], breaks = c(-Inf, C_qt_vector, Inf),
                          labels = 1:(length(C_qt_vector)+1), include.lowest = TRUE)
        C_0_cl[, i] <- C_qt_label

        warning(paste("第", i, "列的分位数不足以进行划分,","现分为",(length(C_qt_vector)+1),"类"))
      }

    }
  }


  C_0_origin <- C_0
  C_0 <- C_0_cl
  C_0 <- as.matrix(C_0)
  C_0_cl <- as.matrix(C_0_cl)


  C_0 <- apply(C_0, 2, function(x) as.numeric(as.character(x)))
  Z_0 <- apply(Z_0, 2, function(x) as.numeric(as.character(x)))
  C_0_origin <- apply(C_0_origin, 2, function(x) as.numeric(as.character(x)))

  data_0 <- list(Z=Z_0,C=C_0,Y=Y_0,beta=beta_0,sd=sd_0,p=p_0,eaf=eaf_0,R2=R2_0,F_stats=F_stats_0)

  # #-----------------------------------若不计算交互矩阵，就到这为止
  #=================================================================================================
  # 并行设置：检测并注册可用核心
  n_cores <- parallel::detectCores() - 2  # 保留一个核心
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # 创建空矩阵
  modifier_matrix <- matrix(0, nrow = p_Z, ncol = p_C)
  modifier_p <- matrix(NA, nrow = p_Z, ncol = p_C)

  # 外层对 k 并行
  results_list <- foreach(k = 1:p_Z, .packages = c("metafor")) %dopar% {
    modifier_row <- numeric(p_C)
    modifier_p_row <- numeric(p_C)

    for (j in 1:p_C) {
      data_chow_test <- data.frame(Z = Z_1[,k], X = X_1, C = C_1[,j])
      data_chow_test <- data_chow_test[order(data_chow_test$C), ]
      C_group_vec <- table(data_chow_test$C)

      end_indice <- cumsum(C_group_vec)
      start_indice <- c(1, head(end_indice, -1) + 1)

      beta_group_vec <- c()
      se_group_vec <- c()

      for (i in 1:length(C_group_vec)) {
        if (as.numeric(C_group_vec[i]) < 2) next
        data_chow_test_sub <- data_chow_test[start_indice[i]:end_indice[i], ]
        model_group_C <- lm(X ~ Z, data = data_chow_test_sub)
        res_model_group <- summary(model_group_C)$coefficients
        if (nrow(res_model_group) < 2) next
        beta_group_vec <- c(beta_group_vec, res_model_group[2,1])
        se_group_vec <- c(se_group_vec, res_model_group[2,2])
      }

      if (length(beta_group_vec) > 1) {
        res_meta_analysis <- metafor::rma.uni(yi = beta_group_vec, sei = se_group_vec, method = "REML")
        p_Q <- res_meta_analysis$QEp
        modifier_p_row[j] <- p_Q
        modifier_row[j] <- as.numeric(p_Q < p_threshold)
      } else {
        modifier_p_row[j] <- NA
        modifier_row[j] <- 0
      }
    }

    list(modifier_row = modifier_row, modifier_p_row = modifier_p_row)
  }

  # 汇总结果
  for (k in 1:p_Z) {
    modifier_matrix[k, ] <- results_list[[k]]$modifier_row
    modifier_p[k, ] <- results_list[[k]]$modifier_p_row
  }

  # 释放资源
  stopCluster(cl)


  C_0_save <- C_0
  C_1_save <- C_1
  C_0_origin_save <- C_0_origin
  C_1_origin_save <- C_1_origin
  #-----------------
  indice_Cnouse <- which(colSums(modifier_matrix) <1)
  #indice_Cnouse <- c(11,12,which(colSums(modifier_matrix) <1))
  C_0 <- C_0_save
  C_1 <- C_1_save
  C_0_origin <- C_0_origin_save
  C_1_origin <- C_1_origin_save

  #indice_Cuse <-
  setdiff(c(1:p_C),indice_Cnouse)
  C_1 <- C_1[,-indice_Cnouse]
  C_1_origin <- C_1_origin[,-indice_Cnouse]
  C_0 <- C_0[,-indice_Cnouse]
  C_0_origin <- C_0_origin[,-indice_Cnouse]

  #---------------------计算过程
  b_Y_list2 <- c()
  sd_Y_list2 <- c()
  # b_Y_list3 <- c()
  # sd_Y_list3 <- c()
  # b_Y_list4 <- c()
  # sd_Y_list4 <- c()

  # modifier_matrix_save <- modifier_matrix
  # modifier_matrix <- modifier_matrix[,setdiff(c(1:22),indice_Cnouse)]

  indice_modify <- which(rowSums(modifier_matrix) > 0)
  indice_notmodify <- setdiff(c(1:p_Z),which(rowSums(modifier_matrix) > 0))
  p_Zm <- length(indice_modify)
  p_Znm <- length(indice_notmodify)

  pval_i <- 0.05
  if(p_Znm>0){

    data3_0 <- data_0
    data3_1 <- data_1

    data3_0$Z <- data3_0$Z[,indice_notmodify,drop=FALSE]
    data3_0$beta <- data3_0$beta[indice_notmodify,drop=FALSE]
    data3_0$sd <- data3_0$sd[indice_notmodify,drop=FALSE]
    data3_0$eaf <- data3_0$eaf[indice_notmodify,drop=FALSE]
    data3_0$R2 <- data3_0$R2[indice_notmodify,drop=FALSE]
    data3_0$F_stats <- data3_0$F_stats[indice_notmodify,drop=FALSE]

    data3_1$Z <- data3_1$Z[,indice_notmodify,drop=FALSE]
    data3_1$beta <- data3_1$beta[indice_notmodify,drop=FALSE]
    data3_1$sd <- data3_1$sd[indice_notmodify,drop=FALSE]
    data3_1$eaf <- data3_1$eaf[indice_notmodify,drop=FALSE]
    data3_1$R2 <- data3_1$R2[indice_notmodify,drop=FALSE]
    data3_1$F_stats <- data3_1$F_stats[indice_notmodify,drop=FALSE]

    beta_hat_res <- Twosamle_package4(p_Z=length(indice_notmodify),p_C=p_C,data_0=data3_0,data_1=data3_1)
    beta_hat_degrade <- beta_hat_res$beta
    se_degrade <- beta_hat_res$se
    #----------------------------------
    if(p_Zm == 0){
      beta_hat_2 <- beta_hat_degrade
      beta_hat_3 <- beta_hat_degrade
      beta_hat_4 <- beta_hat_degrade
      b_Y_list2 <- c(b_Y_list2,beta_hat_degrade)
      b_Y_list3 <- c(b_Y_list3,beta_hat_degrade)
      b_Y_list4 <- c(b_Y_list4,beta_hat_degrade)
      sd_Y_list2 <- c(sd_Y_list2,se_degrade)
      sd_Y_list3 <- c(sd_Y_list3,se_degrade)
      sd_Y_list4 <- c(sd_Y_list4,se_degrade)
      sd_Y_list2_inv <- 1/sd_Y_list2
      sd_Y_list3_inv <- 1/sd_Y_list3
      sd_Y_list4_inv <- 1/sd_Y_list4
    }else{
      Z_0m <- Z_0[,indice_modify,drop=FALSE]
      Z_1m <- Z_1[,indice_modify,drop=FALSE]
      # if(sum(indice_modify) == 1){
      #   Z_0m <- matrix(Z_0m,ncol = 1)
      #   Z_1m <- matrix(Z_1m,ncol = 1)
      # }
      #===============================================================迁移估计量
      #Trans_OR
      X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0_origin,C_1=C_1_origin,X_1=X_1)$X_hat_matrix
      #Trans_IPW
      # X_hat_matrix_IPW <- IPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
      # #Trans_AIPW
      # X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0m,Z_1=Z_1m,C_0=C_0,C_1=C_1,X_1=X_1
      # )$X_hat_matrix
      #-----------------------------------------

      #---------------------------------------------------------------Trans_OR0         21
      if(p_Zm >= 1){
        for (i in 1:p_Zm) {
          model2 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0_origin + Z_0m[,i])  #有Z
          summary_model2 <- summary(model2)
          pval_model2 <- summary_model2$coefficients[,4]
          Z_pval <- pval_model2[length(pval_model2)]
          if(Z_pval < pval_i){
            b_Y_list2[i] <- model2$coefficients[2]
            sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
          }else{
            model2_1 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0_origin)  #无Z
            b_Y_list2[i] <- model2_1$coefficients[2]
            sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
          }
        }
      }
      b_Y_list2 <- c(b_Y_list2, beta_hat_degrade)
      sd_Y_list2 <- c(sd_Y_list2, se_degrade)

      sd_Y_list2_inv <- 1/sd_Y_list2
      beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
      # #---------------------------------------------------------------Trans_IPW0          31
      # if(p_Zm>=1){
      #
      #   for (i in 1:p_Zm) {
      #
      #     model3 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0 + Z_0m[,i])  #有Z
      #     summary_model3 <- summary(model3)
      #     pval_model3 <- summary_model3$coefficients[,4]
      #     Z_pval <- pval_model3[length(pval_model3)]
      #     print(as.numeric(Z_pval))
      #     if(Z_pval < pval_i){
      #       b_Y_list3[i] <- model3$coefficients[2]
      #       sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
      #     }else{
      #       model3_1 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0)  #无Z
      #       b_Y_list3[i] <- model3_1$coefficients[2]
      #       sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
      #     }
      #   }
      # }
      #
      # b_Y_list3 <- c(b_Y_list3, beta_hat_degrade)
      # sd_Y_list3 <- c(sd_Y_list3, se_degrade)
      #
      # sd_Y_list3_inv <- 1/sd_Y_list3
      # beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
      # #---------------------------------------------------------------Trans_AIPW0          41
      # if(p_Zm>=1){
      #
      #   for (i in 1:p_Zm) {
      #     model4 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0 + Z_0m[,i])  #有Z
      #     summary_model4 <- summary(model4)
      #     pval_model4 <- summary_model4$coefficients[,4]
      #     Z_pval <- pval_model4[length(pval_model4)]
      #     if(Z_pval<pval_i){
      #       b_Y_list4[i] <- model4$coefficients[2]
      #       sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
      #     }else{
      #       model4_1 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0)  #无Z
      #       b_Y_list4[i] <- model4_1$coefficients[2]
      #       sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
      #     }
      #   }
      # }
      # b_Y_list4 <- c(b_Y_list4, beta_hat_degrade)
      # sd_Y_list4 <- c(sd_Y_list4, se_degrade)
      #
      # sd_Y_list4_inv <- 1/sd_Y_list4
      # beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
    }
  }else{
    #===============================================================迁移估计量
    #Trans_OR
    X_hat_matrix_OR <- Trans_OLS2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
    # #Trans_IPW
    # X_hat_matrix_IPW <- IPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1)$X_hat_matrix
    # #Trans_AIPW
    # X_hat_matrix_AIPW <- AIPW2(Z_0=Z_0,Z_1=Z_1,C_0=C_0,C_1=C_1,X_1=X_1
    # )$X_hat_matrix
    #---------------------------------------------------------------Trans_OR0         21
    for (i in 1:p_Z) {
      model2 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0 + Z_0[,i])  #有Z
      summary_model2 <- summary(model2)
      pval_model2 <- summary_model2$coefficients[,4]
      Z_pval <- pval_model2[length(pval_model2)]
      print(as.numeric(Z_pval))
      #if(Z_pval < pval_i){
      if(Z_pval < pval_i){
        b_Y_list2[i] <- model2$coefficients[2]
        sd_Y_list2[i] <- summary(model2)$coefficients[,'Std. Error'][2]
      }else{
        model2_1 <- lm(Y_0~X_hat_matrix_OR[,i] + C_0)  #无Z
        b_Y_list2[i] <- model2_1$coefficients[2]
        sd_Y_list2[i] <- summary(model2_1)$coefficients[,'Std. Error'][2]
      }

    }
    # cor(X_hat_matrix_OR[,1],Z_0[,1])
    sd_Y_list2_inv <- 1/sd_Y_list2
    beta_hat_2 <- sum(sd_Y_list2_inv^2 * b_Y_list2)/sum(sd_Y_list2_inv^2)
    # #---------------------------------------------------------------Trans_IPW0          31
    # for (i in 1:p_Z) {
    #
    #   model3 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0 + Z_0[,i])  #有Z
    #   summary_model3 <- summary(model3)
    #   pval_model3 <- summary_model3$coefficients[,4]
    #   Z_pval <- pval_model3[length(pval_model3)]
    #   print(as.numeric(Z_pval))
    #   if(Z_pval < pval_i){
    #     b_Y_list3[i] <- model3$coefficients[2]
    #     sd_Y_list3[i] <- summary(model3)$coefficients[,'Std. Error'][2]
    #   }else{
    #     model3_1 <- lm(Y_0~X_hat_matrix_IPW[,i] + C_0)  #无Z
    #     b_Y_list3[i] <- model3_1$coefficients[2]
    #     sd_Y_list3[i] <- summary(model3_1)$coefficients[,'Std. Error'][2]
    #   }
    #
    #
    # }
    # sd_Y_list3_inv <- 1/sd_Y_list3
    # beta_hat_3 <- sum(sd_Y_list3_inv^2 * b_Y_list3)/sum(sd_Y_list3_inv^2)
    # #---------------------------------------------------------------Trans_AIPW0          41
    # for (i in 1:p_Z) {
    #   model4 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0 + Z_0[,i])  #有Z
    #   summary_model4 <- summary(model4)
    #   pval_model4 <- summary_model4$coefficients[,4]
    #   Z_pval <- pval_model4[length(pval_model4)]
    #   if(Z_pval < pval_i){
    #     b_Y_list4[i] <- model4$coefficients[2]
    #     sd_Y_list4[i] <- summary(model4)$coefficients[,'Std. Error'][2]
    #   }else{
    #     model4_1 <- lm(Y_0~X_hat_matrix_AIPW[,i] + C_0)  #无Z
    #     b_Y_list4[i] <- model4_1$coefficients[2]
    #     sd_Y_list4[i] <- summary(model4_1)$coefficients[,'Std. Error'][2]
    #   }
    # }
    # sd_Y_list4_inv <- 1/sd_Y_list4
    # beta_hat_4 <- sum(sd_Y_list4_inv^2 * b_Y_list4)/sum(sd_Y_list4_inv^2)
  }

  #----------------------------结果保存
  beta_hat_vec <- c(beta_hat_2)
  sd_OR <- sqrt(1/sum(sd_Y_list2_inv^2))
  # sd_IPW <- sqrt(1/sum(sd_Y_list3_inv^2))
  # sd_AIPW <- sqrt(1/sum(sd_Y_list4_inv^2))

  p_power_OR <- 2 * (1 - pnorm(abs(beta_hat_2 / sd_OR)))
  # p_power_IPW <- 2 * (1 - pnorm(abs(beta_hat_3 / sd_IPW)))
  # p_power_AIPW <- 2 * (1 - pnorm(abs(beta_hat_4 / sd_AIPW)))
  # p_power_mm <- c(p_power_OR,p_power_IPW,p_power_AIPW)
  p_power_mm <- c(p_power_OR)
  p_vec <- c(p_power_mm)

  #sd_hat_vec <- c(beta_hat_1_list$se,sd_OR,sd_IPW,sd_AIPW)
  sd_hat_vec <- c(beta_hat_1_list$se,sd_OR)
  p_vec <- c(p_power_mm)
  # uci_vec <- beta_hat_vec + 1.95 * sd_hat_vec
  # lci_vec <- beta_hat_vec - 1.95 * sd_hat_vec

  beta_hat_vec_t <- round(beta_hat_vec,3)
  sd_hat_vec_t <-round(sd_hat_vec,4)
  # uci_vec_t <- round(uci_vec,4)
  # lci_vec_t <- round(lci_vec,4)

  return_list <- list(beta=beta_hat_vec_t, se = sd_hat_vec_t, p = p_vec)

  return(return_list)
}
