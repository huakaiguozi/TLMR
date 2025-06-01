#-------------------------------------------------------------------------------
Total_F <- function(n, k, R2) {
  return(((n - k - 1) / k) * (sum(R2) / (1 - sum(R2))))
}

summary_compute <- function(Z,X,C){
  N <- dim(Z)[1]
  p_Z <- dim(Z)[2]
  indicator_bin <- if(n_distinct(X) < 5){TRUE}else{FALSE}  # TRUE:分类；FALSE:连续

  Z <- as.matrix(Z)
  X <- as.matrix(X)
  C <- as.matrix(C)
  #------------------
  if(!indicator_bin){
    beta_summary <- rep(NA,p_Z)
    beta_summary_p <- rep(NA,p_Z)

    beta_summary_sd <- rep(NA,p_Z)
    R_2 <- rep(NA,p_Z)
    F_stats <- rep(NA,p_Z)

    beta_summary_eaf <- colSums(Z)/ (2 * N)  #也是正确的
    #------------------
    for (z in 1:p_Z){
      data_ivw_x <- as.data.frame(cbind(X,Z[,z],C[,1],C[,2]))
      names(data_ivw_x) <- c('x','z')
      Beta_res <- feols(x ~ z, data=data_ivw_x)
      Beta_res_coef <- coef(Beta_res)  ##报错1
      if ('(Intercept)' %in% names(Beta_res_coef)){
        beta_hat <- Beta_res_coef[2]
      }else{
        beta_hat <- Beta_res_coef[1]
      }
      #------------------summary数据
      beta_summary[z] <- beta_hat
      beta_summary_sd[z] <- Beta_res$coeftable[,2][2]
      beta_summary_p[z] <- Beta_res$coeftable$`Pr(>|t|)`[2]
      # f_1 <- sum(Z[,z]==1)/N  #也可以
      # f_2 <- sum(Z[,z]==2)/N
      # beta_summary_eaf[z] <- f_1/2 + f_2
    }

    for(i in 1:p_Z){
      var_z <- 2 * beta_summary_eaf[i] *  (1-beta_summary_eaf[i])
      term1 <- beta_summary[i]^2 * var_z
      R_2[i] <- term1 / ( term1 + (beta_summary_sd[i]^2 * N * var_z ) )
    }

    for(j in 1:p_Z){
      F_stats[j] <- ( ( N-2 ) * R_2[j] ) / ( 1-R_2[j]  )
    }
  }else{
    beta_summary <- rep(NA,p_Z)
    beta_summary_p <- rep(NA,p_Z)

    beta_summary_sd <- rep(NA,p_Z)
    R_2 <- rep(NA,p_Z)
    F_stats <- rep(NA,p_Z)

    beta_summary_eaf <- colSums(Z)/ (2 * N)  #也是正确的
    #------------------
    for (z in 1:p_Z){
      # data_ivw_x <- as.data.frame(cbind(X,Z[,z]))
      # names(data_ivw_x) <- c('x','z')
      # Beta_res <- feols(x ~ z, data=data_ivw_x)

      model_glm <- glm(X~Z[,z]+C, family = binomial(link = "log"))
      summary_glm <- summary(model_glm)
      Beta_res_coef <- summary_glm$coefficients[,1]
      Beta_res_sd <- summary_glm$coefficients[,2]
      Beta_res_p <- summary_glm$coefficients[,4]
      #Beta_res_coef <- coef(Beta_res)  ##报错1
      if ('(Intercept)' %in% names(Beta_res_coef)){
        beta_hat <- as.numeric(Beta_res_coef[2])
        sd_hat <- as.numeric(Beta_res_sd[2])
        p_hat <- as.numeric(Beta_res_p[2])
      }else{
        beta_hat <- as.numeric(Beta_res_coef[1])
        sd_hat <- as.numeric(Beta_res_sd[1])
        p_hat <- as.numeric(Beta_res_p[1])
      }
      #------------------summary数据
      beta_summary[z] <- beta_hat
      beta_summary_sd[z] <- sd_hat
      beta_summary_p[z] <- p_hat
      # f_1 <- sum(Z[,z]==1)/N  #也可以
      # f_2 <- sum(Z[,z]==2)/N
      # beta_summary_eaf[z] <- f_1/2 + f_2
    }

    for(i in 1:p_Z){
      var_z <- 2 * beta_summary_eaf[i] *  (1-beta_summary_eaf[i])
      term1 <- beta_summary[i]^2 * var_z
      R_2[i] <- term1 / ( term1 + (beta_summary_sd[i]^2 * N * var_z ) )
    }

    for(j in 1:p_Z){
      F_stats[j] <- ( ( N-2 ) * R_2[j] ) / ( 1-R_2[j]  )
    }
  }
  res_list <- list(beta = beta_summary, sd=beta_summary_sd, p=beta_summary_p,
                   eaf=beta_summary_eaf, R2=R_2, F_stats=F_stats )
  return(res_list)
}



logistic_function <- function(x) {
  return(1 / (1 + exp(-x)))
}

Trans_OLS2 <- function(Z_1,C_1,X_1,Z_0,C_0){   ####Trans-OR
  n_0 <- dim(Z_0)[1]
  n_1 <- dim(X_1)[1]
  p_Z <- dim(Z_1)[2]
  p_C <- dim(C_1)[2]

  #----------------------------------------回归
  X_hat3_slices <- matrix(nrow = n_0,ncol = dim(Z_0)[2])
  for(j in 1:dim(Z_0)[2]){
    interaction_term_1 <- matrix(rep(NA,n_1*p_C),nrow = n_1)   #改成分p_z维度：z是标量
    for (i in 1:p_C) {
      interaction_term_1[, i] <- Z_1[,j] * C_1[,i]
    }
    O_x1 <- cbind(1,Z_1[,j],C_1,interaction_term_1)
    Theta_X_1_with_itc <- solve(t(O_x1) %*% O_x1) %*% t(O_x1) %*% X_1
    itc_1 <- Theta_X_1_with_itc[1]
    Theta_X_1 <- Theta_X_1_with_itc[-1]
    #---提取参数估计
    alpha_hat <- Theta_X_1[1]
    eta_1_hat <- Theta_X_1[c((2):(1+p_C))]
    pi_hat <- matrix(Theta_X_1[c((p_C+2):length(Theta_X_1))],nrow = 1)
    #--------
    ### 预测0中的X_hat
    diag_result1 <- rowSums((Z_0[,j] %*% pi_hat) * C_0)
    X_hat3 <- Z_0[,j] * alpha_hat + C_0 %*% eta_1_hat + diag_result1 + itc_1
    #X_hat3 <- Z_0[,j] * alpha_hat + C_0 %*% eta_1_hat + diag(Z_0[,j] %*% pi_hat %*% t(C_0))
    X_hat3_slices[,j] <- X_hat3
  }


  res_list <- list(X_hat_matrix = X_hat3_slices)
  return(res_list)
}

IPW2 <- function(Z_0,Z_1,C_0,C_1,X_1){  #频数估计概率
  n_0 <- dim(Z_0)[1]
  n_1 <- dim(Z_1)[1]
  n <- n_0 + n_1
  P_hat1 <- n_1/(n_0+n_1)
  #-----------------------
  X_hat_slices <- matrix(nrow = n_0,ncol = dim(Z_0)[2])
  for (i_z in 1:dim(Z_0)[2]) {
    Z_slice_0 <- as.matrix(Z_0[,i_z])
    Z_slice_1 <- as.matrix(Z_1[,i_z])

    V_0 <- cbind(Z_slice_0,C_0)
    V_1 <- cbind(Z_slice_1,C_1)
    V_0_string <- matrix(apply(V_0, 1, paste, collapse = ","))
    V_1_string <- matrix(apply(V_1, 1, paste, collapse = ","))

    V_0_unique <- unique(V_0)
    V_0_unique_string <- matrix(apply(V_0_unique, 1, paste, collapse = ","))

    p_V0_unique <- apply(V_0_unique_string, 1, function(row,data){sum(row==data)},data=V_1_string) / n_1


    p_V <- rep(0,n_0)  # 就是P(z,c|S=1)
    X_matrix <- matrix(rep(0,ceiling(max(p_V0_unique)* n_1)*n_0),nrow = n_0 )

    for (i_v in 1:dim(V_0_unique)[1]) {
      p_V[V_0_unique_string[i_v,] == V_0_string] <- p_V0_unique[i_v]

      X_slice <- X_1[V_0_unique_string[i_v,] == V_1_string]
      n_X_slices <- sum(V_0_unique_string[i_v,] == V_0_string)
      if(n_X_slices > 1){
        # 重复行向量n次;转置矩阵，使得每一行都是重复的行向量
        X_slices <- t(replicate(n_X_slices, X_slice))
        X_matrix[V_0_unique_string[i_v,] == V_0_string,][,1:length(X_slice)]  <- X_slices
      }else if(n_X_slices == 1){
        X_matrix[V_0_unique_string[i_v,] == V_0_string,][1:length(X_slice)]  <- X_slice
      }
    }
    p_hat1 <- n_1/(n_0+n_1)
    weight <- 1/(p_V * p_hat1)
    X_hat_slice <- weight * rowSums(X_matrix) / (n_0+n_1)
    X_hat_slices[,i_z] <- X_hat_slice
  }

  #-----------------------

  res_list <- list(X_hat_matrix = X_hat_slices)
  return(res_list)
}

AIPW2 <- function(Z_0,Z_1,C_0,C_1,X_1){  #频数估计概率
  n_0 <- dim(Z_0)[1]
  n_1 <- dim(Z_1)[1]
  n <- n_0 + n_1
  P_hat1 <- n_1/(n_0+n_1)
  p_C <- ncol(C_0)
  #----------------------------------------回归
  X_gap_slices <- matrix(nrow = n_1,ncol = dim(Z_0)[2])
  X_hat3_slices <- matrix(nrow = n_0,ncol = dim(Z_0)[2])
  for(j in 1:dim(Z_0)[2]){
    interaction_term_1 <- matrix(rep(NA,n_1*p_C),nrow = n_1)   #改成分p_z维度：z是标量
    for (i in 1:p_C) {
      interaction_term_1[, i] <- Z_1[,j] * C_1[,i]
    }
    O_x1 <- cbind(1,Z_1[,j],C_1,interaction_term_1)
    Theta_X_1_with_itc <- solve(t(O_x1) %*% O_x1) %*% t(O_x1) %*% X_1
    itc_1 <- Theta_X_1_with_itc[1]
    Theta_X_1 <- Theta_X_1_with_itc[-1]
    #---提取参数估计
    alpha_hat <- Theta_X_1[1]
    eta_1_hat <- Theta_X_1[c((2):(1+p_C))]
    pi_hat <- matrix(Theta_X_1[c((p_C+2):length(Theta_X_1))],nrow = 1)
    #--------
    diag_result1 <- rowSums((Z_1[,j] %*% pi_hat) * C_1)
    X_hat2 <- Z_1[,j] * alpha_hat + C_1 %*% eta_1_hat + diag_result1 + itc_1
    X_gap_slices[,j] <- X_1 - X_hat2
    ### 预测0中的X_hat
    diag_result2 <- rowSums((Z_0[,j] %*% pi_hat) * C_0)
    X_hat3 <- Z_0[,j] * alpha_hat + C_0 %*% eta_1_hat + diag_result2 + itc_1
    X_hat3_slices[,j] <- X_hat3
  }

  #-----------------------
  X_hat_slices <- matrix(nrow = n_0,ncol = dim(Z_0)[2])
  for (i_z in 1:dim(Z_0)[2]) {
    Z_slice_0 <- as.matrix(Z_0[,i_z])
    Z_slice_1 <- as.matrix(Z_1[,i_z])

    V_0 <- cbind(Z_slice_0,C_0)
    V_1 <- cbind(Z_slice_1,C_1)
    V_0_string <- matrix(apply(V_0, 1, paste, collapse = ","))
    V_1_string <- matrix(apply(V_1, 1, paste, collapse = ","))

    V_0_unique <- unique(V_0)
    V_0_unique_string <- matrix(apply(V_0_unique, 1, paste, collapse = ","))

    p_V0_unique <- apply(V_0_unique_string, 1, function(row,data){sum(row==data)},data=V_1_string) / n_1   #每一个V0出现在V1中的的概率：P(z,c|S=1)


    p_V <- rep(0,n_0)  # 空的：n_0个P(z,c|S=1)向量
    X_matrix <- matrix(rep(0,ceiling(max(p_V0_unique)* n_1)*n_0),nrow = n_0 )  #对每一个P(z,c|S=1)中的z、c，出现在dataset1中的X的个数：
    #每一行为所有dataset1中对应z、c的X

    for (i_v in 1:dim(V_0_unique)[1]) {  #对每一个z、c
      p_V[V_0_unique_string[i_v,] == V_0_string] <- p_V0_unique[i_v]  #填补：n_0个P(z,c|S=1)向量

      X_slice <- X_gap_slices[,i_z][V_0_unique_string[i_v,] == V_1_string]  #取出对应(z、c)的所有X向量
      n_X_slices <- sum(V_0_unique_string[i_v,] == V_0_string)  #dataset0中有多少个(z、c)
      if(n_X_slices > 1){
        # 重复行向量n次;转置矩阵，使得每一行都是重复的行向量
        X_slices <- t(replicate(n_X_slices, X_slice))
        X_matrix[V_0_unique_string[i_v,] == V_0_string,][,1:length(X_slice)]  <- X_slices
      }else if(n_X_slices == 1){
        X_matrix[V_0_unique_string[i_v,] == V_0_string,][1:length(X_slice)]  <- X_slice
      }
    }
    p_hat1 <- n_1/(n_0+n_1)
    weight <- 1/(p_V * p_hat1)
    X_hat_slice <- weight * rowSums(X_matrix) / (n_0+n_1)
    X_hat_slices[,i_z] <- X_hat_slice
  }
  X_hat4_slices <- X_hat_slices + X_hat3_slices
  #-----------------------
  res_list <- list(X_hat_matrix = X_hat4_slices)
  return(res_list)
}

# 迭代留一法主函数（针对元分析）
iterative_loo_meta <- function(data, yi_col, vi_col, max_iter = 5, i2_threshold = 25) {
  # 参数说明:
  # data: 输入数据集（需包含效应量(yi)和方差(vi)列）
  # yi_col: 效应量列名（字符类型）
  # vi_col: 方差列名（字符类型）
  # max_iter: 最大迭代次数（默认5）
  # i2_threshold: I²异质性容忍阈值（默认25%，即低于25%时停止）

  current_data <- data              # 当前数据集（动态更新）
  removed_studies <- character(0)   # 记录被剔除的研究
  results_history <- list()         # 记录每次迭代的元分析结果

  if(dim(data)[1]<5){
    max_iter <- floor(dim(data)[1] / 2)
  }

  for (iter in 1:max_iter) {
    # ---------- 步骤1: 执行留一交叉验证 ----------
    loo_results <- leave1out(
      rma(
        yi = get(yi_col),
        vi = get(vi_col),
        data = current_data,
        method = "REML"  # 使用限制性最大似然估计
      )
    )

    # ---------- 步骤2: 计算每个研究剔除后的I²值 ----------
    i2_values <- loo_results$I2
    max_i2_idx <- which.max(i2_values)  # 找到剔除后I²最大的研究

    # ---------- 步骤3: 记录结果并检查终止条件 ----------
    current_model <- rma(
      yi = get(yi_col),
      vi = get(vi_col),
      data = current_data,
      method = "REML"
    )
    results_history[[iter]] <- list(
      model = current_model,
      removed = rownames(current_data)[max_i2_idx]
    )

    # 如果当前I²已低于阈值，提前终止
    if (current_model$I2 < i2_threshold) break

    # ---------- 步骤4: 剔除对异质性贡献最大的研究 ----------
    removed_studies <- c(removed_studies, rownames(current_data)[max_i2_idx])
    current_data <- current_data[-max_i2_idx, ]
  }

  # ---------- 返回结果 ----------
  list(
    final_model = current_model,
    removed_studies = removed_studies,
    history = results_history
  )
}
