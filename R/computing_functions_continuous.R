#--------------------------------------- 计算方法函数
#===============================
#--------------------------------Two-sample MR
MR_lasso<-function(betaYG,betaXG,sebetaYG){
  betaYGw = betaYG/sebetaYG # dividing the association estimates by  sebetaYG is equivalent
  betaXGw = betaXG/sebetaYG # to weighting by sebetaYG^-2
  pleio = diag(rep(1, length(betaXG)))
  l1grid = c(seq(from=0.1, to=5, by=0.1), seq(from=5.2, to=10, by=0.2))
  # values of lambda for grid search
  l1grid_rse = NULL; l1grid_length = NULL; l1grid_beta = NULL;
  l1grid_se = NULL
  for (i in 1:length(l1grid)) {
    l1grid_which = which(attributes(penalized(betaYGw, pleio,
                                              betaXGw, lambda1=l1grid[i], trace=FALSE))$penalized==0)
    l1grid_rse[i] = summary(lm(betaYG[l1grid_which]~betaXG[l1grid_which]-1, weights=sebetaYG[l1grid_which]^-2))$sigma
    l1grid_length[i] = length(l1grid_which)
    l1grid_beta[i] = lm(betaYG[l1grid_which]~betaXG[l1grid_which]-1,
                        weights=sebetaYG[l1grid_which]^-2)$coef[1]
    l1grid_se[i] = summary(lm(betaYG[l1grid_which]~betaXG[l1grid_which]-1, weights=sebetaYG[l1grid_which]^-2))$coef[1,2]/
      min(summary(lm(betaYG[l1grid_which]~betaXG[l1grid_which]-1,weights=sebetaYG[l1grid_which]^-2))$sigma, 1)

  }
  l1which_hetero = c(which(l1grid_rse[1:(length(l1grid)-1)]>1& diff(
    l1grid_rse)>qchisq(0.95, df=1)/l1grid_length[2:length(l1grid)])
    , length(l1grid))[1]
  # heterogeneity criterion for choosing lambda

  l1hetero_beta = l1grid_beta[l1which_hetero]

  l1hetero_se = l1grid_se[l1which_hetero]
  list(ThetaEstimate=l1hetero_beta, ThetaSE=l1hetero_se )

}

Twosamle_package3 <- function(p_Z,p_C,data_0,data_1){
  n_0 <- dim(data_0$Z)[1]
  n_1 <- dim(data_1$Z)[1]
  #beta se eaf
  exp_data <- data.frame(id.exposure = rep(1,p_Z),exposure=rep("x",p_Z),SNP=c(1:p_Z),beta.exposure=data_1$beta,se.exposure=data_1$sd,
                         effect_allele.exposure=rep('A',p_Z),other_allele.exposure=rep('T',p_Z),
                         eaf.exposure=data_1$eaf,samplesize.exposure=rep(n_1,p_Z))
  out_data <- data.frame(id.outcome = rep(2,p_Z),outcome=rep("y",p_Z),SNP=c(1:p_Z),beta.outcome=data_0$beta,se.outcome=data_0$sd,
                         effect_allele.outcome=rep('A',p_Z),other_allele.outcome=rep('T',p_Z),
                         eaf.outcome=data_0$eaf,samplesize.outcome=rep(n_0,p_Z))

  data_merge <- harmonise_data(exposure_dat = exp_data,outcome_dat = out_data,action=2)
  #------------
  data_merge2 <- MendelianRandomization::mr_input(bx=data_merge$beta.exposure,bxse = data_merge$se.exposure
                                                  ,by=data_merge$beta.outcome,byse = data_merge$se.outcome)
  # num_ivs <- length(data_merge2@snps)
  # if(num_ivs<3){
  #   res <- MendelianRandomization::mr_ivw(data_merge2)
  #   beta <- res@Estimate
  #   se <- res@StdError
  # }else{
  #   res <- MendelianRandomization::mr_egger(data_merge2)
  #   beta <- res@Estimate
  #   se <- res@StdError.Est
  # }
  res <- MendelianRandomization::mr_cML(data_merge2,n=n_0)
  beta <- res@Estimate
  se <- res@StdError

  return_list <- list(beta=beta,se=se)
  return(return_list)
}

Twosamle_package2 <- function(p_Z,p_C,data_0,data_1){
  n_0 <- dim(data_0$Z)[1]
  n_1 <- dim(data_1$Z)[1]
  exp_data <- data.frame(id.exposure = rep(1,p_Z),exposure=rep("x",p_Z),SNP=c(1:p_Z),beta.exposure=data_1$beta,se.exposure=data_1$sd,
                         effect_allele.exposure=rep('A',p_Z),other_allele.exposure=rep('T',p_Z),
                         eaf.exposure=data_1$eaf,samplesize.exposure=rep(n_1,p_Z))
  out_data <- data.frame(id.outcome = rep(2,p_Z),outcome=rep("y",p_Z),SNP=c(1:p_Z),beta.outcome=data_0$beta,se.outcome=data_0$sd,
                         effect_allele.outcome=rep('A',p_Z),other_allele.outcome=rep('T',p_Z),
                         eaf.outcome=data_0$eaf,samplesize.outcome=rep(n_0,p_Z))

  data_merge <- harmonise_data(exposure_dat = exp_data,outcome_dat = out_data,action=2)
  #------------
  data_merge2 <- MendelianRandomization::mr_input(bx=data_merge$beta.exposure,bxse = data_merge$se.exposure
                                                  ,by=data_merge$beta.outcome,byse = data_merge$se.outcome)

  res_ivw <- MendelianRandomization::mr_ivw(data_merge2)
  res <- list(Q=res_ivw@Heter.Stat,F=res_ivw@Fstat)
  return(res)
}

Twosamle_package <- function(p_Z,p_C,data_0,data_1,beta){
  n_0 <- dim(data_0$Z)[1]
  n_1 <- dim(data_1$Z)[1]

  methods_select <- c('mr_ivw','mr_egger_regression','mr_weighted_median','mr_weighted_mode')

  exp_data <- data.frame(id.exposure = rep(1,p_Z),exposure=rep("x",p_Z),SNP=c(1:p_Z),beta.exposure=data_1$beta,se.exposure=data_1$sd,
                         effect_allele.exposure=rep('A',p_Z),other_allele.exposure=rep('T',p_Z),
                         eaf.exposure=data_1$eaf,samplesize.exposure=rep(n_1,p_Z))
  out_data <- data.frame(id.outcome = rep(2,p_Z),outcome=rep("y",p_Z),SNP=c(1:p_Z),beta.outcome=data_0$beta,se.outcome=data_0$sd,
                         effect_allele.outcome=rep('A',p_Z),other_allele.outcome=rep('T',p_Z),
                         eaf.outcome=data_0$eaf,samplesize.outcome=rep(n_0,p_Z))

  data_merge <- harmonise_data(exposure_dat = exp_data,outcome_dat = out_data,action=2)

  res <- mr(data_merge,method_list = methods_select)
  #------
  # 初始化 b_4methods 和 se_4methods
  b_4methods <- rep(99999, 4)
  se_4methods <- rep(99999, 4)

  # 对应关系
  method_mapping <- c("Inverse variance weighted" = "mr_ivw",
                      "MR Egger" = "mr_egger_regression",
                      "Weighted median" = "mr_weighted_median",
                      "Weighted mode" = "mr_weighted_mode")

  # 将 res 中的 method 映射到 methods_select 的索引
  for (i in 1:nrow(res)) {
    method_name <- method_mapping[res$method[i]]
    idx <- match(method_name, methods_select) # 找到 method 在 methods_select 中的索引
    if (!is.na(idx)) {
      b_4methods[idx] <- res$b[i]
      se_4methods[idx] <- res$se[i]
    }
  }
  #------
  # res_raps <-mr.raps.shrinkage(b_exp = data_merge$beta.exposure, b_out = data_merge$beta.outcome
  #                              , se_exp = data_merge$se.exposure, se_out = data_merge$se.outcome)
  res_raps <- mr.raps(data_merge,over.dispersion=FALSE)
  beta_raps <- res_raps$beta.hat
  se_raps <- res_raps$beta.se
  p_raps <- 2 * (1 - pnorm(abs(beta_raps/se_raps)))

  #------------
  data_merge2 <- MendelianRandomization::mr_input(bx=data_merge$beta.exposure,bxse = data_merge$se.exposure
                                                  ,by=data_merge$beta.outcome,byse = data_merge$se.outcome)

  #res_ivw <- MendelianRandomization::mr_ivw(data_merge2)

  res_conmix <- MendelianRandomization::mr_conmix(data_merge2)
  #res_robust <- MendelianRandomization::mr_ivw(data_merge2,"random", robust = TRUE)
  #----
  indicator_robust <- FALSE
  res_robust <- tryCatch({
    MendelianRandomization::mr_ivw(data_merge2, "random", robust = TRUE)
  }, error = function(e) {
    indicator_robust <<- TRUE
  })
  #----
  indicator_lasso <- FALSE
  res_lasso <- tryCatch({
    MR_lasso(data_merge$beta.outcome,data_merge$beta.exposure,data_merge$se.outcome)
  }, error = function(e) {
    indicator_lasso <<- TRUE
  })

  #------------
  if(!indicator_robust & !indicator_lasso){
    res_list <- list(beta=c(b_4methods, res_raps$beta.hat, res_conmix@Estimate, res_robust@Estimate, res_lasso$ThetaEstimate)
                     ,se=c(se_4methods, res_raps$beta.se, res_conmix@Psi, res_robust@StdError, res_lasso$ThetaSE))
  }else if(indicator_robust & !indicator_lasso){
    res_list <- list(beta=c(b_4methods, res_raps$beta.hat, res_conmix@Estimate, 99999, res_lasso$ThetaEstimate)
                     ,se=c(se_4methods, res_raps$beta.se, res_conmix@Psi, 99999, res_lasso$ThetaSE))
  }else if(!indicator_robust & indicator_lasso){
    res_list <- list(beta=c(b_4methods, res_raps$beta.hat, res_conmix@Estimate, res_robust@Estimate, 99999)
                     ,se=c(se_4methods, res_raps$beta.se, res_conmix@Psi, res_robust@Estimate, 99999))
  }else{
    res_list <- list(beta=c(b_4methods, res_raps$beta.hat, res_conmix@Estimate, 99999, 99999)
                     ,se=c(se_4methods, res_raps$beta.se, res_conmix@Psi, 99999, 99999))
  }
  p_vec <- 2 * (1 - pnorm(abs(res_list$beta/res_list$se)))
  p_cover_vec <- 2 * (1 - pnorm(abs((res_list$beta-beta)/res_list$se)))

  p_vec[which(res_list$beta == 99999)] <- 99999
  p_cover_vec[which(res_list$beta == 99999)] <- 99999

  res_list[['p']] <- p_vec
  res_list[['p_cover']] <- p_cover_vec
  return(res_list)
  #return(res)
}

#--------------------------------Trans-OLS
Trans_OLS <- function(Z,C,X){

  n <- dim(X)[1]
  p_Z <- dim(Z)[2]
  p_C <- dim(C)[2]

  interaction_term <- matrix(rep(NA,n*p_Z*p_C),nrow = n)
  for (i in 1:p_C) {
    interaction_term[,(((i-1) *p_Z+1):(i*p_Z))] <- Z * C[,i]
  }
  O_x <- cbind(Z,C,interaction_term)
  theta_hat <- solve(t(O_x) %*% O_x) %*% t(O_x) %*% X

  Q_theta <- t(X) %*% X -t(X) %*% (O_x %*% theta_hat)
  sigma_2 <- Q_theta / (n-p_Z-p_C-p_Z*p_C)
  var_theta <- sigma_2[1,1] * solve(t(O_x) %*% O_x)

  res_list <- list(coef = theta_hat, var = var_theta)
  return(res_list)
}


#--------------------------------Trans-Lasso
Trans_Lasso_X <- function(Z,C,X){

  n <- dim(X)[1]
  p_Z <- dim(Z)[2]
  p_C <- dim(C)[2]

  interaction_term <- matrix(rep(NA,n*p_Z*p_C),nrow = n) # 纠错2
  for (i in 1:p_C) {
    interaction_term[,(((i-1) *p_Z+1):(i*p_Z))] <- Z * C[,i]
  }
  O_x <- cbind(Z,C,interaction_term)


  cv_model <- cv.glmnet(O_x, X, alpha = 1,nfolds = 5)

  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  indicator <- 0
  while (indicator ==0 ) {


    best_model <- glmnet(O_x, X, alpha = 1, lambda = best_lambda,intercept = F)

    theta_matrix <- as.matrix(coef(best_model))
    theta_hat <- as.vector(theta_matrix)[2:dim(theta_matrix)[1]]

    alpha_hat <- theta_hat[c(1:p_Z)]
    eta_1_hat <- theta_hat[c((p_Z+1):(p_Z+p_C))]
    pi_hat <- matrix(theta_hat[c((p_Z+p_C+1):length(theta_hat))],nrow = p_Z)

    if( any(pi_hat != 0) & !any(alpha_hat==0) ){
      indicator =1
    }else{
      best_lambda <- best_lambda * 0.7
    }

  }
  var_theta <- 0

  res_list <- list(coef = theta_hat, alpha_hat=alpha_hat, eta_1_hat=eta_1_hat, pi_hat=pi_hat, var = var_theta)
  return(res_list)
}
#--------------------------------两步GMM
GMM_2S <- function(Z_0,C_0,Z_1,C_1,X,Y,theta){  #修改7-整个复制过去
  n_0 <- dim(Z_0)[1]
  n_1 <- dim(Z_1)[1]
  p_Z <- dim(Z_0)[2]
  p_C <- dim(C_0)[2]
  #-----------------
  W <- diag(nrow = p_Z + p_C + p_Z * p_C)
  #---
  O_x0 <- cbind(Z_0,C_0)
  for (i in 1:p_C) {
    O_x0 <- cbind(O_x0,C_0[,i] * Z_0)
  }
  #---
  O_x1 <- cbind(Z_1,C_1)
  for (i in 1:p_C) {
    O_x1 <- cbind(O_x1,C_1[,i] * Z_1)
  }
  #---
  O_y1 <- cbind(Z_1,C_1,X)  #含X的数据集

  #--------
  M_xx0 <- t(O_x0) %*% O_x0 / n_0
  M_xx1 <- t(O_x1) %*% O_x1 / n_1
  M_xy1 <- t(O_x1) %*% O_y1 / n_1
  M_xY0 <- t(O_x0) %*% Y / n_0
  #-----------------
  Omega_hat <- (solve(M_xx0) / n_0) %*% (t(O_x0) %*% diag(diag(Y %*% t(Y))) %*% O_x0 -
                                           t(O_x0) %*% Y %*% t(t(O_x0) %*% Y) / n_0) %*% t(solve(M_xx0)) +
    (solve(M_xx1) / n_1) %*% (t(O_x1) %*% diag(diag(O_y1 %*% theta %*% t(O_y1 %*% theta))) %*% O_x1 -
                                t(O_x1) %*% O_y1 %*% theta %*% t(t(O_x1) %*% (O_y1 %*% theta)) / n_1) %*% t(solve(M_xx1))
  W <- solve(Omega_hat)

  theta <- solve(t(M_xy1) %*% solve(M_xx1) %*% W %*% solve(M_xx1) %*% M_xy1) %*% (
    t(M_xy1) %*% solve(M_xx1) %*% W %*% solve(M_xx0) %*% M_xY0
  )
  res_list <- list(coef = theta)
  return(res_list)
}
#------------------------------------------------- CUE计算
#--------------辅助函数
Objective_function <- function(x,data){
  O_x0 <- data$O_x0
  O_x1<- data$O_x1
  O_y1 <- data$O_y1
  M_xx0 <- data$M_xx0
  M_xx1 <- data$M_xx1
  M_xy1 <- data$M_xy1
  M_xY0 <- data$M_xY0
  Y <- data$Y
  n_0 <- data$n_0
  n_1 <- data$n_1
  theta <- x
  #------------------
  result <- t( (solve(M_xx0) ) %*% M_xY0 - (solve(M_xx1)) %*% M_xy1 %*% theta )%*%
    solve((solve(M_xx0) / n_0) %*% (t(O_x0) %*% diag(diag(Y %*% t(Y))) %*% O_x0 -
                                      t(O_x0) %*% Y %*% t(t(O_x0) %*% Y) / n_0) %*% t(solve(M_xx0)) +
            (solve(M_xx1) / n_1) %*% (t(O_x1) %*% diag(diag(O_y1 %*% theta %*% t(O_y1 %*% theta))) %*% O_x1 -
                                        t(O_x1) %*% O_y1 %*% theta %*% t(t(O_x1) %*% (O_y1 %*% theta)) / n_1) %*%
            t(solve(M_xx1))) %*%((solve(M_xx0) ) %*% M_xY0 - (solve(M_xx1)) %*% M_xy1 %*% theta)
  return(result)
}

#--------------辅助函数(结束)
GMM_CUE <- function(Z_0,C_0,Z_1,C_1,X,Y,theta){
  n_0 <- dim(Z_0)[1]
  n_1 <- dim(Z_1)[1]
  p_Z <- dim(Z_0)[2]
  p_C <- dim(C_0)[2]
  #-----------------
  O_x0 <- cbind(Z_0,C_0)
  for (i in 1:p_C) {
    O_x0 <- cbind(O_x0,C_0[,i] * Z_0)
  }
  #---
  O_x1 <- cbind(Z_1,C_1)
  for (i in 1:p_C) {
    O_x1 <- cbind(O_x1,C_1[,i] * Z_1)
  }
  #---
  O_y1 <- cbind(Z_1,C_1,X)  #含X的数据集

  #--------
  M_xx0 <- t(O_x0) %*% O_x0 / n_0
  M_xx1 <- t(O_x1) %*% O_x1 / n_1
  M_xy1 <- t(O_x1) %*% O_y1 / n_1
  M_xY0 <- t(O_x0) %*% Y / n_0
  #-----------------
  data_CUE <-list(O_x0 = O_x0, O_x1=O_x1, O_y1=O_y1, M_xx0=M_xx0, M_xx1=M_xx1,
                  M_xy1=M_xy1, M_xY0=M_xY0, Y = Y, n_0 = n_0, n_1 = n_1)

  opt <- nloptr(
    x0 = theta,
    #x0 = rep(2,(p_Z+p_C+1)),
    eval_f = function(x) Objective_function(x, data_CUE),
    #opts = list("algorithm" = "NLOPT_LN_NELDERMEAD","xtol_rel"=1.0e-4,"print_level" = 1, "maxeval"=250)
    #opts = list("algorithm" = "NLOPT_LN_NELDERMEAD","xtol_rel"=1.0e-4,"print_level" = 1, "maxeval"=250)
    opts = list("algorithm" = "NLOPT_LN_SBPLX","xtol_rel"=1.0e-4,"print_level" = 1, "maxeval"=100)
  )
  theta_hat <- opt$solution
  res_list <- list(coef = theta_hat)
  # beta_hat <- theta_hat[(p_Z+p_C+1)]
  # res_list <- list(coef = beta_hat)
  return(res_list)
}

IPW <- function(Z_0,Z_1,C_0,C_1,X_1){     ##逻辑回归估计概率
  # #------------------- 调试 得删
  # Z_0_save <- Z_0
  # Z_1_save <- Z_1
  # iii <- 4
  # Z_0 <- as.matrix(Z_0[,iii],ncol=1)
  # Z_1 <- as.matrix(Z_1[,iii],ncol=1)
  # #------------------- 调试 得删

  n_0 <- dim(Z_0)[1]
  n_1 <- dim(Z_1)[1]
  n <- n_0 + n_1
  P_hat1 <- n_1/(n_0+n_1)

  #-----------------------
  C <- rbind(C_0,C_1)
  C_0_cl <- C_0
  C_1_cl <- C_1
  for (i in c(1:dim(C_0)[2])) {
    if( length(unique(C_0[,i])) > 10 ){
      quant_vec <- c(0.3333,0.6667)
      q_vector <- quantile(C[,i],quant_vec)
      label <- cut( C_0[,i],c(-Inf,q_vector,Inf),labels=c(1:(length(quant_vec)+1)) )
      C_0_cl[,i] <- label
      label <- cut( C_1[,i],c(-Inf,q_vector,Inf),labels=c(1:(length(quant_vec)+1)) )
      C_1_cl[,i] <- label
    }
  }
  #-----------------------
  V_0 <- cbind(Z_0,C_0)
  V_1 <- cbind(Z_1,C_1)
  V_0_cl <- cbind(Z_0,C_0_cl)
  V_1_cl <- cbind(Z_1,C_1_cl)
  V_bind <- rbind(V_0_cl,V_1_cl)

  Sample <- c(rep(0,dim(V_0)[1]),rep(1,dim(V_1)[1]))  #样本选择
  logi_data <- data.frame(
    y <- Sample,
    x <- V_bind
  )
  logi_model <- glm(y ~ x, data = logi_data, family = "binomial")
  Pr_Slv <- predict(logi_model,newdata = data.frame(V_bind), type = 'response')[1:n_0]
  #------------------------------------------- V的概率估计：P(V)
  pte_estimand <-  unique(V_0_cl)
  Pr_estimand <- rep(0,dim(pte_estimand)[1])
  for (i in 1:dim(pte_estimand)[1]) {
    fre_pte <- sum(apply(matrix(pte_estimand[i,] == t(V_bind),ncol = dim(V_bind)[2],byrow = TRUE),1,prod))
    #if (fre_pte > 0){
    Pr_estimand[i] <- fre_pte / (n_0 + n_1)
    #}
  }

  Pr_V_vec <- rep(0,dim(V_0)[1])
  for (i in 1:length(Pr_estimand)) {
    fill_V0_index <- apply(V_0_cl, 1, function(row,data) {all(row == data)},data=pte_estimand[i,])
    Pr_V_vec[fill_V0_index] <- Pr_estimand[i]
  }
  Pr_VS <- as.vector(Pr_Slv * Pr_V_vec * P_hat1)   #
  #--------------------------------------------------------
  pte_V0_unique <- unique(V_0_cl)
  V_0_unique_index1 <- outer(1:dim(pte_V0_unique)[1], 1:nrow(V_1_cl), Vectorize(function(i, j) all(pte_V0_unique[i,] == V_1_cl[j,])))
  V_0_unique_index2 <- outer(1:nrow(V_0_cl),1:dim(pte_V0_unique)[1], Vectorize(function(i, j) all(V_0_cl[i,] == pte_V0_unique[j,])))
  X_matrix <- t(apply(V_0_unique_index2, 1, function(r,data){
    data[r,]
  },data=V_0_unique_index1))

  X_sum <- apply(X_matrix, 1, function(r,data){
    sum(data[r])
  },data=X_1)

  #--------------------------------------------------------
  pte_observed <- unique(V_1_cl)

  pte_dismatch <- as.matrix(dplyr::anti_join(unique(V_0_cl) %>%  as.data.frame(), pte_observed %>%as.data.frame()))

  if (dim(pte_dismatch)[1]>0){
    index_dismatch <- FALSE
    for (i in 1:dim(pte_dismatch)[1]) {
      index_dismatch <- index_dismatch | apply(pte_estimand, 1, function(r,data){
        all(r ==data)
      },data=pte_dismatch[i,])
    }
    Pr_dismatch <- Pr_estimand[index_dismatch]
    num_vector <- ceiling( Pr_dismatch * (n_0+n_1)*(n_1/(n_0+n_1)) )
    k_neighbors7 <- max(num_vector)
    D_dist7 <- knn(V_1_cl,pte_dismatch,cl=c(1:dim(V_1_cl)[1]) , k = k_neighbors7, algorithm = "kd_tree")
    D_index7 <- attr(D_dist7, "nn.index")
    D_dist7 <- attr(D_dist7, 'nn.dist')


    dismatch_matrix <- matrix(rep(FALSE,n_1*dim(D_index7)[1]),nrow = dim(D_index7)[1])
    X_sum_dismatch <- rep(0,dim(D_index7)[1])
    for (i in 1:dim(D_index7)[1]) {
      dismatch_matrix[i,][D_index7[i,][1:num_vector[i]]] <- TRUE
      X_sum_dismatch[i] <- sum(X_1[dismatch_matrix[i,]])
    }


    V0_dismatch_index <- which(X_sum == 0)
    V0_dismatch_matrix <- V_0_cl[V0_dismatch_index,]
    if(dim(pte_dismatch)[1]==1){
      V0_dismatch_matrix <- t(as.matrix(V0_dismatch_matrix))
    }

    #X_dismatch_matrix <- matrix(rep(0,length(V0_dismatch_index)*n_1),ncol=n_1)
    X_sum_dismatch_all <- rep(0,length(V0_dismatch_index))
    for (i in 1:dim(D_index7)[1]) {
      X_sum_dismatch_all[apply(V0_dismatch_matrix,1,function(r,data){
        all(r==data)
      },data=pte_dismatch[i,])] <- X_sum_dismatch[i]
    }
    X_sum[V0_dismatch_index] <- X_sum_dismatch_all
  }
  #-----------------------
  weight <- 1/Pr_VS /(n_0+n_1)
  #knn_weight_matrix <- (1/D_dist3)/rowSums(1/D_dist3)   #再试试权重改了对不对
  X_hat <- weight * X_sum

  res_list <- list(X_hat = X_hat)
  return(res_list)
}

IPW3 <- function(Z_0,Z_1,C_0,C_1,X_1,data_type){  #频数估计概率
  n_0 <- dim(Z_0)[1]
  n_1 <- dim(Z_1)[1]
  n <- n_0 + n_1
  p_C <- dim(C_0)[2]
  P_hat1 <- n_1/(n_0+n_1)
  #-----------------------
  C <- rbind(C_0,C_1)
  #-----------------------
  X_hat_slices <- matrix(nrow = n_0,ncol = dim(Z_0)[2])
  p_z_vec <- rep(NA,n_0)
  p_c_vec <- rep(NA,n_0)
  for (i_z in 1:dim(Z_0)[2]) {
    Z_slice_0 <- as.matrix(Z_0[,i_z])
    Z_slice_1 <- as.matrix(Z_1[,i_z])

    p_z1 <- c(sum(Z_slice_1 == 0)/n_1
              ,sum(Z_slice_1 == 1)/n_1
              ,sum(Z_slice_1 == 2)/n_1)

    p_z_vec[Z_slice_0 == 0] <- p_z1[1]
    p_z_vec[Z_slice_0 == 1] <- p_z1[2]
    p_z_vec[Z_slice_0 == 2] <- p_z1[3]

    if(data_type == 1){

      kde_result <- ks::kde(x=C_1,density = TRUE)
      p_c_vec <- predict(kde_result,x = C_0)

    }else if(data_type == 2){   #拉普拉斯平滑


      C_list <- vector("list", p_C)
      for (i_c in 1:p_C) {

        C_list[[i_c]] <- unique(C[, i_c])
      }

      C_combinations <- expand.grid(C_list)
      C_elements <- apply(C_combinations, 1, function(row) paste(row, collapse = ","))
      K_elements <- length(C_elements)
      #p_C_elements <- rep(NA,K_elements)

      C_0_string <- matrix(apply(C_0, 1, paste, collapse = ","))
      C_1_string <- matrix(apply(C_1, 1, paste, collapse = ","))

      for (i_g in 1:K_elements) {
        p_element <- (sum(C_1_string == C_elements[i_g])+1)/(n_1+K_elements)
        #p_C_elements[i_g] <- p_element
        p_c_vec[C_0_string == C_elements[i_g]] <- p_element
      }

    }else if(data_type == 3){
      #----先分好哪些是分类变量
      indice_cat <- rep(NA,p_C) # TRUE为分类，FALSE为连续
      values_thresh <- 7
      for (i_c2 in 1:p_C) {
        values_num <- length(unique(C[,i_c2]))
        if(values_num < values_thresh){
          indice_cat[i_c2] <- TRUE
        }else{
          indice_cat[i_c2] <- FALSE
        }
      }
      #----调整分类变量到factor
      C_0_df <- as.data.frame(C_0)
      for (i_c3 in 1:p_C) {
        if (indice_cat[i_c3]) {
          C_0_df[, i_c3] <- as.factor(C_0_df[, i_c3])
        } else {
          C_0_df[, i_c3] <- C_0[, i_c3]
        }
      }

      C_1_df <- as.data.frame(C_1)
      for (i_c3 in 1:p_C) {
        if (indice_cat[i_c3]) {
          C_1_df[, i_c3] <- as.factor(C_1_df[, i_c3])
        } else {
          C_1_df[, i_c3] <- C_1[, i_c3]
        }
      }
      #----
      bw_C <- np::npudensbw(dat = C_1_df, ftol = 1e-4, nmulti=1)
      density_estimate_C <- np::npudens(bw_C)
      p_c_vec <- predict(density_estimate_C, new_data = C_0_df)

    }

    p_zcs <- p_z_vec * p_c_vec * P_hat1

    # 概率计算完了. 但是后续，I(·)可能找不到取值，后续代码在重叠性不足时，需要近似（后续代码还没写完，后面的代码需要删除）


    V_0 <- cbind(Z_slice_0,C_0_cl)
    V_1 <- cbind(Z_slice_1,C_1_cl)
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




AIPW <- function(Z_0,Z_1,C_0,C_1,X_1){     ##逻辑回归估计概率

  n_0 <- dim(Z_0)[1]
  n_1 <- dim(Z_1)[1]
  n <- n_0 + n_1
  P_hat1 <- n_1/(n_0+n_1)

  #-----------------------
  Z <- rbind(Z_0,Z_1)
  C <- rbind(C_0,C_1)
  C_0_cl <- C_0
  C_1_cl <- C_1
  for (i in c(1:dim(C_0)[2])) {
    if( length(unique(C_0[,i])) > 10 ){
      quant_vec <- c(0.3333,0.6667)
      q_vector <- quantile(C[,i],quant_vec)
      label <- cut( C_0[,i],c(-Inf,q_vector,Inf),labels=c(1:(length(quant_vec)+1)) )
      C_0_cl[,i] <- label
      label <- cut( C_1[,i],c(-Inf,q_vector,Inf),labels=c(1:(length(quant_vec)+1)) )
      C_1_cl[,i] <- label
    }
  }
  #----------------------------------------回归
  interaction_term_1 <- matrix(rep(NA,n_1*p_Z*p_C),nrow = n_1) # 纠错2
  for (i in 1:p_C) {
    interaction_term_1[,(((i-1) *p_Z+1):(i*p_Z))] <- Z_1 * C_1[,i]
  }
  O_x1 <- cbind(Z_1,C_1,interaction_term_1)
  Theta_X_1 <- solve(t(O_x1) %*% O_x1) %*% t(O_x1) %*% X_1
  #---提取参数估计
  alpha_hat <- Theta_X_1[c(1:p_Z)]
  eta_1_hat <- Theta_X_1[c((p_Z+1):(p_Z+p_C))]
  pi_hat <- matrix(Theta_X_1[c((p_Z+p_C+1):length(Theta_X_1))],nrow = p_Z)
  #--------
  X_hat2 <- Z_1 %*% alpha_hat + C_1 %*% eta_1_hat + diag(Z_1 %*% pi_hat %*% t(C_1))
  X_gap <- X_1 - X_hat2
  #X_hat3 <- X_hat1-X_hat2 * (1/Pr_VS) + X_hat2

  #------------------------------
  #-----------------------
  V_0 <- cbind(Z_0,C_0)
  V_1 <- cbind(Z_1,C_1)
  V_0_cl <- cbind(Z_0,C_0_cl)
  V_1_cl <- cbind(Z_1,C_1_cl)
  V_bind <- rbind(V_0_cl,V_1_cl)

  Sample <- c(rep(0,dim(V_0)[1]),rep(1,dim(V_1)[1]))  #样本选择
  logi_data <- data.frame(
    y <- Sample,
    x <- V_bind
  )
  logi_model <- glm(y ~ x, data = logi_data, family = "binomial")
  Pr_Slv <- predict(logi_model,newdata = data.frame(V_bind), type = 'response')[1:n_0]
  #------------------------------------------- V的概率估计：P(V)
  pte_estimand <-  unique(V_0_cl)
  Pr_estimand <- rep(0,dim(pte_estimand)[1])
  for (i in 1:dim(pte_estimand)[1]) {
    fre_pte <- sum(apply(matrix(pte_estimand[i,] == t(V_bind),ncol = dim(V_bind)[2],byrow = TRUE),1,prod))
    #if (fre_pte > 0){
    Pr_estimand[i] <- fre_pte / (n_0 + n_1)
    #}
  }

  Pr_V_vec <- rep(0,dim(V_0)[1])
  for (i in 1:length(Pr_estimand)) {
    fill_V0_index <- apply(V_0_cl, 1, function(row,data) {all(row == data)},data=pte_estimand[i,])
    Pr_V_vec[fill_V0_index] <- Pr_estimand[i]
  }
  Pr_VS <- as.vector(Pr_Slv * Pr_V_vec * P_hat1)   #
  #--------------------------------------------------------
  pte_V0_unique <- unique(V_0_cl)
  # pte_V0_unique中的点（行），在V1哪些是一样的
  V_0_unique_index1 <- outer(1:dim(pte_V0_unique)[1], 1:nrow(V_1_cl), Vectorize(function(i, j) all(pte_V0_unique[i,] == V_1_cl[j,])))
  V_0_unique_index2 <- outer(1:nrow(V_0_cl),1:dim(pte_V0_unique)[1], Vectorize(function(i, j) all(V_0_cl[i,] == pte_V0_unique[j,])))
  X_matrix <- t(apply(V_0_unique_index2, 1, function(r,data){
    data[r,]
  },data=V_0_unique_index1))

  X_sum <- apply(X_matrix, 1, function(r,data){
    sum(data[r])
  },data=X_gap)

  #--------------------------------------------------------
  pte_observed <- unique(V_1_cl)

  pte_dismatch <- as.matrix(dplyr::anti_join(unique(V_0_cl) %>%  as.data.frame(), pte_observed %>% as.data.frame()))

  #-----------------------
  weight <- 1/Pr_VS /(n_0+n_1)
  X_hat1 <- weight * X_sum
  #----------------------------------------AIPW
  X_hat3 <- Z_0 %*% alpha_hat + C_0 %*% eta_1_hat + diag(Z_0 %*% pi_hat %*% t(C_0))

  X_hat4 <- X_hat1 + X_hat3

  #------------------------------
  res_list <- list(X_hat = X_hat4)
  return(res_list)
}

# 示例：自定义岭回归的 Bootstrap 函数
ridge_reg <- function(data, indices, lambda) {
  sample_data <- data[indices, ]
  best_model <- glmnet::glmnet(sample_data[,-1], sample_data[,1], alpha = 0, lambda = lambda,intercept = F)
  return(best_model$beta[1])
}
