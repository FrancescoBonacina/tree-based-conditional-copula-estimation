################################################################################
#
#   Tree & Clayton copula estimation.
#
#   - COVARIATES: stepwise combination of x1,x2
#   - OBS: estimation on true observations (u1,u2) and on pseudo-observations
#          (v1,v2), (w1,w2)
#   - Likelihood: iTau estimation
#   - Pruning: OOB logL maximization
#
################################################################################

rm(list=ls())

library(rpart)
library(rpart.plot)
library(rgl)     # 3D plot library
library(MASS)
library(copula)
library(sigmoid)
library(caret)

### set directory for output files
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
wd <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
outdir <- paste0('30_out')
dir.create(file.path(wd, outdir), showWarnings = FALSE)
outdir_path <- paste0(wd, '/', outdir)

### simulation parameters   =========================================================================================
copule <- 'clayton'             # frank, gumbel, clayton
ndataset <- 1                 # n. of datasets
N <- 1000                       # n. obs
# simulating copula samples
taus <- c(0.3,0.5,0.7,0.9)      # taus (~thetas), copula parameters 
cov_combination <- 'stepwise'   # stepwise, gentle_sigmoid, steep_sigmoid
th_x1 <- 0.4                    # threshold positions along x1 and x2
th_x2 <- 0.75
# splitting and evaluation functions
verbose <- FALSE                # to know when these functions are called by rpart
scaling <- 0.001                # scaling factor for logL estimation
min_leaf_size <- 10             # estimate logL only if the node contains at least min_leaf_size observations
# CV procedure
cp <- 0.01                      # complexity parameter of the tree
kfolds <- 4                     # n. of folds
times_fold <- 1                 # m. of repetitions of k-fold --> total n. of folds will be kfold*times_fold
max_nleaves <- 8                # do not compute likelihood for trees with >max_nleaves

### define splitting methods based on copulas
source('./trees&copulas_init-split-eval_functions.R')
clayton_method_uu <- list(eval=eval_clayton_uu, split=split_clayton_uu, init=init)
clayton_method_vv <- list(eval=eval_clayton_vv, split=split_clayton_vv, init=init)
clayton_method_ww <- list(eval=eval_clayton_ww, split=split_clayton_ww, init=init)

### run simulations   =================================================================================================
mat_output <- matrix(NA, nrow=ndataset, ncol=21)
mat_columns <- c('logL_bench_uu','MSE_tau_bench_uu','MSE_C_bench_uu','logL_uu','MSE_tau_uu','MSE_C_uu','nsplits_uu',
                 'logL_bench_vv','MSE_tau_bench_vv','MSE_C_bench_vv','logL_vv','MSE_tau_vv','MSE_C_vv','nsplits_vv',
                 'logL_bench_ww','MSE_tau_bench_ww','MSE_C_bench_ww','logL_ww','MSE_tau_ww','MSE_C_ww','nsplits_ww')
start_time <- Sys.time()
for (D in 1:ndataset){
  print(D)
  #D <- 1
  
  ###   create data   =================================================================================================
  # covariates
  set.seed(D)
  x1 <- runif(N)
  x2 <- runif(N)
  y_idx <- 1:N
  
  # simulate copula samples: tau(x1,x2), theta(tau), copula obs (u1,u2) and marginals (y1,y2) 
  tau_x <- (taus[2]-taus[1])*(1+sign(x1-th_x1))*1/2 + (taus[3]-taus[1])*(1+sign(x2-th_x2))*1/2 + taus[1] 
  theta_x <- iTau(claytonCopula(dim=2), tau_x)
  uu <- t(sapply(theta_x, function(i){rCopula(1, claytonCopula(param=i, dim=2))}))
  y1 <- qnorm(uu[,1], mean = 1 + 0.2*x1 + 0.05*x2)
  y2 <- qnorm(uu[,2], mean = -1 - 0.1*x1 + 0.2*x2)

  # simulate 2 types of pseudo-obs: (v1,v2) and (w1,w2)
  # pseudo-obs (v1,v2): we assume to know that margins are Gaussian distributions
  #                     with std=1 and mean depending on the covariates. 
  #                     Thus, we use a linear model to estimate the mean of the distrib.
  lm1 <- lm(y1~x1+x2, data = data.frame(y1=y1, y2=y2, x=x1, x2=x2))
  a1_hat <- lm1$coefficients[1]
  b1_hat <- lm1$coefficients[2]
  c1_hat <- lm1$coefficients[3]
  lm2 <- lm(y2~x1+x2, data = data.frame(y1=y1, y2=y2, x=x1, x2=x2))
  a2_hat <- lm2$coefficients[1]
  b2_hat <- lm2$coefficients[2]
  c2_hat <- lm2$coefficients[3]
  v1 <- pnorm(y1, mean = a1_hat + b1_hat*x1 + c1_hat*x2)
  v2 <- pnorm(y2, mean = a2_hat + b2_hat*x1 + c2_hat*x2)
  # pseudo-obs (w1,w2): we assume that we don't know the marginal distributions.
  #                     We use a kernel method conditional to the covariates:
  #                     in practice, we compute an ECDF, by giving more weight to
  #                     the obs. which are close in terms of covariates.
  # vec_h <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.3, 0.33, 0.36, 0.4, 0.43, 0.46, 0.5, 1)
  # err <- c()
  # vec_logL <- c()
  # for (h in vec_h){
  #   w1 <- apply(matrix(c(x1,x2,y1), ncol=3), 1, function(obs){a<-obs[1]; b<-obs[2]; c<-obs[3]; sum( dnorm((x1-a)/h)*dnorm((x2-b)/h)*(y1<=c) ) / sum( dnorm((x1-a)/h)*dnorm((x2-b)/h) ) })
  #   w2 <- apply(matrix(c(x1,x2,y2), ncol=3), 1, function(obs){a<-obs[1]; b<-obs[2]; c<-obs[3]; sum( dnorm((x1-a)/h)*dnorm((x2-b)/h)*(y2<=c) ) / sum(dnorm( (x1-a)/h)*dnorm((x2-b)/h) ) })
  #   w1 <- w1-0.0001
  #   w2 <- w2-0.0001
  #   err <- c(err, mean(sqrt((w1-uu[,1])^2 + (w2-uu[,2])^2)))
  #   yy <- matrix(c(w1,w2), ncol=2)
  #   tau <- cor(yy, method="kendall")[1,2]
  #   theta <- iTau(claytonCopula(dim=2), tau)
  #   vec_logL <- c(vec_logL, loglikCopula(param=theta, u=yy, copula=claytonCopula(dim=2)))
  # }
  # vec_h[which.min(err)]
  # vec_h[which.max(vec_logL)]
  h <- 0.4   # bandwidth of the kernel estimator
  w1 <- apply(matrix(c(x1,x2,y1), ncol=3), 1, function(obs){a<-obs[1]; b<-obs[2]; c<-obs[3]; sum( dnorm((x1-a)/h)*dnorm((x2-b)/h)*(y1<=c) ) / sum( dnorm((x1-a)/h)*dnorm((x2-b)/h) ) }) - 0.0001
  w2 <- apply(matrix(c(x1,x2,y2), ncol=3), 1, function(obs){a<-obs[1]; b<-obs[2]; c<-obs[3]; sum( dnorm((x1-a)/h)*dnorm((x2-b)/h)*(y2<=c) ) / sum(dnorm( (x1-a)/h)*dnorm((x2-b)/h) ) }) - 0.0001
  
  df_cop <- data.frame(y_idx=y_idx,
                       x1=x1,
                       x2=x2,
                       tau_x=tau_x,
                       theta_x=theta_x,
                       u1=uu[,1],
                       u2=uu[,2],
                       v1=v1,
                       v2=v2,
                       w1=w1,
                       w2=w2)
  #str(df_cop)
  #nrow(df_cop[df_cop$tau_x==taus[1],]) 
  #nrow(df_cop[df_cop$tau_x==my_taus[2],]) 
  #nrow(df_cop[df_cop$tau_x==my_taus[3],])
  #nrow(df_cop[df_cop$tau_x>my_taus[3],])
  
  ### CV procedure for tree pruning: we need to decide the number of leaves   =========================================
  # define folds for cross-validation
  folds <- createMultiFolds(1:N, k=kfolds, times=times_fold)   
  CV_results_uu <- matrix(NA, nrow=length(folds), ncol=max_nleaves)
  CV_results_vv <- matrix(NA, nrow=length(folds), ncol=max_nleaves)
  CV_results_ww <- matrix(NA, nrow=length(folds), ncol=max_nleaves)
  for (f in 1:length(folds)){
    #print (f)
    #f <- 1
    ### test and train
    fold <- folds[[f]]
    train <- df_cop[df_cop$y_idx %in% fold, ]
    test <- df_cop[!(df_cop$y_idx %in% fold),]
    
    ### train tree
    # true observations
    tree_complete_uu <- rpart(data = train,
                              formula = y_idx ~ x1 + x2,
                              method = clayton_method_uu,
                              control = rpart.control(cp=cp),
                              xval = 0)
    # pseudo-observations
    tree_complete_vv <- rpart(data = train,
                              formula = y_idx ~ x1 + x2,
                              method = clayton_method_vv,
                              control = rpart.control(cp=cp),
                              xval = 0)
    tree_complete_ww <- rpart(data = train,
                              formula = y_idx ~ x1 + x2,
                              method = clayton_method_ww,
                              control = rpart.control(cp=cp),
                              xval = 0)

    ### assess goodness of trees, from 0 to max number of splits
    cps_uu <- tree_complete_uu$cptable[,c('CP')]
    nsplits_uu <- tree_complete_uu$cptable[,c('nsplit')]
    cps_vv <- tree_complete_vv$cptable[,c('CP')]
    nsplits_vv <- tree_complete_vv$cptable[,c('nsplit')]
    cps_ww <- tree_complete_ww$cptable[,c('CP')]
    nsplits_ww <- tree_complete_ww$cptable[,c('nsplit')]
    for (c in 1:max_nleaves){
      
      # true observations
      if (!c>length(cps_uu)){
        mycp <- cps_uu[[c]]
        mynsplit <- nsplits_uu[[c]]
        if (!(is.na(mycp) | mynsplit+1>max_nleaves)){   # if number of leaves is defined and not > of max n. of leaves
          # pruned tree
          mytree <- prune(tree_complete_uu, cp=mycp)
          # predictions
          taus_hat <- predict(mytree, test)
          df_pred_ <- data.frame(thetas_hat=iTau(copula=claytonCopula(dim=2), tau=taus_hat))
          df_pred_$y_idx <- as.numeric(rownames(df_pred_))
          df_pred <- merge(df_pred_, df_cop[,c('y_idx','u1','u2')], by='y_idx', all=FALSE)
          # compute logL
          logL <- 0
          for (theta in unique(df_pred$thetas_hat)){
            df <- df_pred[df_pred$thetas_hat==theta,]
            logL <- logL + loglikCopula(param=theta, u=as.matrix(df[, c('u1', 'u2')], ncol=2), copula=claytonCopula(dim=2))
          }
          # store nsplist and logL
          CV_results_uu[f, mynsplit+1] <- scaling*logL
        }
      }
      
      # pseudo-obs
      if (!c>length(cps_vv)){
        mycp <- cps_vv[[c]]
        mynsplit <- nsplits_vv[[c]]
        if (!(is.na(mycp) | mynsplit+1>max_nleaves)){   # if number of leaves is defined and not > of max n. of leaves
          # pruned tree
          mytree <- prune(tree_complete_vv, cp=mycp)
          # predictions
          taus_hat <- predict(mytree, test)
          df_pred_ <- data.frame(thetas_hat=iTau(copula=claytonCopula(dim=2), tau=taus_hat))
          df_pred_$y_idx <- as.numeric(rownames(df_pred_))
          df_pred <- merge(df_pred_, df_cop[,c('y_idx','u1','u2')], by='y_idx', all=FALSE)
          # compute logL
          logL <- 0
          for (theta in unique(df_pred$thetas_hat)){
            df <- df_pred[df_pred$thetas_hat==theta,]
            logL <- logL + loglikCopula(param=theta, u=as.matrix(df[, c('u1', 'u2')], ncol=2), copula=claytonCopula(dim=2))
          }
          # store nsplist and logL
          CV_results_vv[f, mynsplit+1] <- scaling*logL
        }
      }
      
      if (!c>length(cps_ww)){
        mycp <- cps_ww[[c]]
        mynsplit <- nsplits_ww[[c]]
        if (!(is.na(mycp) | mynsplit+1>max_nleaves)){   # if number of leaves is defined and not > of max n. of leaves
          # pruned tree
          mytree <- prune(tree_complete_ww, cp=mycp)
          # predictions
          taus_hat <- predict(mytree, test)
          df_pred_ <- data.frame(thetas_hat=iTau(copula=claytonCopula(dim=2), tau=taus_hat))
          df_pred_$y_idx <- as.numeric(rownames(df_pred_))
          df_pred <- merge(df_pred_, df_cop[,c('y_idx','u1','u2')], by='y_idx', all=FALSE)
          # compute logL
          logL <- 0
          for (theta in unique(df_pred$thetas_hat)){
            df <- df_pred[df_pred$thetas_hat==theta,]
            logL <- logL + loglikCopula(param=theta, u=as.matrix(df[, c('u1', 'u2')], ncol=2), copula=claytonCopula(dim=2))
          }
          # store nsplist and logL
          CV_results_ww[f, mynsplit+1] <- scaling*logL
        }
      }
      
    }
  }
  
  ### Find optimal number of leaves
  #   I look at values of average logL for increasing number of leaves.
  #   At the beginning logL increases, than it sort of stabilizes, oscillating around a maximum value.
  #   I follow a forward selection: I stop at the first step at which logL starts decreasing.
  CV_means_uu <- apply(CV_results_uu, 2, function(x) mean(x, na.rm = TRUE))
  increasing_logL_uu <- CV_means_uu[2:max_nleaves] > CV_means_uu[1:(max_nleaves-1)]
  best_nleaves_uu <- if (0 %in% increasing_logL_uu) {which.min(increasing_logL_uu)} else {max_nleaves}
  CV_means_vv <- apply(CV_results_vv, 2, function(x) mean(x, na.rm = TRUE))
  increasing_logL_vv <- CV_means_vv[2:max_nleaves] > CV_means_vv[1:(max_nleaves-1)]
  best_nleaves_vv <- if (0 %in% increasing_logL_vv) {which.min(increasing_logL_vv)} else {max_nleaves}
  CV_means_ww <- apply(CV_results_ww, 2, function(x) mean(x, na.rm = TRUE))
  increasing_logL_ww <- CV_means_ww[2:max_nleaves] > CV_means_ww[1:(max_nleaves-1)]
  best_nleaves_ww <- if (0 %in% increasing_logL_ww) {which.min(increasing_logL_ww)} else {max_nleaves}
  
  ### boxplots of CV
  # boxplot to compare logL for true-obs vs. pseudo-obs
  #out_file_CV <- sprintf('%s/CV_stepwise_clayton_N%d_4x1folds.pdf', outdir_path, N)
  #pdf(file=out_file_CV, width=12, height=4)
  #par(mfrow=c(1,3))
  #nleaves <- as.character(rep(1:ncol(CV_results_uu), nrow(CV_results_uu)*2))
  #obs_category <- c(rep('true-obs', ncol(CV_results_uu)*nrow(CV_results_uu)), rep('pseudo-obs', ncol(CV_results_vv)*nrow(CV_results_vv)))
  #vec_logL <- c(c(t(CV_results_uu)), c(t(CV_results_vv)))
  #df_boxplot <- data.frame(nleaves=nleaves, obs_category=obs_category, logL=vec_logL)
  #ggplot(df_boxplot, aes(x=nleaves, y=logL, fill=obs_category)) + geom_boxplot()
  
  #CV_sd <- apply(CV_results, 2, function(x) sd(x, na.rm = TRUE))
  #plot(1:length(CV_means), CV_means)
  #arrows(x0=1:length(CV_means), y0=CV_means-CV_sd, x1=1:length(CV_means), y1=CV_means+CV_sd, code=3, angle=90, length=0.1)
  #points(best_nleaves, CV_means[best_nleaves], pch='X', cex=2, col='red')
  #CV_sd_pseudo <- apply(CV_results_pseudo, 2, function(x) sd(x, na.rm = TRUE))
  #plot(1:length(CV_means_pseudo), CV_means_pseudo)
  #arrows(x0=1:length(CV_means_pseudo), y0=CV_means_pseudo-CV_sd_pseudo, x1=1:length(CV_means_pseudo), y1=CV_means_pseudo+CV_sd_pseudo, code=3, angle=90, length=0.1)
  #points(best_nleaves_pseudo, CV_means_pseudo[best_nleaves_pseudo], pch='X', cex=2, col='red')
  #dev.off()
  
  ### fit tree on all obs/pseudo-obs   ================================================================================
  tree_max_uu <- rpart(data = df_cop,
                       formula = y_idx ~ x1 + x2,
                       method = clayton_method_uu,
                       control = rpart.control(cp=cp),
                       xval = 0)
  tree_max_vv <- rpart(data = df_cop,
                       formula = y_idx ~ x1 + x2,
                       method = clayton_method_vv,
                       control = rpart.control(cp=cp),
                       xval = 0)
  tree_max_ww <- rpart(data = df_cop,
                       formula = y_idx ~ x1 + x2,
                       method = clayton_method_ww,
                       control = rpart.control(cp=cp),
                       xval = 0)

  # params for pruning
  cptab_uu <- data.frame(tree_max_uu$cptable)
  best_nsplit_uu <- max(cptab_uu[cptab_uu$nsplit<best_nleaves_uu, c('nsplit')])
  best_cp_uu <- cptab_uu[cptab_uu$nsplit==best_nsplit_uu, c('CP')]
  cptab_vv <- data.frame(tree_max_vv$cptable)
  best_nsplit_vv <- max(cptab_vv[cptab_vv$nsplit<best_nleaves_vv, c('nsplit')])
  best_cp_vv <- cptab_vv[cptab_vv$nsplit==best_nsplit_vv, c('CP')]
  cptab_ww <- data.frame(tree_max_ww$cptable)
  best_nsplit_ww <- max(cptab_ww[cptab_ww$nsplit<best_nleaves_ww, c('nsplit')])
  best_cp_ww <- cptab_ww[cptab_ww$nsplit==best_nsplit_ww, c('CP')]
  
  tree_uu <- prune(tree_max_uu, cp=best_cp_uu)
  tree_vv <- prune(tree_max_vv, cp=best_cp_vv)
  tree_ww <- prune(tree_max_ww, cp=best_cp_ww)

  ### compute goodness of fit: MSE(tau), MSE(C), LogL   =======================================================================
  #   compute predictions
  df_res_uu <- as.data.frame(tree_uu$where)                                     # leaf of each obs
  df_frame_uu <- as.data.frame(tree_uu$frame)                                   # tau estimation by leaf
  df_frame_uu$theta_est_uu <- iTau(claytonCopula(), tau=df_frame_uu$yval)       # theta estimation
  df_res_uu$y <- as.integer(rownames(df_res_uu))                                # observation index
  df_res_uu$tau_est_uu <- sapply(df_res_uu$`tree_uu$where`, function(n){df_frame_uu$yval[n]})           # tau est. for each obs
  df_res_uu$theta_est_uu <- sapply(df_res_uu$`tree_uu$where`, function(n){df_frame_uu$theta_est_uu[n]}) # theta ''          ''
  df_summary_uu <- df_res_uu[order(df_res_uu$y),]                               # sort observations

  df_res_vv <- as.data.frame(tree_vv$where)                                     # leaf of each obs
  df_frame_vv <- as.data.frame(tree_vv$frame)                                   # tau estimation by leaf
  df_frame_vv$theta_est_vv <- iTau(claytonCopula(), tau=df_frame_vv$yval)       # theta estimation
  df_res_vv$y <- as.integer(rownames(df_res_vv))                                # observation index
  df_res_vv$tau_est_vv <- sapply(df_res_vv$`tree_vv$where`, function(n){df_frame_vv$yval[n]})           # tau est. for each obs
  df_res_vv$theta_est_vv <- sapply(df_res_vv$`tree_vv$where`, function(n){df_frame_vv$theta_est_vv[n]}) # theta ''          ''
  df_summary_vv <- df_res_vv[order(df_res_vv$y),]                               # sort observations
  
  df_res_ww <- as.data.frame(tree_ww$where)                                     # leaf of each obs
  df_frame_ww <- as.data.frame(tree_ww$frame)                                   # tau estimation by leaf
  df_frame_ww$theta_est_ww <- iTau(claytonCopula(), tau=df_frame_ww$yval)       # theta estimation
  df_res_ww$y <- as.integer(rownames(df_res_ww))                                # observation index
  df_res_ww$tau_est_ww <- sapply(df_res_ww$`tree_ww$where`, function(n){df_frame_ww$yval[n]})           # tau est. for each obs
  df_res_ww$theta_est_ww <- sapply(df_res_ww$`tree_ww$where`, function(n){df_frame_ww$theta_est_ww[n]}) # theta ''          ''
  df_summary_ww <- df_res_ww[order(df_res_ww$y),]                               # sort observations

  df_summary <- cbind(df_cop[order(df_cop$y_idx),], 
                      df_summary_uu[,c('tau_est_uu','theta_est_uu')],
                      df_summary_vv[,c('tau_est_vv','theta_est_vv')], 
                      df_summary_ww[,c('tau_est_ww','theta_est_ww')])

  #   compute cumulative copulas
  df_summary$C <- apply(df_summary[,c('u1','u2','theta_x')], 1, function(x){pCopula(u=c(x[1],x[2]), copula=claytonCopula(param=x[3], dim=2))})
  df_summary$C_hat_uu <- apply(df_summary[,c('u1','u2','theta_est_uu')], 1, function(x){pCopula(u=c(x[1],x[2]), copula=claytonCopula(param=x[3], dim=2))})
  df_summary$C_hat_vv <- apply(df_summary[,c('u1','u2','theta_est_vv')], 1, function(x){pCopula(u=c(x[1],x[2]), copula=claytonCopula(param=x[3], dim=2))})
  df_summary$C_hat_ww <- apply(df_summary[,c('u1','u2','theta_est_ww')], 1, function(x){pCopula(u=c(x[1],x[2]), copula=claytonCopula(param=x[3], dim=2))})

  #   compute MSE(tau) and MSE(C)
  MSE_tau_uu <- sqrt(mean((df_summary$tau_x - df_summary$tau_est_uu)^2)) 
  MSE_tau_vv <- sqrt(mean((df_summary$tau_x - df_summary$tau_est_vv)^2)) 
  MSE_tau_ww <- sqrt(mean((df_summary$tau_x - df_summary$tau_est_ww)^2)) 
  MSE_C_uu <- sqrt(mean((df_summary$C - df_summary$C_hat_uu)^2)) 
  MSE_C_vv <- sqrt(mean((df_summary$C - df_summary$C_hat_vv)^2)) 
  MSE_C_ww <- sqrt(mean((df_summary$C - df_summary$C_hat_ww)^2)) 
  
  #   compute logL
  logL_uu <- 0
  for (theta_est in unique(df_summary$theta_est_uu)){
    df <- df_summary[df_summary$theta_est_uu==theta_est,]
    logL_uu <- logL_uu + loglikCopula(param=theta_est, u=as.matrix(df[, c('u1', 'u2')], ncol=2), copula=claytonCopula(dim=2))
  }
  logL_vv <- 0
  for (theta_est in unique(df_summary$theta_est_vv)){
    df <- df_summary[df_summary$theta_est_vv==theta_est,]
    logL_vv <- logL_vv + loglikCopula(param=theta_est, u=as.matrix(df[, c('v1', 'v2')], ncol=2), copula=claytonCopula(dim=2))
  }
  logL_ww <- 0
  for (theta_est in unique(df_summary$theta_est_ww)){
    df <- df_summary[df_summary$theta_est_ww==theta_est,]
    logL_ww <- logL_ww + loglikCopula(param=theta_est, u=as.matrix(df[, c('w1', 'w2')], ncol=2), copula=claytonCopula(dim=2))
  }

  ### Fit the benchmark model   =======================================================================================
  #   estimations
  tau_bench_uu <- df_frame_uu[df_frame_uu$n==N, c('yval')]
  theta_bench_uu <- df_frame_uu[df_frame_uu$n==N, c('theta_est_uu')]
  tau_bench_vv <- df_frame_vv[df_frame_vv$n==N, c('yval')]
  theta_bench_vv <- df_frame_vv[df_frame_vv$n==N, c('theta_est_vv')]
  tau_bench_ww <- df_frame_ww[df_frame_ww$n==N, c('yval')]
  theta_bench_ww <- df_frame_ww[df_frame_ww$n==N, c('theta_est_ww')]
  df_benchmark_ <- data.frame(tau_bench_uu=rep(tau_bench_uu, N),
                              theta_bench_uu=rep(theta_bench_uu, N),
                              tau_bench_vv=rep(tau_bench_vv, N),
                              theta_bench_vv=rep(theta_bench_vv, N),
                              tau_bench_ww=rep(tau_bench_ww, N),
                              theta_bench_ww=rep(theta_bench_ww, N))
  df_benchmark <- cbind(df_cop[order(df_cop$y_idx),], df_benchmark_)

  #   cumulative copulas
  df_benchmark$C <- apply(df_benchmark[,c('u1','u2','theta_x')], 1, function(x){pCopula(u=c(x[1],x[2]), copula=claytonCopula(param=x[3], dim=2))})
  df_benchmark$C_hat_uu <- apply(df_benchmark[,c('u1','u2','theta_bench_uu')], 1, function(x){pCopula(u=c(x[1],x[2]), copula=claytonCopula(param=x[3], dim=2))})
  df_benchmark$C_hat_vv <- apply(df_benchmark[,c('u1','u2','theta_bench_vv')], 1, function(x){pCopula(u=c(x[1],x[2]), copula=claytonCopula(param=x[3], dim=2))})
  df_benchmark$C_hat_ww <- apply(df_benchmark[,c('u1','u2','theta_bench_ww')], 1, function(x){pCopula(u=c(x[1],x[2]), copula=claytonCopula(param=x[3], dim=2))})
  #   MSE(tau), MSE(C) and logL
  MSE_tau_bench_uu <- sqrt(mean((df_benchmark$tau_x - df_benchmark$tau_bench_uu)^2)) 
  MSE_tau_bench_vv <- sqrt(mean((df_benchmark$tau_x - df_benchmark$tau_bench_vv)^2)) 
  MSE_tau_bench_ww <- sqrt(mean((df_benchmark$tau_x - df_benchmark$tau_bench_ww)^2)) 
  MSE_C_bench_uu <- sqrt(mean((df_benchmark$C - df_benchmark$C_hat_uu)^2)) 
  MSE_C_bench_vv <- sqrt(mean((df_benchmark$C - df_benchmark$C_hat_vv)^2)) 
  MSE_C_bench_ww <- sqrt(mean((df_benchmark$C - df_benchmark$C_hat_ww)^2)) 
  logL_bench_uu <- loglikCopula(param=theta_bench_uu, u=as.matrix(df_cop[, c('u1', 'u2')], ncol=2), copula=claytonCopula(dim=2))
  logL_bench_vv <- loglikCopula(param=theta_bench_vv, u=as.matrix(df_cop[, c('v1', 'v2')], ncol=2), copula=claytonCopula(dim=2))
  logL_bench_ww <- loglikCopula(param=theta_bench_ww, u=as.matrix(df_cop[, c('w1', 'w2')], ncol=2), copula=claytonCopula(dim=2))

  ### store results   =================================================================================================
  mat_output[D,] <- c(logL_bench_uu, MSE_tau_bench_uu, MSE_C_bench_uu, logL_uu, MSE_tau_uu, MSE_C_uu, best_nleaves_uu,
                      logL_bench_vv, MSE_tau_bench_vv, MSE_C_bench_vv, logL_vv, MSE_tau_vv, MSE_C_vv, best_nleaves_vv,
                      logL_bench_ww, MSE_tau_bench_ww, MSE_C_bench_ww, logL_ww, MSE_tau_ww, MSE_C_ww, best_nleaves_ww)
}
end_time <- Sys.time()
end_time - start_time

### save results
name_fileout <- sprintf('%s/%s_%s_D%d_N%d_taus%0.0f-%0.0f-%0.0f-%0.0f-CV%dx%d_max%dleaves.csv',outdir_path, copule, cov_combination, ndataset, N, taus[1]*100, taus[2]*100, taus[3]*100, taus[4]*100, kfolds, times_fold, max_nleaves)
df_output <- data.frame(mat_output)
colnames(df_output) <- mat_columns
write.table(df_output, name_fileout, sep = ',', dec = '.', row.names = F, col.names = T)



# ############################ per capire cosa fa pobs:   --------------------------------------------------
# set.seed(123456)
# nn <- 4
# mat1 <- matrix(rnorm(nn*2), ncol=2)
# mat2 <- matrix(rnorm(nn*2), ncol=2)
# pobs1 <- pobs(mat1)
# pobs2 <- pobs(mat2)
# par(mfrow=c(1,2))
# plot(pobs1)
# plot(pobs2)
# # columns of the 2 matrices xontain the same elements. That is, values in [0,1] spaced by 1/nn.
# # but they are in different order, so matrices pobs1 and pobs2 are different
# 
# 
# ######################### la log-likelihood è proporzionale al numero di osservazioni:  ---------------------------
# theta <- iTau(copula=claytonCopula(dim=2), tau=0.5)
# # ex: 
# uv100 <- rCopula(100, copula = claytonCopula(param=theta, dim=2))
# uv1000 <- rCopula(1000, copula = claytonCopula(param=theta, dim=2))
# loglikCopula(param=theta, u=uv100, copula=claytonCopula(dim=2))
# loglikCopula(param=theta, u=uv1000, copula=claytonCopula(dim=2))
# # ripetiamo l'esempio tante volte, vediamo che la logL per 1000 punti è 10 volte quella calcolata per 100 punti
# logL100 <- c()
# logL1000 <- c()
# for (i in 1:100){
#   uv100 <- rCopula(100, copula = claytonCopula(param=theta, dim=2))
#   uv1000 <- rCopula(1000, copula = claytonCopula(param=theta, dim=2))
#   logL100 <- c(logL100, loglikCopula(param=theta, u=uv100, copula=claytonCopula(dim=2)))
#   logL1000 <- c(logL1000, loglikCopula(param=theta, u=uv1000, copula=claytonCopula(dim=2)))
# }
# hist((logL1000)/10 - logL100)
# 
# 
# ###################### understanding the deviance score in the evaluation function: ------------------------------
# df_copL <- df_cop[df_cop$x2<0.75,]
# df_copR <- df_cop[df_cop$x2>=0.75,]
# nL <- nrow(df_copL)
# nR <- nrow(df_copR)
# devL <- eval_claytonCop(df_copL$y_idx)$deviance
# devR <- eval_claytonCop(df_copR$y_idx)$deviance
# dev_parent <- eval_claytonCop(df_cop$y_idx)$deviance
# true_gain_clayton <- (devL/nL + devR/nR - dev_parent/N)
# 
# true_gain_clayton <- (devL/nL + devR/nR - dev_parent/N)
# 
# split1 <- split_claytonCop(df_cop[order(df_cop$x2),]$y_idx)
# 
# plot(1:(N-1), split1$goodness, col='red', main='clayton', ylim=c(-0.005, 0.15)) #ylim=c(-0.005, clayton_gain+0.05))
# abline(v=which.max(split1$goodness), lty='solid', col='black')
# abline(v=nL, lty='dotted', col='blue')
# points(nL, true_gain_clayton, pch='X', cex=4, col='blue')
# 
