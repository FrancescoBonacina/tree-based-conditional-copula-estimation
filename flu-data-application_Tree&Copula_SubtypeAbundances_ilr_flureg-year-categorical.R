################################################################################
#
#   Tree & copulas to fit subtype abundances
#
################################################################################

rm(list=ls())

library(copula)
library(caret)
library(rpart)
library(rpart.plot)
library(glmnet)
library(fastDummies)

### additional sources
source('CART_PersonalizedSplits_CategoricalCovariates.R')

### functions to compute log-likelihood of copulas functions -------------------
logL_clayton <- function(yy){
  # Log-likelihood of a Clayton copula.
  # Computation is done via the inverse tau method.
  
  theta <- iTau(claytonCopula(dim=2), cor(yy, method="kendall")[1,2])
  logL <- loglikCopula(param=theta, u=yy, copula=claytonCopula(dim=2))
  return(logL)
}
logL_frank <- function(yy){
  # Log-likelihood of a Frank copula.
  # Computation is done via the inverse tau method.
  
  theta <- iTau(frankCopula(dim=2), cor(yy, method="kendall")[1,2])
  logL <- loglikCopula(param=theta, u=yy, copula=frankCopula(dim=2))
  if (is.na(logL)){logL <- -999999}   # if computation went wrong, fix logL to an arbitrary large negative value
  return(logL)
}
logL_gumbel <- function(yy){
  # Log-likelihood of a Gumbel copula.
  # Computation is done via the inverse tau method.
  
  theta <- iTau(gumbelCopula(dim=2), cor(yy, method="kendall")[1,2])
  logL <- loglikCopula(param=theta, u=yy, copula=gumbelCopula(dim=2))
  if (is.na(logL)){logL <- -999999}   # if computation went wrong, fix logL to an arbitrary large negative value
  return(logL)
}
logL_joe <- function(yy){
  # Log-likelihood of a Joe copula.
  # Computation is done via the inverse tau method.
  
  theta <- iTau(joeCopula(dim=2), cor(yy, method="kendall")[1,2])
  logL <- loglikCopula(param=theta, u=yy, copula=joeCopula(dim=2))
  if (is.na(logL)){logL <- -999999}   # if computation went wrong, fix logL to an arbitrary large negative value
  return(logL)
}
param_estimation <- function(yy){
  # estimation of the parameter of interest.
  
  tau <- cor(yy, method="kendall")[1,2] # 
  return (tau)
}
### ----------------------------------------------------------------------------

### set directory for output files
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
wd <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
outdir <- paste0('40_out')
dir.create(file.path(wd, outdir), showWarnings = FALSE)
outdir_path <- paste0(wd, '/', outdir)
main_directory <- gsub('/40_other_attempts_tree&copula_subtype_abundances', '', wd)



### 1) Data on subtype abundances   --------------------------------------------
data_years <- '2010-2019'    # '2010-2019', '2009-2021' 
th_annual_cases <- 50        # 50, 500
infile <- sprintf('%s/20_trees&copulas_categorical_covariates/22_out/df_subtype_compositions_%s_th%dcases_ilr-map.csv', main_directory, data_years, th_annual_cases)
df_flu <- read.csv(infile)
df_flu <- df_flu[, !(colnames(df_flu) %in% c('u','v'))]  # remove columns (u,v). They are not correct, we need to re-define them
str(df_flu)
#plot(df_flu[,'p'], df_flu[,'q'])

# define short labels for flu regions
fluregs <- sort(unique(df_flu$fluregion))
fluregs_new <- c('C-Am', 'C-As','E-Af','E-As','E-Eu','M-Af','N-Am','N-Af','N-Eu','O',
                 'SW-Eu','SE-As','S-Af','S-As','TempS-Am','TropS-Am','W-Af','W-As')
dict_fluregs <- data.frame(old_label=fluregs, new_label=fluregs_new)
df_flu$fluregion <- sapply(df_flu$fluregion, function(r){dict_fluregs[dict_fluregs$old_label==r,'new_label']})
#    how many countries per fluregion ?
# unique(df_flu[df_flu$fluregion=='E-Af', 'country'])
# unique_flureg <- c()
# for (r in fluregs_new){
#   df <- df_flu[df_flu$fluregion==r,]
#   unique_flureg <- c(unique_flureg, length(unique(df$country)))
# }
# df_flureg <- data.frame(flureg=fluregs_new, nb_countries=unique_flureg)

# define short labels for W.H.O regions
whoregs <- c("European Region of WHO", "Region of the Americas of WHO", "Western Pacific Region of WHO",
             "South-East Asia Region of WHO", "African Region of WHO", 
             "Eastern Mediterranean Region of WHO")
whoregs_new <- c('Eu', 'Am', 'W-P', 'SE-As', 'Af', 'E-M')
dict_whoregs <- data.frame(old_label=whoregs, new_label=whoregs_new)
df_flu$whoregion <- sapply(df_flu$whoregion, function(r){dict_whoregs[dict_whoregs$old_label==r,'new_label']})

#   define new regions, to reduce the number of modalities
reg_Eu <- c('E-Eu','N-Eu','SW-Eu')
reg_Am <- c('C-Am', 'N-Am','TempS-Am','TropS-Am')
reg_SWANA <- c('C-As','N-Af','S-As','W-As')                     # South-West Asia and North Africa
reg_AsAfO <- c('E-Af','E-As','M-Af','O','SE-As','S-Af','W-Af')
my_regs <- list(reg_Eu=reg_Eu, reg_Am=reg_Am, reg_SWANA=reg_SWANA, reg_AsAfO=reg_AsAfO)
my_regs_labels <- c('Eu','Am','SWANA','AsAfO')
df_flu$my_macroreg <- sapply(df_flu$fluregion, function(r){idx <- -1; for (i in 1:4){if (r %in% my_regs[[i]]){idx <- my_regs_labels[i]}}; idx})

# set covariates as factors (categorical variables)
df_flu$fluregion <- as.factor(df_flu$fluregion)
df_flu$whoregion <- as.factor(df_flu$whoregion)
df_flu$my_macroreg <- as.factor(df_flu$my_macroreg)
df_flu$year <- as.factor(df_flu$year)

# add index column
df_flu$y_idx <- 1:nrow(df_flu)
str(df_flu)

### 2) compute the uniform r.v. conditionally on the covariates ----------------
set.seed(231101)
df_p <- df_flu[, c('p', 'fluregion', 'year')]
df_q <- df_flu[, c('q', 'fluregion', 'year')]
#   We will consider 2 models: non-parametric and parametric
#   ____________________
#   A - regression tree 
#   ____________________
#
#   cross-validation, to set the optimal nb of leaves
max_nsplits <- 16
folds <- createMultiFolds(1:nrow(df_flu), k=3, times=50)   
CV_results_p <- matrix(NA, nrow=length(folds), ncol=max_nsplits)
CV_results_q <- matrix(NA, nrow=length(folds), ncol=max_nsplits)
for (f in 1:length(folds)){
  #f <- 1
  ### test and train
  fold <- folds[[f]]
  train <- df_flu[df_flu$y_idx %in% fold, ]
  test <- df_flu[!(df_flu$y_idx %in% fold),]
  
  ### predicting p
  # train the tree
  tree_complete_p <- rpart(data = train,
                           formula = p ~ fluregion + year,
                           control = rpart.control(cp=0.01, minsplit = 81),
                           xval = 0)
  # assess goodness of trees, from 0 to max number of splits
  cps_p <- tree_complete_p$cptable[,c('CP')]
  nsplits_p <- tree_complete_p$cptable[,c('nsplit')]

  for (i in 1:length(nsplits_p)){
    mynsplit <- nsplits_p[[i]]
    mycp <- cps_p[[i]]
    if (mynsplit<max_nsplits){  # if not too many split, consider the subtree and perform prediction on the test set
      # prune the tree
      mytree <- prune(tree_complete_p, cp=mycp, minsplit = 81)
      # predictions
      p_pred <- predict(mytree, test)
      # prediction error
      mse <- mean((test$p - p_pred)^2)
      # store error
      CV_results_p[f, mynsplit+1] <- mse
    }
  }
  
  ### predicting q
  # train the tree
  tree_complete_q <- rpart(data = train,
                           formula = q ~ fluregion + year,
                           control = rpart.control(cp=0.01, minsplit = 81),
                           xval = 0)
  # assess goodness of trees, from 0 to max number of splits
  cps_q <- tree_complete_q$cptable[,c('CP')]
  nsplits_q <- tree_complete_q$cptable[,c('nsplit')]
  
  for (i in 1:length(nsplits_q)){
    mynsplit <- nsplits_q[[i]]
    mycp <- cps_q[[i]]
    if (mynsplit<max_nsplits){  # if not too many split, consider the subtree and perform prediction on the test set
      # prune the tree
      mytree <- prune(tree_complete_q, cp=mycp, minsplit = 81)
      # predictions
      q_pred <- predict(mytree, test)
      # prediction error
      mse <- mean((test$q - q_pred)^2)
      # store error
      CV_results_q[f, mynsplit+1] <- mse
    }
  }
}
#   find best nb of split, create the optimal tree and store classification of obs. according to the optimal tree
min_nb_values <- 50
# p)
CV_mse_mean_p <- apply(CV_results_p, 2, function(x){if (sum(!is.na(x))>min_nb_values) mean(x, na.rm = TRUE) else NA})
while (is.na(CV_mse_mean_p[length(CV_mse_mean_p)])){CV_mse_mean_p <- head(CV_mse_mean_p,-1)}
CV_mse_sd_p <- apply(CV_results_p, 2, function(x) sd(x, na.rm = TRUE))[1:length(CV_mse_mean_p)] # remove all the NAs at the end
breiman_th_p <- CV_mse_mean_p[which.min(CV_mse_mean_p)] + CV_mse_sd_p[which.min(CV_mse_mean_p)]
best_nsplits_p <- which(CV_mse_mean_p<breiman_th_p)[1] # the smallest subtree that has a good performance (following the Breiman's rule)
plot(1:length(CV_mse_mean_p), CV_mse_mean_p, xlab = 'nb. of leaves', ylab = 'MSE', ylim=c(-5,10)) #ylim=c(0, max(CV_mse_mean_p+CV_mse_sd_p)+1))
arrows(x0=1:length(CV_mse_mean_p), y0=CV_mse_mean_p-CV_mse_sd_p, x1=1:length(CV_mse_mean_p), y1=CV_mse_mean_p+CV_mse_sd_p, code=3, angle=90, length=0.1)
points(best_nsplits_p, CV_mse_mean_p[best_nsplits_p], pch='X', cex=2, col='red')
abline(h=breiman_th_p, lty=2)
complete_tree_p <- rpart(data = df_flu,
                         formula = p ~ fluregion + year,
                         control = rpart.control(cp=0.01, minsplit = 81),
                         xval = 0)
cptab_p <- as.data.frame(complete_tree_p$cptable)
best_cp_p <- tail(cptab_p[cptab_p$nsplit<=best_nsplits_p, 'CP'], 1)
optimal_tree_p <- prune(complete_tree_p, cp=best_cp_p, minsplit = 81)
plot_tree_p <- sprintf('%s/optimal_tree_p(fluregion,year)_%s_th%dcases.png', outdir_path, data_years, th_annual_cases)
png(plot_tree_p, width=1500, height=600)
rpart.plot(optimal_tree_p)
dev.off()
df_flu$leaf_p <- optimal_tree_p$where
tree_pred_p <- predict(optimal_tree_p, newdata=df_flu, type='vector')
# q)
CV_mse_mean_q <- apply(CV_results_q, 2, function(x){if (sum(!is.na(x))>min_nb_values) mean(x, na.rm = TRUE) else NA})
while (is.na(CV_mse_mean_q[length(CV_mse_mean_q)])){CV_mse_mean_q <- head(CV_mse_mean_q,-1)}
CV_mse_sd_q <- apply(CV_results_q, 2, function(x) sd(x, na.rm = TRUE))[1:length(CV_mse_mean_q)] # remove all the NAs at the end
breiman_th_q <- CV_mse_mean_q[which.min(CV_mse_mean_q)] + CV_mse_sd_q[which.min(CV_mse_mean_q)]
best_nsplits_q <- which(CV_mse_mean_q<breiman_th_q)[1] # the smallest subtree that has a good performance (following the Breiman's rule)
complete_tree_q <- rpart(data = df_flu,
                         formula = q ~ fluregion + year,
                         control = rpart.control(cp=0.01, minsplit = 81),
                         xval = 0)
cptab_q <- as.data.frame(complete_tree_q$cptable)
best_cp_q <- tail(cptab_q[cptab_q$nsplit<=best_nsplits_q, 'CP'], 1)
optimal_tree_q <- prune(complete_tree_q, cp=best_cp_q, minsplit = 81)
plot_tree_q <- sprintf('%s/optimal_tree_q(fluregion,year)_%s_th%dcases.png', outdir_path, data_years, th_annual_cases)
png(plot_tree_q, width=1500, height=600)
rpart.plot(optimal_tree_q)
dev.off()
df_flu$leaf_q <- optimal_tree_q$where
tree_pred_q <- predict(optimal_tree_q, newdata=df_flu, type='vector')

# tree performance
r2_tree_p <- 1 - mean((df_flu$p - tree_pred_p)^2)/var(df_flu$p)
r2_tree_q <- 1 - mean((df_flu$q - tree_pred_q)^2)/var(df_flu$q)

#   compute the uniforms (u1,u2) from tree classifications.
#   uniforms are calculated as a mixture of ECDFs: in practice, for each
#   leaf of the tree, we consider the observations y_i within the leaf and 
#   we define u_i = ECDF(y_i)
df_flu$u1 <- NA
for (leaf_p in unique(df_flu$leaf_p)){
  df <- df_flu[df_flu$leaf_p==leaf_p,]
  ECDF <- ecdf(df$p)
  df_flu[df_flu$leaf_p==leaf_p, 'u1'] <- ECDF(df$p)
}
df_flu$u2 <- NA
for (leaf_q in unique(df_flu$leaf_q)){
  df <- df_flu[df_flu$leaf_q==leaf_q,]
  ECDF <- ecdf(df$q)
  df_flu[df_flu$leaf_q==leaf_q, 'u2'] <- ECDF(df$q)
}
df_flu$u1 <- df_flu$u1 - 0.0001  # to avoid value 1, which should appear with probability=0
df_flu$u2 <- df_flu$u2 - 0.0001

#   __________________________________________________________________
#   B - gaussian linear model with dummies variables for each modality
#   __________________________________________________________________

# create the matrix of covariates with dummies variables
X <- dummy_cols(df_flu[,c('my_macroreg','year')],  # all covariates we consider
                select_columns=c('my_macroreg','year'),   # covariates to be transformed into dummy variables
                remove_first_dummy = TRUE, 
                remove_selected_columns = TRUE)

# predict p and q, by using a linear model with lasso regularization
lm_cv_p <- cv.glmnet(x=as.matrix(X), y=df_flu$p, type.measure = "mse", nfolds = 10)
lm_p <- glmnet(x=as.matrix(X), y=df_flu$p) #, family = gamma(link="log")) 
lm_pred_p <- predict(lm_p, newx = as.matrix(X), type = "response", s = lm_cv_p$lambda.min) # lm_cv_p$lambda.1se
lm_cv_q <- cv.glmnet(x=as.matrix(X), y=df_flu$q, type.measure = "mse", nfolds = 10)
lm_q <- glmnet(x=as.matrix(X), y=df_flu$q)
lm_pred_q <- predict(lm_q, newx = as.matrix(X), type = "response", s = lm_cv_q$lambda.min)
#coef(lm_p, s = lm_cv_p$lambda.min)  # coefficients of the selected model
#coef(lm_q, s = lm_cv_q$lambda.min)

# performance of the linear model
r2_lm_p <- 1 - mean((df_flu$p - lm_pred_p)^2)/var(df_flu$p)
r2_lm_q <- 1 - mean((df_flu$q - lm_pred_q)^2)/var(df_flu$q)
#plot(1:length(df_flu$p), ((df_flu$p - lm_pred_p))/sd(df_flu$p))

# compute uniforms
v1 <- pnorm(df_flu$p, mean = lm_pred_p, sd = mean(abs(df_flu$p - lm_pred_p)))   
v2 <- pnorm(df_flu$q, mean = lm_pred_q, sd = mean(abs(df_flu$q - lm_pred_q)))   
df_flu$v1 <- 0.0001 + 0.9998*(v1 - min(v1))/(max(v1) - min(v1))   # rescale values in order to avoid zeros and ones
df_flu$v2 <- 0.0001 + 0.9998*(v2 - min(v2))/(max(v2) - min(v2))

# compute the log-likelihood of the root node (all observations together)
uu <- as.matrix(df_flu[,c('u1','u2')])
vv <- as.matrix(df_flu[,c('v1','v2')])
df_logL <- data.frame(copula=c('clayton','frank','gumbel','joe'),
                      logL_uu_tree = c(logL_clayton(uu), logL_frank(uu), logL_gumbel(uu), logL_joe(uu)),
                      logL_vv_lm = c(logL_clayton(vv), logL_frank(vv), logL_gumbel(vv), logL_joe(vv)))
outplot <- sprintf('%s/scatter_margins_uu_vv.png',outdir_path)
png(outplot, width = 750, height = 750)
par(mfrow=c(1,2), pty="s")
plot(df_flu$u1, df_flu$u2, asp=1, main='TREE estimation of the margins')
plot(df_flu$v1, df_flu$v2, asp=1, main='LM estimation of the margins')
dev.off()

# save dataframe with the pseudo-observations (u1,u2) and (v1,v2)
outfile_df_flu <- sprintf('%s/df_flu_pseudo-obs_%s_th%dcases.csv',outdir_path, data_years, th_annual_cases)
write.table(df_flu, outfile_df_flu, sep = ',', dec = '.', row.names = F, col.names = T)

### 3) CV procedure for tree pruning   -----------------------------------------
copula <- 'Gumbel'
if (copula=='Clayton'){logL_estimation <- logL_clayton}
if (copula=='Frank'){logL_estimation <- logL_frank}
if (copula=='Gumbel'){logL_estimation <- logL_gumbel}

#      re-set column type
df_flu$fluregion <- as.character(df_flu$fluregion)
#df_flu$my_macroreg <- as.character(df_flu$my_macroreg)
df_flu$year <- as.numeric(as.character(df_flu$year))
str(df_flu)

#      we need to decide the number of leaves
max_depth <- 5
response_var_uu <- c('u1','u2')
response_var_vv <- c('v1','v2')
covariates <- c('fluregion', 'year')
categorical_cov <- c(TRUE, TRUE)
N <- nrow(df_flu)

kfolds <- 3
cv_repetitions <- 50
folds <- createMultiFolds(1:N, k=kfolds, times=cv_repetitions)
max_nb_leaves <- 2^max_depth
CV_results_uu <- matrix(NA, nrow=length(folds), ncol=max_nb_leaves)
CV_results_vv <- matrix(NA, nrow=length(folds), ncol=max_nb_leaves)
for (f in 1:length(folds)){
  #f <- 13
  print (f)
  
  ### test and train
  fold <- folds[[f]]
  train <- df_flu[df_flu$y_idx %in% fold, ]
  test <- df_flu[!(df_flu$y_idx %in% fold),]
  
  ### train tree
  tree_maximale_uu <- build_max_tree(data=train, response_var = response_var_uu, covariates = covariates,
                                     categorical_cov = categorical_cov, max_depth = max_depth, 
                                     logL_func = logL_estimation, param_estimation = param_estimation,
                                     minbucket = 5)
  tree_maximale_vv <- build_max_tree(data=train, response_var = response_var_vv, covariates = covariates,
                                     categorical_cov = categorical_cov, max_depth = max_depth, 
                                     logL_func = logL_estimation, param_estimation = param_estimation,
                                     minbucket = 5)
  #labels <- sapply(tree_maximale_uu, function(t){t$label})
  #status <- sapply(tree_maximale_uu, function(t){t$status})
  #gains <- sapply(tree_maximale_uu, function(t){t$gain})
  #data.frame(label=labels, status=status, gain=gains)
  
  ### apply the tree on the test set
  predicted_tree_uu <- apply_max_tree(trained_tree = tree_maximale_uu, data_test = test, 
                                      response_var = response_var_uu, covariates = covariates,
                                      categorical_cov = categorical_cov, logL_func = logL_estimation,
                                      param_estimation = param_estimation, minbucket = 5)
  predicted_tree_vv <- apply_max_tree(trained_tree = tree_maximale_vv, data_test = test, 
                                      response_var = response_var_vv, covariates = covariates,
                                      categorical_cov = categorical_cov, logL_func = logL_estimation,
                                      param_estimation = param_estimation, minbucket = 5)
  #labels <- sapply(predicted_tree_uu, function(t){t$label})
  #status <- sapply(predicted_tree_uu, function(t){t$status})
  #gains <- sapply(predicted_tree_uu, function(t){t$gain})
  #data.frame(label=labels, status=status, gain=gains)
  
  ### tree pruning
  pruning_uu <- tree_pruning(max_tree = predicted_tree_uu, nleaves = 1)
  pruning_vv <- tree_pruning(max_tree = predicted_tree_vv, nleaves = 1)
  
  ### store subtree gains
  CV_results_uu[f, pruning_uu$subtree_nleaves] <- pruning_uu$subtree_gain
  CV_results_vv[f, pruning_vv$subtree_nleaves] <- pruning_vv$subtree_gain
}

# count nb of NA
#perc_na_uu <- c()
#perc_na_vv <- c()
#for (i in 1:6){
#  col_uu <- CV_results_uu[,i]
#  col_vv <- CV_results_vv[,i]
#  na_uu <- sum(is.na(col_uu))/length(col_uu)
#  na_vv <- sum(is.na(col_vv))/length(col_vv)
#  perc_na_uu <- c(perc_na_uu, na_uu)
#  perc_na_vv <- c(perc_na_vv, na_vv)
#}
#1 - perc_na_uu
#1 - perc_na_vv

# compare gains in logL per nb of leaves
min_nb_values <- cv_repetitions*kfolds/3

#View(CV_results_uu)
CV_results_uu[is.infinite(CV_results_uu)] <- NA
CV_means_uu <- apply(CV_results_uu, 2, function(x){if (sum(!is.na(x))>min_nb_values) mean(x, na.rm = TRUE) else NA})
while (is.na(CV_means_uu[length(CV_means_uu)])){CV_means_uu <- head(CV_means_uu,-1)}
CV_sd_uu <- apply(CV_results_uu, 2, function(x) sd(x, na.rm = TRUE))[1:length(CV_means_uu)] # remove all the NAs at the end
breiman_th_uu <- CV_means_uu[which.max(CV_means_uu)] - CV_sd_uu[which.max(CV_means_uu)]
best_nleaves_uu <- which(CV_means_uu>=breiman_th_uu)[1] # the smallest subtree that has a good performance (following the Breiman's rule)
#pruning_plot_uu <- sprintf('%s/Tree&%s_uu_X=%s-%s_CV_%dfoldx%d_%s_th%dcases.png', outdir_path, copula, covariates[1], covariates[2], kfolds, cv_repetitions, data_years, th_annual_cases)
#png(file=pruning_plot_uu, width=900, height=600)
plot(1:length(CV_means_uu), CV_means_uu, xlab = 'nb. of leaves', ylab = 'gain in logL', ylim=c(-15,12)) #, ylim=c(-5, max(CV_means_uu+CV_sd_uu)+5))
arrows(x0=1:length(CV_means_uu), y0=CV_means_uu-CV_sd_uu, x1=1:length(CV_means_uu), y1=CV_means_uu+CV_sd_uu, code=3, angle=90, length=0.1)
points(best_nleaves_uu, CV_means_uu[best_nleaves_uu], pch='X', cex=2, col='red')
abline(h=breiman_th_uu, lty=2)
#dev.off()

#View(CV_results_vv)
CV_results_vv[is.infinite(CV_results_vv)] <- NA
CV_means_vv <- apply(CV_results_vv, 2, function(x){if (sum(!is.na(x))>min_nb_values) mean(x, na.rm = TRUE) else NA})
while (is.na(CV_means_vv[length(CV_means_vv)])){CV_means_vv <- head(CV_means_vv,-1)}
CV_sd_vv <- apply(CV_results_vv, 2, function(x) sd(x, na.rm = TRUE))[1:length(CV_means_vv)] # remove all the NAs at the end
breiman_th_vv <- CV_means_vv[which.max(CV_means_vv)] - CV_sd_vv[which.max(CV_means_vv)]
best_nleaves_vv <- which(CV_means_vv>=breiman_th_vv)[1] # the smallest subtree that has a good performance (following the Breiman's rule)
#pruning_plot_vv <- sprintf('%s/Tree&%s_vv_X=%s-%s_CV_%dfoldx%d_%s_th%dcases.png', outdir_path, copula, covariates[1], covariates[2], kfolds, cv_repetitions, data_years, th_annual_cases)
#png(file=pruning_plot_vv, width=900, height=600)
plot(1:length(CV_means_vv), CV_means_vv, xlab = 'nb. of leaves', ylab = 'gain in logL', ylim=c(-30, max(CV_means_vv+CV_sd_vv)+30))
arrows(x0=1:length(CV_means_vv), y0=CV_means_vv-CV_sd_vv, x1=1:length(CV_means_vv), y1=CV_means_vv+CV_sd_vv, code=3, angle=90, length=0.1)
points(best_nleaves_vv, CV_means_vv[best_nleaves_vv], pch='X', cex=2, col='red')
abline(h=breiman_th_vv, lty=2)
#dev.off()

### optimal tree
complete_tree_uu <- build_max_tree(data=df_flu, response_var = response_var_uu, covariates = covariates,
                                   categorical_cov = categorical_cov, max_depth = max_depth, 
                                   logL_func = logL_clayton, param_estimation = param_estimation_clayton,
                                   minbucket = 5)
final_pruning_uu <- tree_pruning(max_tree = complete_tree_uu, nleaves = best_nleaves_uu)
df_optimal_tree_uu <- final_pruning_uu$final_subtree
name_fileout_uu <- sprintf('%s/df_%s_optimal_tree_uu_x=%s-%s_%s_th%dcases.csv',outdir_path, copula, covariates[1], covariates[2], data_years, th_annual_cases)
write.table(df_optimal_tree_uu, name_fileout_uu, sep = ',', dec = '.', row.names = F, col.names = T)

complete_tree_vv <- build_max_tree(data=df_flu, response_var = response_var_vv, covariates = covariates,
                                   categorical_cov = categorical_cov, max_depth = max_depth, 
                                   logL_func = logL_clayton, param_estimation = param_estimation_clayton,
                                   minbucket = 5)
final_pruning_vv <- tree_pruning(max_tree = complete_tree_vv, nleaves = best_nleaves_vv)
df_optimal_tree_vv <- final_pruning_vv$final_subtree
name_fileout_vv <- sprintf('%s/df_%s_optimal_tree_vv_x=%s-%s_%s_th%dcases.csv',outdir_path, copula, covariates[1], covariates[2], data_years, th_annual_cases)
write.table(df_optimal_tree_vv, name_fileout_vv, sep = ',', dec = '.', row.names = F, col.names = T)

sum(df_optimal_tree_uu[df_optimal_tree_uu$status %in% c('terminal', 'parent'), 'gain'])
sum(df_optimal_tree_vv[df_optimal_tree_vv$status %in% c('terminal', 'parent'), 'gain'])

