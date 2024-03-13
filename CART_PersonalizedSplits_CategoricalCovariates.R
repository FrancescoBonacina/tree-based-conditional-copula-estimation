################################################################################
#
#   Implementation of the CART algorithm.
#   The splits are performed such as to maximize a log-likelihood 
#   function. The splitting procedure allows the use of categorical covariates.
# 
################################################################################

logL_bigauss <- function(yy, Sig=matrix(c(0.002,0,0,0.002), 2)){
  #   Log-Likelihood of a bivariate normal distribution,
  #   with fixed covariance matrix (Sig).
  
  Inv <- chol2inv(chol(Sig))
  logdet <- determinant(Sig, logarithm = TRUE)$modulus[1]
  n <- nrow(yy)
  mu <- colMeans(yy)
  resid <- sweep(yy, 2, mu)
  logL <- -0.5*n*(2*log(2*pi) + logdet + mean(sapply(1:n, function(i){t(resid[i,]) %*% Inv %*% resid[i,]})))
  return(logL)
}


param_estimation_bigauss <- function(yy){
  # estimation of the parameter of interest assuming that response variables
  # (yy) are distributed according to a bivariate normal distribution of 
  # center (mu,mu) and covariance matrix Sig=matrix(c(0.002,0,0,0.002), 2).
  
  mu <- mean(colMeans(yy))
  return (mu)
}


logL_clayton <- function(yy){
  # Log-likelihood of a Clayton copula.
  # Compuation is done via the inverse tau method.
  
  theta <- iTau(claytonCopula(dim=2), cor(yy, method="kendall")[1,2])
  logL <- loglikCopula(param=theta, u=yy, copula=claytonCopula(dim=2))
  return(logL)
}


param_estimation_clayton <- function(yy){
  # estimation of the parameter of interest, assuming that the response
  # variables (yy) are distributed according to a Clayton copula of parameter 
  # tau (i.e. the correlation coef. that is directly associated to the 
  # copula parameter theta via a bijective relationship) 

  tau <- cor(yy, method="kendall")[1,2] # 
  return (tau)
}


ordering_step <- function(data, x_categorical, response_var){
  # compute an auxiliary variable to be used in the splitting procedure 
  # in place of the covariate x_categorical.
  # Following suggestion in section 5 of the "User-written split functions 
  # for RPART" vignette, the classes of the categorical covariate are sorted
  # according to the parameter of interest. Thus, the auxiliary covariate 
  # will contain values from 1 to nb of classes, corresponding to this ordering
  # of the classes. In the splitting procedure, (nb of classes -1) splits will 
  # be testes along the auxiliary covariate.
  
  # all existing classes
  classes <- sort(unique(data[, x_categorical]))
  
  # parameter estimation per class
  param_by_class <- c()
  for (x in classes){
    dfx <- data[data[,x_categorical]==x,]
    param_by_class <- c(param_by_class, mean(as.matrix(dfx[, response_var], ncol=2)))   # parameter estimation
  }
  
  # define the auxiliary covariate, with classes sorted by the estimated param.
  df_dict <- data.frame(class=classes, param_by_class=param_by_class)
  df_dict <- df_dict[order(df_dict$param_by_class),]
  df_dict$idx <- 1:nrow(df_dict)
  x1_cat_aux <- sapply(data[,x_categorical], function(cl) df_dict[df_dict$class==cl,c('idx')])
  
  # save list of sorted classes
  classes_sorted <- paste0(c(x_categorical, df_dict$class), collapse='.')
  
  return(list(aux_covariate=x1_cat_aux, classes_sorted=classes_sorted))
}


create_triangle <- function(parent_triangle=NULL, brench_direction=NULL, parent_data=NULL, 
                           response_var=NULL, covariates=NULL, categorical_cov=NULL, 
                           logL_func=NULL, param_estimation=NULL, minbucket=1){
  # Function to create a triangle to be hanged to the parent triangle.
  # The child triangle will be hanged on the left (L) or on the right (R), according
  # to the brench_direction.
  # To create a child triangle we consider data from the parent node and we perform
  # a split research in order to infer the response variables, by considering the covariates.
  # Some covariates might be categorical, in this case an additional ordering-step is 
  # performed before the split research.
  # If the parent triangle exists and the brench direction is set, the information about
  # the parent node are taken from the parent triangle. Otherwise, the parent_data argument is
  # used.
  # 
  # Inputs: 
  # parent_triangle: list of information defining the parent triangle, as returned by
  #                  the function create_triangle().
  # brench_direction: 'L' or 'R', defining if the new triangle will be hanged on the
  #                    left/right vertex of the parent triangle
  # parent_data: dataframe containing the row observations. To be used instead of the 
  #              parent_triangle and brench_direction arguments.
  # response_var: list of labels, to identify the columns of data which correspond to 
  #               the response variables
  # covariates: list of columns to be used as covariates to build the tree
  # categorical_cov: list of booleans to specify if each covariate is categorical or not
  # logL_func: function to compute the logL that we want maximize
  # param_estimation: function to estimate the parameter of interest
  # minbucket: integer, the minimum nb. of obs. in any leaf (as in RPART)
  #
  # Output:
  # new_triangle: list defining the new triangle.
  
  ### info of the parent node
  if (!is.null(parent_triangle)){
    depth <- parent_triangle$depth + 1
    label <- paste0(parent_triangle$label, brench_direction)
    if (brench_direction=='L'){
      data.P <- parent_triangle$data.L
      n.P <- parent_triangle$n.L
      logL.P <- parent_triangle$logL.L
      param.P <- parent_triangle$param.L
    } else if (brench_direction=='R'){
      data.P <- parent_triangle$data.R
      n.P <- parent_triangle$n.R
      logL.P <- parent_triangle$logL.R
      param.P <- parent_triangle$param.R
    }
  } else {
    depth <- 1
    label <- 'T'
    data.P <- data.frame(parent_data)
    n.P <- nrow(data.P)
    yy.P <- as.matrix(data.P[, response_var])
    logL.P <- logL_func(yy.P)
    param.P <- param_estimation(yy.P)
  }
  
  ### ordering step: sort categorical covariates
  classes_sorted <- c()
  covariates_to_test <- covariates   # list of covariates that will be considered for the split
  for (i in 1:length(covariates)){
    if (categorical_cov[i]){
      cov_label <- paste0(covariates[i], '_aux')
      covariates_to_test[i] <- cov_label
      sorted_cov <- ordering_step(data=data.P, x_categorical=covariates[i], response_var = response_var)
      data.P[,cov_label] <- sorted_cov$aux_covariate
      classes_sorted <- c(classes_sorted, sorted_cov$classes_sorted)
    }
  }
  
  ### split research
  gain <- 0             # store the gain in logL
  best_cov <- -999      # covariate used for the split
  best_th <- -999       # threshold, i.e. the location of the split for quantitative covariates or the list of classes that goes to the left for the categorical covariates
  direction <- -999     # used to keep on the left obs for which the estimated parameter is smaller
  data.A <- NULL        # info of the children nodes
  data.B <- NULL
  param.A <- NULL
  param.B <- NULL
  logL.A <- NULL
  logL.B <- NULL
  # test all the covariates
  for (cov in covariates_to_test){
    # sort data
    data.P <- data.P[order(data.P[,cov]),]
    unique_cov <- unique(data.P[,cov])
    
    # test all the possible splits
    goodness <- c(0)
    if (length(unique_cov)>1){
      for (t in unique_cov[2:length(unique_cov)]){
        # groups
        yyA <- as.matrix(data.P[data.P[,cov]<t, response_var], ncol=2)
        yyB <- as.matrix(data.P[data.P[,cov]>=t, response_var], ncol=2)
        
        # if children do not contain a minimum nb. of observations, set the gain
        # to zero, otherwise compute the gain in logL
        if (nrow(yyA)<minbucket | nrow(yyB)<minbucket){
          logL_gain <- 0
        } else {
          logL_gain <- (logL_func(yyA) + logL_func(yyB) - logL.P)
        }
        
        goodness <- c(goodness, logL_gain)
      }
    }
    
    # if there is a convenient split, store info
    #cat(parent_triangle$label, ':   MAX= ', max(goodness), 'goodness= ', goodness)    ###
    if (max(goodness)>gain){
      gain <- max(goodness)
      best_cov <- cov
      best_th <- unique_cov[which.max(goodness)]
      
      # info of the children nodes
      data.A <- data.P[data.P[,cov]<best_th,]
      data.B <- data.P[data.P[,cov]>=best_th,]
      yyA <- as.matrix(data.A[,response_var], ncol=2)
      yyB <- as.matrix(data.B[,response_var], ncol=2)
      param.A <- param_estimation(yyA)
      param.B <- param_estimation(yyB)
      direction <- sign(sign(param.A-param.B)-0.5) # repeated use of the sign() function to avoid zeros
      logL.A <- logL_func(yyA)
      logL.B <- logL_func(yyB)
    }
  }
  
  ### assign info to the left/right child
  if (direction<0){
    data.L <- data.A
    data.R <- data.B
    param.L <- param.A
    param.R <- param.B
    logL.L <- logL.A
    logL.R <- logL.B
  } else if (direction>0) {
    data.L <- data.B
    data.R <- data.A
    param.L <- param.B
    param.R <- param.A
    logL.L <- logL.B
    logL.R <- logL.A
  }
  n.L <- nrow(data.L)
  n.R <- nrow(data.R)
  
  ### fix the threshold format for the categorical covariate (only if the split occurred)
  if (gain>0){
    if (categorical_cov[match(best_cov, covariates_to_test)]){
      classes_ <- unlist(strsplit(classes_sorted[match(best_cov, covariates_to_test)], '.', fixed=TRUE))
      vec_classes <- classes_[2:length(classes_)]
      best_th <- vec_classes[1:(best_th-1)]
    }
  }
  
  ### set outputs
  status <- 'non_convenient'
  has_children <- c(0,0)
  if (gain>0){
    status <- 'terminal'
  }
  new_triangle <- list(label=label,
                       depth=depth,
                       has_children=has_children,
                       status=status,
                       gain=gain,
                       best_cov=best_cov,
                       best_th=best_th,
                       direction=direction,
                       data.P=data.P,
                       data.L=data.L,
                       data.R=data.R,
                       param.P=param.P,
                       param.L=param.L,
                       param.R=param.R,
                       logL.P=logL.P,
                       logL.L=logL.L,
                       logL.R=logL.R,
                       n.P=n.P,
                       n.L=n.L,
                       n.R=n.R)
  
  return(new_triangle)
}


create_empty_triangle <- function(parent_triangle, brench_direction){
  # create an empty triangle:
  # when there is no convenient split, we simply create a triangle set to a 
  # 'non_convenient' status, without computing further info.
  # Specifically, an empty triangle keeps info about the label, the depth,
  # the status, the presence of children triangles and it has gain=0. The rest 
  # is set to NULL.
  
  new_triangle <- list(label=paste0(parent_triangle$label, brench_direction),
                       depth=parent_triangle$depth + 1,
                       has_children=c(0,0),
                       status='non_convenient',
                       gain=0,
                       best_cov=NULL,
                       best_th=NULL,
                       direction=NULL,
                       data.P=NULL,
                       data.L=NULL,
                       data.R=NULL,
                       param.P=NULL,
                       param.L=NULL,
                       param.R=NULL,
                       logL.P=NULL,
                       logL.L=NULL,
                       logL.R=NULL,
                       n.P=NULL,
                       n.L=NULL,
                       n.R=NULL)
  return(new_triangle)
}


build_max_tree <- function(data, response_var, covariates, categorical_cov, max_depth, logL_func,
                           param_estimation, minbucket=1){
  # Build a regression tree developed until the maximal depth.
  # The tree consists in a list of triangles, which are the building blocks of the tree.
  # The maximal depth is referred to the triangles, as follow:
  #    max_depth=1 : 1 triangles, 2 leaves, 3 nodes
  #    max_depth=2 : 3 triangles, 4 leaves, 7 nodes
  #    max_depth=3 : 7 triangles, 8 leaves, 15 nodes
  #    ...
  #
  # Inputs:
  # data: dataframe of data, included response variables and covariates
  # response_var: list of labels, to identify the columns of data which correspond to 
  #               the response variables
  # covariates: list of columns to be used as covariates to build the tree
  # categorical_cov: list of booleans to specify if each covariate is categorical or not
  # max_depth: integer, must be >0
  # logL_func: function to compute the log-likelihood that we want maximize
  # param_estimation: function to estimate the parameter of interest
  # minbucket: integer, the minimum nb. of obs. in any leaf (as in RPART)
  #
  # Output:
  # tree: list of (2^max_depth -1) triangles
  
  tree <- list()
  count_triangles <- 0
  for (depth in 1:max_depth){
    
    # for depth=1, create the first triangle from the row data
    if (depth==1){
      count_triangles <- count_triangles + 1
      tree[[count_triangles]] <- create_triangle(parent_data = data, response_var = response_var,
                                                covariates = covariates, categorical_cov = categorical_cov,
                                                logL_func = logL_func, param_estimation = param_estimation,
                                                minbucket = minbucket)
    } else {
      # indexes of triangles of the previous depth level, i.e. indexes of the parent triangles
      parent_indexes <- (2^(depth-2)):(2^(depth-1)-1)
      
      # for each parent triangle, create children
      for (idx in parent_indexes){
        
        # get parent triangle
        parent_triangle <- tree[[idx]]

        # if parent is a terminal node, create children, else create 'empty' children
        if (parent_triangle$status=='terminal'){
          triangle.L <- create_triangle(parent_triangle = parent_triangle, brench_direction = 'L',
                                       response_var = response_var, covariates = covariates, 
                                       categorical_cov = categorical_cov, logL_func = logL_func, 
                                       param_estimation = param_estimation, minbucket = minbucket)
          triangle.R <- create_triangle(parent_triangle = parent_triangle, brench_direction = 'R',
                                       response_var = response_var, covariates = covariates, 
                                       categorical_cov = categorical_cov, logL_func = logL_func,
                                       param_estimation = param_estimation, minbucket = minbucket)
        } else {
          triangle.L <- create_empty_triangle(parent_triangle = parent_triangle, brench_direction = 'L')
          triangle.R <- create_empty_triangle(parent_triangle = parent_triangle, brench_direction = 'R')
        }
        
        # add children triangles to the tree
        count_triangles <- count_triangles + 1
        tree[[count_triangles]] <- triangle.L
        count_triangles <- count_triangles + 1
        tree[[count_triangles]] <- triangle.R
        
        # update info of the parent triangle
        if (triangle.L$gain>0 | triangle.R$gain>0){
          tree[[idx]]$status <- 'parent'
          tree[[idx]]$has_children <- c(1*(triangle.L$gain>0),1*(triangle.R$gain>0))
        }
      }
    }
  }
  return(tree)
}


apply_max_tree <- function(trained_tree, data_test, response_var, covariates,
                            categorical_cov, logL_func, param_estimation, minbucket=1){
  # Apply the trained (maximal) tree to the observations of the test set.
  #
  # Inputs:
  # trained_tree: list of triangles that defines the (maximal) tree trained on 
  #               the training set
  # data_test: dataframe with obs of the test set
  # response_var: list of labels, to identify the columns of data which 
  #               correspond to the response variables
  # covariates: list of columns to be used as covariates to build the tree
  # categorical_cov: list of booleans to specify if each covariate is categorical or not
  # logL_func: function to compute the log-likelihood than we want to maximize
  # param_estimation: function to compute the parameter of interest
  # minbucket: integer, the minimum nb. of obs. in any leaf (as in RPART)
  #
  # Output:
  # predicted_tree: list of triangles, of the same length of the trained tree
  
  predicted_tree <- list()
  nb_triangles <- length(trained_tree)
  
  for (t in 1:nb_triangles){
    
    ### general info of the triangle
    triangle_train <- trained_tree[[t]]
    status <- triangle_train$status
    depth <- triangle_train$depth
    label <- triangle_train$label
    best_cov <- triangle_train$best_cov
    best_th <- triangle_train$best_th
    
    ### if status=='non-convenient', we won't perform the split, otherwise we will split the obervations
    if (status=='non_convenient'){
      current_status <- 'non_convenient'        # if the triangle status of the trained tree is 'not_convenient', also the triangle status in the predicted tree will be 'non_convenient'
      triangle_test <- list(label=label,
                            depth=depth,
                            has_children=c(0,0),
                            status=current_status,
                            gain=0,
                            best_cov=NULL,
                            best_th=NULL,
                            direction=NULL,
                            data.P=NULL,
                            data.L=NULL,
                            data.R=NULL,
                            param.P=NULL,
                            param.L=NULL,
                            param.R=NULL,
                            logL.P=NULL,
                            logL.L=NULL,
                            logL.R=NULL,
                            n.P=NULL,
                            n.L=NULL,
                            n.R=NULL)
    } else if (status=='parent' | status=='terminal') {
      
      ### get data to be splitted
      if (depth==1){              # to create the first triangle, consider all the test set
        data.P <- data_test
        yyP <- as.matrix(data.P[,response_var], ncol=2)
        param.P <- param_estimation(yyP)
        logL.P <- logL_func(yyP)
        n.P <- nrow(data.P)
      } else {                    # otherwise consider the obs of the parent node (only if the parent node is a terminal/parent node)
        label_parent <- substr(label, 1, nchar(label)-1)
        brench_direction <- substr(label, nchar(label), nchar(label))
        idx_parent <- which(sapply(predicted_tree, function(t){t$label})==label_parent)
        status_parent <- predicted_tree[[idx_parent]]$status
        if (status_parent=='terminal' | status_parent=='parent'){
          if (brench_direction=='L'){
            data.P <- predicted_tree[[idx_parent]]$data.L
            n.P <- predicted_tree[[idx_parent]]$n.L
            param.P <- predicted_tree[[idx_parent]]$param.L
            logL.P <- predicted_tree[[idx_parent]]$logL.L
          } else if (brench_direction=='R'){
            data.P <- predicted_tree[[idx_parent]]$data.R
            n.P <- predicted_tree[[idx_parent]]$n.R
            param.P <- predicted_tree[[idx_parent]]$param.R
            logL.P <- predicted_tree[[idx_parent]]$logL.R
          }
        } else if (status_parent=='not_enough_obs'){    # if the parent triangle is 'not_enough_obs', stop the split procedure
          data.P <- NULL
        }
      } 
      
      ### perform the split
      current_status <- 'not_enough_obs'
      if (!is.null(data.P)){
        if (grepl('aux', best_cov)){                        # if the covariate is categorical
          cov <- substr(best_cov, 1, nchar(best_cov)-4)     # go back to the original label of the covariate
          data.L <- data.P[data.P[,cov] %in% best_th,]
          if (nrow(data.L)>=minbucket & nrow(data.L)<=(n.P-minbucket)){ # if the children nodes contain at least 'minbucket' obs, perform the split
            data.R <- data.P[!data.P[,cov] %in% best_th,]
            yyL <- as.matrix(data.L[,response_var], ncol=2)
            yyR <- as.matrix(data.R[,response_var], ncol=2)
            param.L <- param_estimation(yyL)
            param.R <- param_estimation(yyR)
            direction <- sign(sign(param.L-param.R)-0.5)   
            logL.L <- logL_func(yyL)
            logL.R <- logL_func(yyR)
            n.L <- nrow(data.L)
            n.R <- nrow(data.R)
            gain <- logL.L + logL.R - logL_func(as.matrix(data.P[,response_var], ncol=2))
            current_status <- 'terminal'
          }
        } else {                                            # if the covariate is quantitative
          data.A <- data.P[data.P[,best_cov]<best_th,]
          if (nrow(data.A)>=minbucket & nrow(data.A)<(n.P-minbucket)){ # if the children nodes contain at least minbucket obs, perform the split
            data.B <- data.P[data.P[,best_cov]>=best_th,]
            yyA <- as.matrix(data.A[,response_var], ncol=2)
            yyB <- as.matrix(data.B[,response_var], ncol=2)
            param.A <- param_estimation(yyA)
            param.B <- param_estimation(yyB)
            direction <- sign(sign(param.A-param.B)-0.5) # repeated use of the sign() function to avoid zeros
            logL.A <- logL_func(yyA)
            logL.B <- logL_func(yyB)
            if (direction<0){
              data.L <- data.A
              data.R <- data.B
              param.L <- param.A
              param.R <- param.B
              logL.L <- logL.A
              logL.R <- logL.B
            } else if (direction>0) {
              data.L <- data.B
              data.R <- data.A
              param.L <- param.B
              param.R <- param.A
              logL.L <- logL.B
              logL.R <- logL.A
            }
            n.L <- nrow(data.L)
            n.R <- nrow(data.R)
            gain <- logL.L + logL.R - logL_func(as.matrix(data.P[,response_var], ncol=2))
            current_status <- 'terminal'
          }
        }
      }
      
      ### define triangle and store it
      if (current_status=='terminal'){
        triangle_test <- list(label=label,
                              depth=depth,
                              has_children=c(0,0),
                              status=current_status,
                              gain=gain,
                              best_cov=best_cov,
                              best_th=best_th,
                              direction=direction,
                              data.P=data.P,
                              data.L=data.L,
                              data.R=data.R,
                              param.P=param.P,
                              param.L=param.L,
                              param.R=param.R,
                              logL.P=logL.P,
                              logL.L=logL.L,
                              logL.R=logL.R,
                              n.P=n.P,
                              n.L=n.L,
                              n.R=n.R)
      } else if (current_status=='not_enough_obs'){
        triangle_test <- list(label=label,
                              depth=depth,
                              has_children=c(0,0),
                              status=current_status,
                              gain=0,
                              best_cov=NULL,
                              best_th=NULL,
                              direction=NULL,
                              data.P=NULL,
                              data.L=NULL,
                              data.R=NULL,
                              param.P=NULL,
                              param.L=NULL,
                              param.R=NULL,
                              logL.P=NULL,
                              logL.L=NULL,
                              logL.R=NULL,
                              n.P=NULL,
                              n.L=NULL,
                              n.R=NULL)
      }
    }
    
    ### add triangle to the tree
    predicted_tree[[t]] <- triangle_test
    
    ### update info of the parent triangle
    if (current_status=='terminal' & depth>1){
      predicted_tree[[idx_parent]]$status <- 'parent'
      if (brench_direction=='L'){
        predicted_tree[[idx_parent]]$has_children <- predicted_tree[[idx_parent]]$has_children + c(1,0)
      } else if (brench_direction=='R'){
        predicted_tree[[idx_parent]]$has_children <- predicted_tree[[idx_parent]]$has_children + c(0,1)
      }
    }
  }
  
  return(predicted_tree)
}


tree_pruning <- function(max_tree, nleaves){
  # function to prune the maximal tree, in order to get the optimal subtree with
  # number of leaves equals to 'nleaves'.
  # The tree is pruned iteratively, one triangle at the time. Following the Breiman 
  # procedure, at each iteration the pair of leaves which provides the smallest 
  # gain in log-likelihood is removed.
  # 
  # Inputs:
  # max_tree: list of triangles defining the maximal tree, as returned by the 
  #           functions 'build_max_tree()' and 'build_predicted_max_tree()'
  # nleaves: integer, the number of leaves of the final subtree
  #
  # Output, list containing 3 objects:
  # subtree_nleaves: vector of integers = nleaves of the maximal_tree:nleaves
  # subtree_gain: vector of the gains in logL of the optimal subtrees identified
  #               by the pruning procedure. That is, all the subtrees from the 
  #               maximal_tree itself, to the final subtree.
  # final_subtree: dataframe summarizing the information of the final subtree. 
  #                Rows of the dataframe correspond to the triangles that 
  #                constitute that tree
  
  # organize info of the tree in a dataframe
  df_subtree <- data.frame(label=sapply(max_tree, function(t){t$label}),
                           status=sapply(max_tree, function(t){t$status}),
                           has_child_L=sapply(max_tree, function(t){t$has_children[1]}),
                           has_child_R=sapply(max_tree, function(t){t$has_children[2]}),
                           gain=sapply(max_tree, function(t){t$gain}),
                           depth=sapply(max_tree, function(t){t$depth}),
                           split_cov=sapply(max_tree, function(t){cov=t$best_cov; if(!is.null(cov)) cov else '/'}),
                           split_th=sapply(max_tree, function(t){th=t$best_th; if(!is.null(th)) paste0(th, collapse=';') else '/'}),
                           split_direction=sapply(max_tree, function(t){dir=t$direction; if(!is.null(dir)) dir else '/'}),
                           n.P=sapply(max_tree, function(t){n=t$n.P; if(!is.null(n)) n else '/'}),
                           n.L=sapply(max_tree, function(t){n=t$n.L; if(!is.null(n)) n else '/'}),
                           n.R=sapply(max_tree, function(t){n=t$n.R; if(!is.null(n)) n else '/'}),
                           param.P=sapply(max_tree, function(t){p=t$param.P; if(!is.null(p)) p else '/'}),
                           param.L=sapply(max_tree, function(t){p=t$param.L; if(!is.null(p)) p else '/'}),
                           param.R=sapply(max_tree, function(t){p=t$param.R; if(!is.null(p)) p else '/'}),
                           logL.P=sapply(max_tree, function(t){ll=t$logL.P; if(!is.null(ll)) ll else '/'}),
                           logL.L=sapply(max_tree, function(t){ll=t$logL.L; if(!is.null(ll)) ll else '/'}),
                           logL.R=sapply(max_tree, function(t){ll=t$logL.R; if(!is.null(ll)) ll else '/'}))
  nleaves_max <- nrow(df_subtree[df_subtree$status %in% c('terminal','parent'),]) + 1
  
  # pruning step
  subtree_nleaves <- nleaves_max:nleaves
  subtree_gain <- c(sum(df_subtree[df_subtree$status %in% c('terminal','parent'), 'gain']))
  #               N.B. the gain in logL given by a tree is the sum of the gains provided by each triangle
  
  if (nleaves_max>1){
    for (nl in subtree_nleaves[2:length(subtree_nleaves)]){
      #message <- sprintf('*** max nleaves = %d,    nl = %d', nleaves_max, nl)
      #print (message)
      
      # identify the triangle to be pruned and change its status
      df_terminals <- df_subtree[df_subtree$status=='terminal',]
      triangle_to_prune <- df_terminals[order(df_terminals$gain), 'label'][1]
      df_subtree[df_subtree$label==triangle_to_prune, "status"] <- 'pruned'
      
      # update info of the parent triangle (except if I'm already at the root node)
      if (nl>1){ 
        #label <- df_subtree[(df_subtree$status=='pruned' & df_subtree$gain==min_gain), "label"]
        parent_label <- substr(triangle_to_prune, 1, nchar(triangle_to_prune)-1)
        has_child <- paste0('has_child_', substr(triangle_to_prune, nchar(triangle_to_prune), nchar(triangle_to_prune)))
        df_subtree[df_subtree$label==parent_label, has_child] <- 0
        if (sum(df_subtree[df_subtree$label==parent_label, c('has_child_L', 'has_child_R')])==0){
          df_subtree[df_subtree$label==parent_label, 'status'] <- 'terminal'
        }
      }
      
      # compute the subtree gain
      if (nl != nrow(df_subtree[df_subtree$status %in% c('terminal','parent'),]) + 1){print('ERROR: wrong computation of the nb. of leaves')}
      subtree_gain <- c(subtree_gain, sum(df_subtree[df_subtree$status %in% c('terminal','parent'), 'gain']))
    }
  }
  
  return(list(subtree_nleaves=subtree_nleaves, subtree_gain=subtree_gain, final_subtree=df_subtree))
}


# ------------------------------------------------------------------------------
# apply_max_tree_OLD <- function(trained_tree, data_test, response_var, covariates,
#                                categorical_cov, logL_func, param_estimation, minbucket=1){
#   # Apply the trained (maximal) tree to the observations of the test set.
#   #
#   # *** OLD IMPLEMENTATION: here I apply the trained tree on the data of the test set (that is, 
#   #     I apply the rules that I've learned), but I perform the splits only if there is a gain
#   #     in logL. However, what I should do (and I do in the NEW version of the function) is to 
#   #     split the data anyway and then, a posteriori, evaluate all the nested subtrees in a 
#   #     backward fashion, in order to find the optimal subtree.
#   #     In fact, it can happen that a split results not convenient, but the following splits do.
#   #
#   # Inputs:
#   # trained_tree: list of triangles that defines the (maximal) tree trained on 
#   #               the training set
#   # data_test: dataframe with obs of the test set
#   # response_var: list of labels, to identify the columns of data which 
#   #               correspond to the response variables
#   # covariates: list of columns to be used as covariates to build the tree
#   # categorical_cov: list of booleans to specify if each covariate is categorical or not
#   # logL_func: function to compute the log-likelihood than we want to maximize
#   # param_estimation: function to compute the parameter of interest
#   # minbucket: integer, the minimum nb. of obs. in any leaf (as in RPART)
#   #
#   # Output:
#   # predicted_tree: list of triangles, of the same length of the trained tree
#   
#   predicted_tree <- list()
#   nb_triangles <- length(trained_tree)
#   for (t in 1:nb_triangles){
#     
#     ### general info of the triangle
#     triangle_train <- trained_tree[[t]]
#     depth <- triangle_train$depth
#     label <- triangle_train$label
#     best_cov <- triangle_train$best_cov
#     best_th <- triangle_train$best_th
#     
#     ### get data to be splitted
#     if (depth==1){              # to create the first triangle, consider all the test set
#       data.P <- data_test
#       yyP <- as.matrix(data.P[,response_var], ncol=2)
#       param.P <- param_estimation(yyP)
#       logL.P <- logL_func(yyP)
#       n.P <- nrow(data.P)
#     } else {                    # otherwise consider the obs of the parent node (only if the parent node is a terminal/parent node)
#       label_parent <- substr(label, 1, nchar(label)-1)
#       brench_direction <- substr(label, nchar(label), nchar(label))
#       idx_parent <- which(sapply(predicted_tree, function(t){t$label})==label_parent)
#       status_parent <- predicted_tree[[idx_parent]]$status
#       if (status_parent=='terminal' | status_parent=='parent'){
#         if (brench_direction=='L'){
#           data.P <- predicted_tree[[idx_parent]]$data.L
#           n.P <- predicted_tree[[idx_parent]]$n.L
#           param.P <- predicted_tree[[idx_parent]]$param.L
#           logL.P <- predicted_tree[[idx_parent]]$logL.L
#         } else if (brench_direction=='R'){
#           data.P <- predicted_tree[[idx_parent]]$data.R
#           n.P <- predicted_tree[[idx_parent]]$n.R
#           param.P <- predicted_tree[[idx_parent]]$param.R
#           logL.P <- predicted_tree[[idx_parent]]$logL.R
#         }
#       } else {                  # if the parent triangle is 'non_convenient', stop the split procedure
#         data.P <- NULL
#       }
#     } 
#     
#     ### perform the split
#     gain <- 0
#     if (!is.null(data.P)){
#       if (grepl('aux', best_cov)){                        # if the covariate is categorical
#         cov <- substr(best_cov, 1, nchar(best_cov)-4)     # go back to the original label of the covariate
#         data.L <- data.P[data.P[,cov] %in% best_th,]
#         if (nrow(data.L)>=minbucket & nrow(data.L)<=(n.P-minbucket)){ # if the children nodes contain at least 'minbucket' obs, perform the split
#           data.R <- data.P[!data.P[,cov] %in% best_th,]
#           yyL <- as.matrix(data.L[,response_var], ncol=2)
#           yyR <- as.matrix(data.R[,response_var], ncol=2)
#           param.L <- param_estimation(yyL)
#           param.R <- param_estimation(yyR)
#           direction <- sign(sign(param.L-param.R)-0.5)   
#           logL.L <- logL_func(yyL)
#           logL.R <- logL_func(yyR)
#           n.L <- nrow(data.L)
#           n.R <- nrow(data.R)
#           gain <- logL.L + logL.R - logL_func(as.matrix(data.P[,response_var], ncol=2))
#         }
#       } else {                                            # if the covariate is quantitative
#         data.A <- data.P[data.P[,best_cov]<best_th,]
#         if (nrow(data.A)>=minbucket & nrow(data.A)<(n.P-minbucket)){ # if the children nodes contain at least minbucket obs, perform the split
#           data.B <- data.P[data.P[,best_cov]>=best_th,]
#           yyA <- as.matrix(data.A[,response_var], ncol=2)
#           yyB <- as.matrix(data.B[,response_var], ncol=2)
#           param.A <- param_estimation(yyA)
#           param.B <- param_estimation(yyB)
#           direction <- sign(sign(param.A-param.B)-0.5) # repeated use of the sign() function to avoid zeros
#           logL.A <- logL_func(yyA)
#           logL.B <- logL_func(yyB)
#           if (direction<0){
#             data.L <- data.A
#             data.R <- data.B
#             param.L <- param.A
#             param.R <- param.B
#             logL.L <- logL.A
#             logL.R <- logL.B
#           } else if (direction>0) {
#             data.L <- data.B
#             data.R <- data.A
#             param.L <- param.B
#             param.R <- param.A
#             logL.L <- logL.B
#             logL.R <- logL.A
#           }
#           n.L <- nrow(data.L)
#           n.R <- nrow(data.R)
#           gain <- logL.L + logL.R - logL_func(as.matrix(data.P[,response_var], ncol=2))
#         }
#       }
#     }
#     
#     ### define triangle and store it
#     if (gain>0){
#       triangle_test <- list(label=label,
#                             depth=depth,
#                             has_children=c(0,0),
#                             status='terminal',
#                             gain=gain,
#                             best_cov=best_cov,
#                             best_th=best_th,
#                             direction=direction,
#                             data.P=data.P,
#                             data.L=data.L,
#                             data.R=data.R,
#                             param.P=param.P,
#                             param.L=param.L,
#                             param.R=param.R,
#                             logL.P=logL.P,
#                             logL.L=logL.L,
#                             logL.R=logL.R,
#                             n.P=n.P,
#                             n.L=n.L,
#                             n.R=n.R)
#     } else {
#       triangle_test <- list(label=label,
#                             depth=depth,
#                             has_children=c(0,0),
#                             status='non_convenient',
#                             gain=0,
#                             best_cov=NULL,
#                             best_th=NULL,
#                             direction=NULL,
#                             data.P=NULL,
#                             data.L=NULL,
#                             data.R=NULL,
#                             param.P=NULL,
#                             param.L=NULL,
#                             param.R=NULL,
#                             logL.P=NULL,
#                             logL.L=NULL,
#                             logL.R=NULL,
#                             n.P=NULL,
#                             n.L=NULL,
#                             n.R=NULL)
#     }
#     
#     ### add triangle to the tree
#     predicted_tree[[t]] <- triangle_test
#     
#     ### update info of the parent triangle
#     if (gain>0 & depth>1){
#       predicted_tree[[idx_parent]]$status <- 'parent'
#       if (brench_direction=='L'){
#         predicted_tree[[idx_parent]]$has_children <- predicted_tree[[idx_parent]]$has_children + c(1,0)
#       } else if (brench_direction=='R'){
#         predicted_tree[[idx_parent]]$has_children <- predicted_tree[[idx_parent]]$has_children + c(0,1)
#       }
#     }
#   }
#   
#   return(predicted_tree)
# }
