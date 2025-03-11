#deal with survival time and censoring indicator
preprocessing <- function(hcc) {
  hcc$hss.day <- as.Date(hcc$Last_Date) - as.Date(hcc$OR_Date)  
  hcc$hss.day <- as.numeric(hcc$hss.day) 
  for (i in 1:nrow(hcc)) {
    if (is.na(hcc$Recur_Date[i])) {
      hcc$rfs.day[i] <- as.Date(hcc$Last_Date[i]) - as.Date(hcc$OR_Date[i])
    } else {
      hcc$rfs.day[i] <- as.Date(hcc$Recur_Date[i]) - as.Date(hcc$OR_Date[i])
    }
  }
  
  hcc$hss.censor <- 1
  hcc$hss.censor[which(is.na(hcc$Death_Date))] <- 0
  
  hcc$rfs.censor <- 1
  hcc$rfs.censor[which(is.na(hcc$Recur_Date))] <- 0 
  hcc$rfs.censor[which(!is.na(hcc$Death_Date))] <- 1
  return(hcc)
}


#calculate the expected survival time E[T|X,A]
expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))
  return(c(cbind(1, S.hat) %*% grid.diff))
}


#augmented-owl
owl.md <- function(X, R, A, testX) {
  
  # Function for Implement the proposed method OWL-MD
  # Input:
  # X: design matrix (do not include intercept)
  # R: outcome vector
  # A: treatment vector
  # testX: test data matrix (do not include intercept)
  
  # Output:
  # trt.owlmd: predict optimal treatments on the test data
  
  X <- as.matrix(X)  # do not include intercept
  X.total <- cbind(1, X)  # add constant (intercept)
  n <- length(R)  # sample size
  p <- ncol(X)  # covariate dimension
  K <- length(unique(A))  # treatment number
  trt <- 1:K  # treatment space {1,...,K}
  # estimate the Q-functions
  Q.hat <- matrix(data = 0, nrow = n, ncol = K)
  reg <- lm(formula = R ~ X * factor(A))  # fit regression model
  for (k in 1:K) {
    Q.hat[, k] <- predict(object = reg, newdata = data.frame(X, A = rep(trt[k], n)))
  }
  # estimate the sample propensity score
  #pr <- as.numeric(table(A)/n)
  #Ps <- matrix(data = pr, nrow = n, ncol = K, byrow = TRUE)
  fit <- multinom(A ~ X, data = data.frame(X, A), trace = F)
  Ps <- fitted(fit)
  
  # estimate weights matrix
  W <- matrix(data = NA, nrow = n, ncol = K)
  for (k in 1:K) {
    W[, k] <- (A == k) * R/Ps[, k] + (1 - (A == k)/Ps[, k]) * Q.hat[, k]
  }
  # row summation of weight matrix
  v <- rowSums(W)
  # row summation of negative part of weight matrix
  u <- rowSums(abs(W) * (W < 0))
  
  # initial value of 
  beta.old <- rep(0, (p + 1) * (K - 1))
  for (i in 1:100) {
    # the (p+1)*(K-1) matrix
    Beta <- matrix(data = beta.old, nrow = p + 1, ncol = K - 1, byrow = FALSE)
    # linear part in logit transformation
    Eta <- X.total %*% cbind(Beta, 0)
    expEta <- exp(Eta)
    # calculation sample probabilities
    Prob <- expEta/rowSums(expEta)  #pk
    
    # gradient vector
    ng <- Reduce(f = "+", x = sapply(X = 1:n, FUN = function(i) {
      kronecker(X = as.matrix(W[i, -K] + u[i] - (v[i] + K * u[i]) * Prob[i, -K]), Y = as.matrix(X.total[i, ]))
    }, simplify = FALSE))
    # Hessian matrix
    rowPr <- diag(K - 1)
    nH <- -Reduce(f = "+", x = sapply(X = 1:n, FUN = function(i) {
      diag(rowPr) <- Prob[i, -K]
      (v[i] + K * u[i]) * kronecker(X = (rowPr - tcrossprod(Prob[i, -K])), Y = tcrossprod(X.total[i, ]))
    }, simplify = FALSE))
    # update beta
    beta.new <- try(expr = beta.old - solve(nH, ng), silent = TRUE)
    if ("try-error" %in% class(beta.new)) {
      beta.new <- rep(NA, (p + 1) * (K - 1))
      cat("Netwon mathed failed.", "\n")
      break
    }
    # stop rule
    beta.new <- as.numeric(beta.new)
    if (sum(is.na(beta.new)))
      break  # divergence
    if (sum(abs(beta.new - beta.old)) < 10^(-7))
      break  # convergence
    beta.old <- beta.new
  }
  # coefficient matrix
  Beta <- matrix(data = beta.old, nrow = p + 1, ncol = K - 1, byrow = FALSE)
  Beta <- cbind(Beta, 0)  #fK=0
  colnames(Beta) <- paste("beta", 1:K, sep = "")
  # predict treatments
  trt.owlmd <- apply(X = cbind(1, testX) %*% Beta, MARGIN = 1, FUN = which.max)
  
  # return predicted treatments
  return(as.numeric(trt.owlmd))
}

#Matching based on covariates 
MultiMatch <- function(reference_group, simu_GPS, calipernum = NULL) {
  dim_simuGPS <- dim(simu_GPS)[2]
  if (reference_group == 1) {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 1"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 3"
    trueT1 <- 1
    trueT2 <- 2
    trueT3 <- 3
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  } else if (reference_group == 2) {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 1"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 3"
    trueT1 <- 2
    trueT2 <- 1
    trueT3 <- 3
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  } else if (reference_group == 3) {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 3"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 1"
    trueT1 <- 3
    trueT2 <- 2
    trueT3 <- 1
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  }
  
  # T1 reference treatment
  temp12 <- dplyr::filter(simu_GPS, treat != paste("Treatment", trueT3))
  temp13 <- dplyr::filter(simu_GPS, treat != paste("Treatment", trueT2))
  
  # matching T1 with the remaining groups
  match12 <- Match(Y = NULL, Tr = temp12$treat == paste("Treatment", trueT1), X = temp12[, 1:p], caliper = calipernum, ties = FALSE)
  
  match13 <- Match(Y = NULL, Tr = temp13$treat == paste("Treatment", trueT1), X = temp13[, 1:p], caliper = calipernum, ties = FALSE)
  
  simu_GPS$id <- 1:nrow(simu_GPS)
  
  #units in T1 that matched to all the remaining groups
  simu_GPS$both.1 <- simu_GPS$id %in% match12$index.treated & simu_GPS$id %in% match13$index.treated
  
  temp <- simu_GPS[simu_GPS$both.1 == "TRUE", ]
  m12 <- cbind(match12$index.treated, match12$index.control)
  m13 <- cbind(match13$index.treated, match13$index.control + sum(simu_GPS$T2 == T))
  
  m12 <- m12[m12[, 1] %in% rownames(temp), ]
  m13 <- m13[m13[, 1] %in% rownames(temp), ]
  
  matchset <- cbind(m12[order(m12[, 1]), ], m13[order(m13[, 1]), ])
  matchset <- as.matrix(matchset[, c(1, 2, 4)])  #Matched set
  n.trip <- nrow(matchset)
  
  ############ Deal the matching result #########
  #In the same match set, Bi denote argmax Ri, R_Bi denote maxRi, Ci denote argmax Ri, R_Ci denote minRi, g_weight denote g1() in paper 
  simu_GPS <- cbind(simu_GPS[, 1:dim_simuGPS], B = rep(0, dim(simu_GPS)[1]), R_Bi = rep(0, dim(simu_GPS)[1]), R_Ci = rep(0, dim(simu_GPS)[1]), g_weight = rep(0, dim(simu_GPS)[1]))
  for (i in 1:n.trip) {
    index <- matchset[i, ]
    temp <- data.frame(index = c(index[1], index[2], index[3]), impute = c(simu_GPS$impute[index[1]], simu_GPS$impute[index[2]], simu_GPS$impute[index[3]]), trt = c(trueT1, trueT2, trueT3))
    maxvalue <- apply(temp, 2, max)[2]
    minvalue <- apply(temp, 2, min)[2]
    argmax <- temp[which(temp$impute == maxvalue), 3]
    argmin <- temp[which(temp$impute == minvalue), 3]
    if (length(argmax) != 1) {
      argmax <- min(argmax)
    }
    if (length(argmin) != 1) {
      argmin <- min(argmin)
    }
    simu_GPS$B[index[1]] <- argmax
    simu_GPS$R_Bi[index[1]] <- maxvalue
    simu_GPS$R_Ci[index[1]] <- minvalue
    
    simu_GPS$g_weight[index[1]] <- abs(maxvalue - temp$impute[1]) + abs(maxvalue - temp$impute[2]) + abs(maxvalue - temp$impute[3])
    
  }
  
  Match.result <- simu_GPS[which(simu_GPS$B != 0), ]
  
  #matched set
  colnames(matchset) <- c(paste("Treatment", trueT1), paste("Treatment", trueT2), paste("Treatment", trueT3))
  #adjust the order
  for (i in 1:K) {
    matchset[, i] <- simu_GPS$origin_order[matchset[, i]]
  }
  
  return <- list(`reference group` = reference_group, result = Match.result, `match set` = matchset)
}

#Matching based on estimated generalized propensity score
MultiMatchGps <- function(reference_group, simu_GPS, calipernum = NULL) {
  dim_simuGPS <- dim(simu_GPS)[2]
  if (reference_group == 1) {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 1"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 3"
    trueT1 <- 1
    trueT2 <- 2
    trueT3 <- 3
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  } else if (reference_group == 2) {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 1"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 3"
    trueT1 <- 2
    trueT2 <- 1
    trueT3 <- 3
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  } else if (reference_group == 3) {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 3"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 1"
    trueT1 <- 3
    trueT2 <- 2
    trueT3 <- 1
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  }
  
  # T1 reference treatment
  temp12 <- dplyr::filter(simu_GPS, treat != paste("Treatment", trueT3))
  temp13 <- dplyr::filter(simu_GPS, treat != paste("Treatment", trueT2))
  
  # matching T1 with the remaining groups
  match12 <- Match(Y = NULL, Tr = temp12$treat == paste("Treatment", trueT1), X = temp12[, 1:3], caliper = calipernum, ties = FALSE)
  
  match13 <- Match(Y = NULL, Tr = temp13$treat == paste("Treatment", trueT1), X = temp13[, 1:3], caliper = calipernum, ties = FALSE)
  
  simu_GPS$id <- 1:nrow(simu_GPS)
  
  #units in T1 that matched to all the remaining groups
  simu_GPS$both.1 <- simu_GPS$id %in% match12$index.treated & simu_GPS$id %in% match13$index.treated
  
  temp <- simu_GPS[simu_GPS$both.1 == "TRUE", ]
  m12 <- cbind(match12$index.treated, match12$index.control)
  m13 <- cbind(match13$index.treated, match13$index.control + sum(simu_GPS$T2 == T))
  
  m12 <- m12[m12[, 1] %in% rownames(temp), ]
  m13 <- m13[m13[, 1] %in% rownames(temp), ]
  
  matchset <- cbind(m12[order(m12[, 1]), ], m13[order(m13[, 1]), ])
  matchset <- as.matrix(matchset[, c(1, 2, 4)])  #Matched set
  n.trip <- nrow(matchset)
  
  ############ Deal the matching result#########
  #In the same match set, Bi denote argmax Ri, R_Bi denote maxRi, Ci denote argmax Ri, R_Ci denote minRi, g_weight denote g1() in paper 
  simu_GPS <- cbind(simu_GPS[, 1:dim_simuGPS], B = rep(0, dim(simu_GPS)[1]), R_Bi = rep(0, dim(simu_GPS)[1]), R_Ci = rep(0, dim(simu_GPS)[1]), g_weight = rep(0, dim(simu_GPS)[1]))
  for (i in 1:n.trip) {
    index <- matchset[i, ]
    temp <- data.frame(index = c(index[1], index[2], index[3]), impute = c(simu_GPS$impute[index[1]], simu_GPS$impute[index[2]], simu_GPS$impute[index[3]]), trt = c(trueT1, trueT2, trueT3))
    maxvalue <- apply(temp, 2, max)[2]
    minvalue <- apply(temp, 2, min)[2]
    argmax <- temp[which(temp$impute == maxvalue), 3]
    argmin <- temp[which(temp$impute == minvalue), 3]
    if (length(argmax) != 1) {
      argmax <- min(argmax)
    }
    if (length(argmin) != 1) {
      argmin <- min(argmin)
    }
    simu_GPS$B[index[1]] <- argmax
    simu_GPS$R_Bi[index[1]] <- maxvalue
    simu_GPS$R_Ci[index[1]] <- minvalue
    if (argmax == reference_group) {
      simu_GPS$R_Zi[index[1]] <- minvalue
    } else {
      simu_GPS$R_Zi[index[1]] <- simu_GPS$impute[index[1]]
    }
    
    simu_GPS$g_weight[index[1]] <- abs(maxvalue - temp$impute[1]) + abs(maxvalue - temp$impute[2]) + abs(maxvalue - temp$impute[3])
    
  }
  
  Match.result <- simu_GPS[which(simu_GPS$B != 0), ]
  
  colnames(matchset) <- c(paste("Treatment", trueT1), paste("Treatment", trueT2), paste("Treatment", trueT3))
  #matched set 
  for (i in 1:K) {
    matchset[, i] <- simu_GPS$origin_order[matchset[, i]]
  }
  
  return <- list(`reference group` = reference_group, result = Match.result, `match set` = matchset)
}


#misclassification rate
error.rate <- function(y, fit.class) return(sum(fit.class != y)/length(fit.class))

#cross validation for tuning penalty parameter lambda
cvfun <- function(inputX, inputY, originR, originA, fold = 3, lambda, kernel = "linear", kparam = 1, weight = NA) {
  #X covariate, Y label for MSVM, A original treatment, R: outcome
  set.seed(2021)
  if (kernel == "linear") {
    kparam <- 1
  }
  
  folds <- createFolds(inputY, k = fold)
  ValueFun <- matrix(0, ncol = length(kparam), nrow = length(lambda))
  for (i in 1:fold) {
    cat("Leaving subset[", i, "] out in", fold, "fold CV:", "\n")
    fold_testx <- inputX[folds[[i]], ]  #folds[[i]] for validation
    fold_testy <- inputY[folds[[i]], ]
    fold_trainx <- inputX[-folds[[i]], ]  #the remaining for training
    fold_trainy <- inputY[-folds[[i]], ]
    
    # A and R used for valuefun
    fold_testR <- originR[folds[[i]], ]
    fold_testA <- originA[folds[[i]], ]
    fold_testing <- data.frame(cbind(A = fold_testA, R = fold_testR, fold_testx))
    
    fold_weight <- if (sum(is.na(weight)) != 0) {
      NA
    } else {
      weight[-folds[[i]], ]
    }
    
    row.index <- 0
    for (j in lambda) {
      col.index <- 0
      row.index <- row.index + 1
      for (k in kparam) {
        col.index <- col.index + 1
        if (sum(is.na(fold_weight)) != 0) {
          ramsvm.out <- try(ramsvm(fold_trainx, fold_trainy, lambda = j, kernel = kernel, kparam = k), TRUE)
          if (is.character(ramsvm.out)) {
            print("cv ramsvm failed")
            next
          }
          fit.class <- predict(ramsvm.out, fold_testx)
          j <- paste(j)
          fit.class <- fit.class[[j]]
        } else {
          ramsvm.out <- try(ramsvm(fold_trainx, fold_trainy, lambda = j, kparam = k, kernel = kernel, weight = as.vector(fold_weight)), TRUE)
          if (is.character(ramsvm.out)) {
            print("cv ramsvm failed")
            next
          }
          fit.class <- predict(ramsvm.out, fold_testx)
          j <- paste(j)
          fit.class <- fit.class[[j]]
        }
        ValueFun[row.index, col.index] <- ValueFun[row.index, col.index] + valuefun(fold_testing, est_ITR = fit.class)/fold
      }
      cat("*")
    }
    cat("\n")
  }
  optIndex <- which(ValueFun == max(ValueFun), arr.ind = TRUE)[1, ]
  return(list(lambda = lambda[optIndex[1]], kparam = kparam[optIndex[2]], value = ValueFun))
}

g <- function(x, y) {
  return(x - y)
}

propensity <- function(testing, A) {
  PS <- switch(as.numeric(A), testing$p1, testing$p2, testing$p3)
  return(PS)
}

#calculate value function
valuefun <- function(testing, est_ITR, if_test = FALSE, testing_id = 1) {
  # Input: covariate, treatment, outcome, estimated ITR
  # Output: value function
  # When if_test=FLASE, estimate the propensity score by multinomial logistics
  # When calculate the testing performance, set if_test=TRUE and use the estimated generalized propensity score by all sample
  # value function = sum(I[Ai=D(Xi)]Ri/p(Ai,Xi)) / sum(I[Ai=D(Xi)] /p(Ai,Xi))
  # Since in real data we don't know the true survival time, we use the empirical pseudo value function that replace censored observations by imputation
  id <- which(est_ITR == testing$A)
  if (if_test == FALSE) {
    logit_formula <- paste("A  ~", paste(covariate_name, collapse = "+"))
    
    fit <- multinom(logit_formula, data = testing, trace = F)
    Rx <- fitted(fit)
    colnames(Rx) <- c("p1", "p2", "p3")
    test_PS <- cbind(testing, Rx)
    
    denominator <- numerator <- 0
    for (i in id) denominator <- denominator + 1/propensity(test_PS[i, ], test_PS$A[i])
    for (i in id) numerator <- numerator + test_PS$R[i]/propensity(test_PS[i, ], test_PS$A[i])
  } else {
    R1 <- testing$R
    all_sample_PS <- all_sample_PS[testing_id, ]
    denominator <- numerator <- 0
    for (i in id) denominator <- denominator + 1/propensity(all_sample_PS[i, ], all_sample_PS$treat[i])
    for (i in id) numerator <- numerator + R1[i]/propensity(all_sample_PS[i, ], all_sample_PS$treat[i])
  }
  return(numerator/denominator)
}




########### all the comparing methods ##############
paraPredict <- function(ii) {
  set.seed(ii)
  
  result <- matrix(0, fold, 10)
  colnames(result) <- c("g1-cov", "gweight1-cov", "gweight2-cov", "g1-gps", "gweight1-gps", "gweight2-gps", "aug-owl", "owl", "cox", "weight-cox")
  
  performance <- matrix(0, 1, 10)
  rownames(performance) <- "value fun"
  colnames(performance) <- c("g1-cov", "gweight1-cov", "gweight2-cov", "g1-gps", "gweight1-gps", "gweight2-gps", "aug-owl", "owl", "cox", "weight-cox")
  
  folds <- createFolds(hss_impute$impute, k = fold)
  
  for (f in 1:fold) {
    cat("Overall Leaving subset[", f, "] out in", fold, "fold CV:", "\n")
    
    ####### partition training and testing ########
    train <- hss_impute[-folds[[f]], ]
    testing <- hss_impute[folds[[f]], ]
    n_train <- nrow(train)
    
    testing <- cbind(testing[, 1:p], testing$treat, testing$censor, testing$impute)
    colnames(testing)[(p + 1):length(testing)] <- c("A", "censor", "R")
    
    ###############Standard Cox model#####################
    Coxtrain <- cbind(train[, 1:p], train$treat, train$censor, train$truncateR)
    colnames(Coxtrain)[(p + 1):length(Coxtrain)] <- c("A", "censor", "R")
    
    Coxtrain$A <- factor(Coxtrain$A)
    testing$A <- factor(testing$A)
    
    #Cox regression with (X,A,AX)
    cox_formula <- formula(paste("Surv(R,censor) ~ A +", paste(covariate_name, collapse = "+"), "+", paste(covariate_name, rep("*A", length(covariate_name)), sep = "", collapse = "+")))
    cox_model <- coxph(cox_formula, data = Coxtrain, method = "breslow", eps = 1e-07, iter.max = 20)
    
    testA1 <- data.frame(testing[, 1:p], A = factor(rep(1, nrow(testing))))
    testA2 <- data.frame(testing[, 1:p], A = factor(rep(2, nrow(testing))))
    testA3 <- data.frame(testing[, 1:p], A = factor(rep(3, nrow(testing))))
    
    cox_result <- matrix(rep(0, nrow(testing) * 3), ncol = 3)
    dimnames(cox_result) <- list(1:nrow(testing), 1:3)
    cox_result[, 1] <- predict(cox_model, testA1, type = "risk")
    cox_result[, 2] <- predict(cox_model, testA2, type = "risk")
    cox_result[, 3] <- predict(cox_model, testA3, type = "risk")
    
    predict_ITR <- as.numeric(apply(cox_result, 1, function(t) colnames(cox_result)[which.min(t)]))
    
    performance[1, 9] <- valuefun(testing, est_ITR = predict_ITR, if_test = TRUE, testing_id = folds[[f]])
    testing$cox <- predict_ITR
    
    cat("standard cox over!", "\n")
    
    ###############Weighted Cox model#####################
    #not displayed in the paper
    fit <- multinom(A ~ . - R - censor, data = Coxtrain, trace = F)
    Rx <- fitted(fit)
    colnames(Rx) <- c("p1", "p2", "p3")
    Coxtrain_GPS <- cbind(Coxtrain, Rx)
    
    for (i in 1:nrow(Coxtrain_GPS)) {
      if (Coxtrain_GPS[i, ]$A == "1")
        Coxtrain_GPS$weight[i] <- Coxtrain_GPS$p1[i]
      if (Coxtrain_GPS[i, ]$A == "2")
        Coxtrain_GPS$weight[i] <- Coxtrain_GPS$p2[i]
      if (Coxtrain_GPS[i, ]$A == "3")
        Coxtrain_GPS$weight[i] <- Coxtrain_GPS$p3[i]
    }
    
    cox_formula <- formula(paste("Surv(R,censor) ~ A +", paste(covariate_name, collapse = "+"), "+", paste(covariate_name, rep("*A", length(covariate_name)), sep = "", collapse = "+")))
    cox_GPS_model <- coxph(cox_formula, data = Coxtrain_GPS, weights = weight)
    
    cox_result <- matrix(rep(0, nrow(testing) * 3), ncol = 3)
    dimnames(cox_result) <- list(1:nrow(testing), 1:3)
    cox_result[, 1] <- predict(cox_GPS_model, testA1, type = "risk")
    cox_result[, 2] <- predict(cox_GPS_model, testA2, type = "risk")
    cox_result[, 3] <- predict(cox_GPS_model, testA3, type = "risk")
    
    predict_ITR <- as.numeric(apply(cox_result, 1, function(t) colnames(cox_result)[which.min(t)]))
    
    performance[1, 10] <- valuefun(testing, est_ITR = predict_ITR, if_test = TRUE, testing_id = folds[[f]])
    testing$weight_cox <- predict_ITR
    cat("weighted cox owl over!", "\n")
    
    ############## augmented owl ################
    Owltrain <- cbind(train[, 1:p], train$treat, train$censor, train$impute)
    colnames(Owltrain)[(p + 1):length(Owltrain)] <- c("A", "censor", "R")
    
    trt.owlmd <- owl.md(X = as.matrix(Owltrain[, 1:p]), R = Owltrain$R, A = Owltrain$A, testX = as.matrix(testing[, 1:p]))
    performance[1, 7] <- valuefun(testing, est_ITR = trt.owlmd, if_test = TRUE, testing_id = folds[[f]])
    testing$aug_owl <- trt.owlmd
    
    cat("augmented owl over!", "\n")
    
    
    ################## Match+ramsvm based methods#################
    
    ################## CovMatch #############################
    train_GPS <- cbind(train, origin_order = 1:nrow(train))
    
    train_GPS[which(train_GPS$treat == 1), ]$treat <- "Treatment 1"
    train_GPS[which(train_GPS$treat == 2), ]$treat <- "Treatment 2"
    train_GPS[which(train_GPS$treat == 3), ]$treat <- "Treatment 3"
    
    #Matching based on covariates
    MultiMtrt1 <- MultiMatch(1, train_GPS, calipernum)
    MultiMtrt2 <- MultiMatch(2, train_GPS, calipernum)
    MultiMtrt3 <- MultiMatch(3, train_GPS, calipernum)
    cat("Match over!", "\n")
    
    #for calculate valuefun
    originR <- rbind(as.matrix(train$impute[MultiMtrt1$`match set`[, 1]]), as.matrix(train$impute[MultiMtrt2$`match set`[, 1]]), as.matrix(train$impute[MultiMtrt3$`match set`[, 1]]))
    
    originA <- as.matrix(rep(c(1, 2, 3), c(length(MultiMtrt1$result$impute), length(MultiMtrt2$result$impute), length(MultiMtrt3$result$impute))))
    
    #########################gaussian g1##########################################
    inputX <- rbind(MultiMtrt1$result[, 1:p], MultiMtrt2$result[, 1:p], MultiMtrt3$result[, 1:p])
    inputX <- as.matrix(inputX)
    inputY <- rbind(as.matrix(MultiMtrt1$result$B), as.matrix(MultiMtrt2$result$B), as.matrix(MultiMtrt3$result$B))
    
    cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian")
    ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian"), TRUE)
    if (is.character(ramsvm.out)) {
      return("gaussian g1 failed")
    }
    
    #prediction on test data
    fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
    
    performance[1, 1] <- valuefun(testing, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
    testing$g1_cov <- fit.class[[paste(cv_param$lambda)]]
    
    cat("gaussian g1 over!", "\n")
    
    
    #########################gaussian gweight2##########################################
    #weight
    msvm_weight <- as.matrix(c(g(MultiMtrt1$result$R_Bi, MultiMtrt1$result$R_Ci), g(MultiMtrt2$result$R_Bi, MultiMtrt2$result$R_Ci), g(MultiMtrt3$result$R_Bi, MultiMtrt3$result$R_Ci)))
    
    cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
    ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
    if (is.character(ramsvm.out)) {
      return("gaussian gweight2 failed")
    }
    #prediction on test data
    fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
    
    performance[1, 3] <- valuefun(testing, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
    testing$gx_v2_cov <- fit.class[[paste(cv_param$lambda)]]
    
    cat("gaussian gweight2 over!", "\n")
    
    ##########################gaussian gweight1##########################################
    #train
    msvm_weight <- as.matrix(c(MultiMtrt1$result$g_weight, MultiMtrt2$result$g_weight, MultiMtrt3$result$g_weight))
    cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
    ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
    if (is.character(ramsvm.out)) {
      return("gaussian gweight1 failed")
    }
    #prediction on test data
    fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
    
    performance[1, 2] <- valuefun(testing, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
    testing$gweight_cov <- fit.class[[paste(cv_param$lambda)]]
    
    cat("gaussian gweight1!", "\n")
    
    
    ########################GpsMatch################################################
    
    p_MatchGps <- dim(hss_impute)[2]
    
    #Estimated generalize propensity scores for matching in training data, but still using covariates for fitting ITR
    #train_MatchGps includes gps1 gps2 gps3 treatment censor truncateR imputeR
    train_MatchGps <- all_sample_PS[-folds[[f]], c(p_MatchGps + 1, p_MatchGps + 2, p_MatchGps + 3, p_MatchGps - 3, p_MatchGps - 2, p_MatchGps - 1, p_MatchGps)]
    train_GPS <- cbind(train_MatchGps, origin_order = 1:nrow(train_MatchGps))
    
    train_GPS[which(train_GPS$treat == 1), ]$treat <- "Treatment 1"
    train_GPS[which(train_GPS$treat == 2), ]$treat <- "Treatment 2"
    train_GPS[which(train_GPS$treat == 3), ]$treat <- "Treatment 3"
    
    #Matching based on gps
    MultiMtrt1 <- MultiMatchGps(1, train_GPS, calipernum)
    MultiMtrt2 <- MultiMatchGps(2, train_GPS, calipernum)
    MultiMtrt3 <- MultiMatchGps(3, train_GPS, calipernum)
    cat("Match over!", "\n")
    
    
    originR <- rbind(as.matrix(train$impute[MultiMtrt1$`match set`[, 1]]), as.matrix(train$impute[MultiMtrt2$`match set`[, 1]]), as.matrix(train$impute[MultiMtrt3$`match set`[, 1]]))
    
    originA <- as.matrix(rep(c(1, 2, 3), c(length(MultiMtrt1$result$impute), length(MultiMtrt2$result$impute), length(MultiMtrt3$result$impute))))
    
    #########################gaussian g1##########################################
    #train 
    inputX <- rbind(train[MultiMtrt1$`match set`[, 1], 1:p], train[MultiMtrt2$`match set`[, 1], 1:p], train[MultiMtrt3$`match set`[, 1], 1:p])
    inputX <- as.matrix(inputX)
    inputY <- rbind(as.matrix(MultiMtrt1$result$B), as.matrix(MultiMtrt2$result$B), as.matrix(MultiMtrt3$result$B))
    
    cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian")
    ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian"), TRUE)
    if (is.character(ramsvm.out)) {
      return("gaussian g1 failed")
    }
    
    #prediction on test data
    fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
    
    performance[1, 4] <- valuefun(testing, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
    testing$g1_gps <- fit.class[[paste(cv_param$lambda)]]
    
    cat("gaussian g1 over!", "\n")
    
    
    #########################gaussian gweight2##########################################
    #train
    msvm_weight <- as.matrix(c(g(MultiMtrt1$result$R_Bi, MultiMtrt1$result$R_Ci), g(MultiMtrt2$result$R_Bi, MultiMtrt2$result$R_Ci), g(MultiMtrt3$result$R_Bi, MultiMtrt3$result$R_Ci)))
    
    cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
    ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
    if (is.character(ramsvm.out)) {
      return("gaussian gweight2 failed")
    }
    #prediction on test data
    fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
    
    performance[1, 6] <- valuefun(testing, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
    testing$gx_v2_gps <- fit.class[[paste(cv_param$lambda)]]
    
    cat("gaussian gweight2 over!", "\n")
    
    ##########################gaussian gweight1##########################################
    #train
    msvm_weight <- as.matrix(c(MultiMtrt1$result$g_weight, MultiMtrt2$result$g_weight, MultiMtrt3$result$g_weight))
    cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
    ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
    if (is.character(ramsvm.out)) {
      return("gaussian gweight1 failed")
    }
    #prediction on test data
    fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
    
    performance[1, 5] <- valuefun(testing, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
    testing$gweight_gps <- fit.class[[paste(cv_param$lambda)]]
    
    cat("gaussian gweight1 over!", "\n")
    
    
    ###############gaussian owl########################
    #train
    inputX <- as.matrix(train[, 1:p])
    inputY <- as.matrix(train$treat)
    
    train_df <- data.frame(inputX, A = inputY, R = train$impute)
    fit <- multinom(A ~ . - R, data = train_df, trace = F)
    Rx <- fitted(fit)
    colnames(Rx) <- c("p1", "p2", "p3")
    train_df_GPS <- cbind(train_df, Rx)
    
    msvm_weight <- rep(0, n_train)
    
    #weighted by Yi/πi
    msvm_weight[1:(n_train/K)] <- train_df_GPS$R[1:(n_train/K)]/train_df_GPS$p1[1:(n_train/K)]
    msvm_weight[(n_train/K + 1):(2 * n_train/K)] <- train_df_GPS$R[(n_train/K + 1):(2 * n_train/K)]/train_df_GPS$p2[(n_train/K + 1):(2 * n_train/K)]
    msvm_weight[(2 * n_train/K + 1):(3 * n_train/K)] <- train_df_GPS$R[(2 * n_train/K + 1):(3 * n_train/K)]/train_df_GPS$p3[(2 * n_train/K + 1):(3 * n_train/K)]
    
    msvm_weight <- as.matrix(msvm_weight)
    
    cv_param <- cvfun(inputX, inputY, originR = as.matrix(train_df$R), originA = inputY, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
    ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = msvm_weight), TRUE)
    if (is.character(ramsvm.out)) {
      return("gaussian owl failed")
    }
    
    #prediction on test data
    fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
    
    performance[1, 8] <- valuefun(testing, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
    testing$gau_owl <- fit.class[[paste(cv_param$lambda)]]
    
    cat("gaussian owl over!", "\n")
    
    result[f, ] <- performance  
  }

  return(result)
}


