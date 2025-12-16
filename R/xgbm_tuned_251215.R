################################################################################
##### xgbm_tuned_yymmdd.R ######################################################
################################################################################
## SIMPLE XGBoost fit ==========================================================
#' Get a simple XGBoost model fit (no tuning) 
#' 
#' @description This fits a gradient boosting machine model using the XGBoost
#' platform.  If uses a single set of hyperparameters that have sometimes been 
#' reasonable so runs very fast.  For a better fit one can use xgb.tuned()
#' which searches for a set of hyperparameters using the mlrMBO package
#' which will generally provide a better fit but take much longer.  See xgb.tuned()
#' for a description of the data format required for input.  
#' 
#' @param train.xgb.dat The data to be used for training the XGBoost model
#' @param booster for now just "gbtree" (default) 
#' @param objective one of "survival:cox" (default), "binary:logistic" or "reg:squarederror"
#' @param eval_metric one of "cox-nloglik" (default), "auc", "rmse" or NULL.  Default
#' of NULL will select an appropriate value based upon the objective value.    
#' @param minimize whether the eval_metric is to be minimized or maximized
#' @param seed a seed for set.seed() to assure one can get the same results twice.  If NULL 
#' the program will generate a random seed.  Whether specified or NULL, the seed is stored in the output
#' object for future reference.  
#' @param folds an optional list where each element is a vector of indexes for a 
#' test fold.  Default is NULL.  If specified then doxgb$nfold is ignored as in xgb.cv().
#' @param doxgb a list with parameters for passed to xgb.cv() including $nfold, $nrounds,
#' and $early_stopping_rounds.  If not provided defaults will be used.  Defaults
#' can be seen form the output object$doxgb element, again a list. In case not NULL, 
#' the seed and folds option values override the $seed and $folds values in doxgb.   
#' @param track 0 (default) to not track progress, 2 to track progress.   
#' 
#' @return a XGBoost model fit 
#' 
#' @seealso
#   \code{\link{predict_nested_xgb}} , 
#'   \code{\link{xgb.tuned}} , \code{\link{nested.glmnetr}} 
#' 
#' @author Walter K Kremers with contributions from Nicholas B Larson
#'
#' @importFrom xgboost xgb.cv xgb.train xgb.DMatrix getinfo 
#' 
#' @export
#' 
xgb.simple = function(train.xgb.dat, 
                      booster     = "gbtree",
                      objective   = "survival:cox",
                      eval_metric = NULL,
                      minimize = NULL,
                      seed = NULL, 
                      folds=NULL, 
                      doxgb = NULL, 
                      track = 2 ) {
  
  if (is.null(doxgb)) { doxgb = list() }
  
  if  (!is.null(seed)) { doxgb$seed = seed   
  } else if (is.null( doxgb$seed )) { seed = round(runif(1)*1e9)  ; doxgb$seed = seed }
  if (seed == 0) { seed = round(runif(1)*1e9) ; doxgb$seed = seed }
  
  #  if (!is.null(folds)) { doxgb$folds = folds } 
  
  if ( (is.null(folds)) & (!is.null(doxgb$folds)) ) { folds = doxgb$folds } 
  
  if (!is.null(folds)) { doxgb$folds = folds ; doxgb$nfold = length(doxgb$folds) } 
  
  if ((is.null(doxgb$folds)) & (is.null(doxgb$nfold)))  { doxgb$nfold = 5 } 
  
  if ((is.null(doxgb$folds)) & (!is.null(doxgb$nfold)))  { doxgb$nfold = max(doxgb$nfold, 3) } 
  
  if (is.null( doxgb$nrounds )) { doxgb$nrounds = 1000 } 
  if (is.null( doxgb$early_stopping_rounds )) { doxgb$early_stopping_rounds = min(100,max(10,doxgb$nrounds/5)) } 
  
  if (objective == "survival:cox") { family = "cox" 
  } else if (objective == "binary:logistic") { family = "binomial" 
  } else { family = "gaussian" }
  
  if (track >= 3) { print("XGB Simple - about to get folds")}
  
  stratified = 1 
  if (is.null(folds)) { 
    if  (is.null(seed))   { seed = round(runif(1)*1e9) 
    } else if (seed == 0) { seed = round(runif(1)*1e9) 
    }
    set.seed(seed) 
    label = getinfo(train.xgb.dat,'label')
    if (objective == "survival:cox") {
      y_ = abs( label )
      event = ( label > 0) * 1 
    } else {
      y_ = label  
      event = NULL
    }
    folds_v = get.foldid(y_, event, family, doxgb$nfold, stratified) 
    folds = list()
    for (i1 in c(1:doxgb$nfold)) { folds[[i1]] = c(1:length(y_))[folds_v == i1] } 
  }
  
  if (track >= 3) { print("XGB Simple - about to get eval_metric")}
  
  if        (objective == "survival:cox"    ) { minimize = TRUE ; if (is.null(eval_metric)) { eval_metric = "cox-nloglik" }
  } else if (objective == "binary:logistic" ) { if (is.null(minimize)) { minimize = FALSE } ; if (is.null(eval_metric)) { eval_metric = "auc" }
  } else if (objective == "reg:squarederror") { minimize = TRUE ;  if (is.null(eval_metric)) { eval_metric = "rmse" }
  } else if (objective == "count:poisson"   ) { minimize = TRUE ;  if (is.null(eval_metric)) { eval_metric = "poisson-nloglik" }
  } else { minimize = TRUE }
  
  if (track >= 3) { print("XGB Simple - set param.final") }
  
  param.final = list(booster = booster, objective = objective, eval_metric = eval_metric,
                     nthread = 2, eta = 1, max_depth = 2, nrounds = 2 )
  param.final = list(booster = booster, objective = objective, eval_metric = eval_metric,
                     nthread = 16, eta = .2, max_depth = 10, gamma=0.1, min_child_weight=5, 
                     colsample_bytree=0.75, alpha=0.1, subsample=0.6 )
  
  if (track >= 3) { print("XGB Simple - set early stopping rounds") }
  
  if (track >= 3) { print("XGB Simple - get xgb.cv.fit") }

  set.seed(seed)                                                                ## 240315 
  
  xgb.cv.fit = xgb.cv(param.final,
                      data = train.xgb.dat,
                      nrounds = doxgb$nrounds,
                      showsd = TRUE,
                      folds = folds, 
                      early_stopping_rounds = doxgb$early_stopping_rounds,
                      verbose = 0)
  
  if (track >= 3) { print("XGB Simple - set up return object") }
  
  set.seed(seed)                                                                ## 240315 
  
  xgb.simple.fit = xgb.train(params = param.final, data = train.xgb.dat, nrounds = xgb.cv.fit$best_iteration)
  doxgb$nrounds.final = xgb.cv.fit$best_iteration
  xgb.simple.fit$doxgb = doxgb   
  xgb.simple.fit$param.final = param.final
##  xgb.simple.fit$data=train.xgb.dat                                             ##<<<<<<<<<<<<<--------------------
##  xgb.simple.fit$seed=seed                                                      ##<<<<<<<<<<<<<--------------------
  return(xgb.simple.fit)
  if (track >= 3) { print("XGB Simple - return") }
}

################################################################################################################
## TUNED XGBoost fit ===========================================================
#' Get a tuned XGBoost model fit
#' 
#' @description This fits a gradient boosting machine model using the XGBoost
#' platform.  It uses the mlrMBO mlrMBO package to search for a well fitting set of 
#' hyperparameters and will generally provide a better fit than xgb.simple(). 
#' Both this program and xgb.simple() require data to be provided in a 
#' xgb.DMatrix() object.  This object can be constructed with a command like 
#' data.full <- xgb.DMatrix( data=myxs, label=mylabel), where myxs object contains the 
#' predictors (features) in a numerical matrix format with no missing 
#' values, and mylabel is the outcome or dependent variable.  For logistic regression
#' this would typically be a vector of 0's and 1's.  For linear regression this would be 
#' vector of numerical values. For a Cox proportional hazards model this would be 
#' in a format required for XGBoost, which is different than for the survival package 
#' or glmnet package.  For the Cox model a vector is used where observations 
#' associated with an event are assigned the time of event, and observations 
#' associated with censoring are assigned the NEGATIVE of the time of censoring.  In    
#' this way information about time and status are communicated in a single vector
#' instead of two vectors.  The xgb.tuned() function does not handle (start,stop) 
#' time, i.e. interval, data.  To tune the xgboost model we use the mlrMBO package
#' which "suggests" the DiceKriging and rgenoud packages, but doe not install 
#' these.  Still, for xgb.tuned() to run it seems at the time of this writing
#' that one should install the DiceKriging and rgenoud packages.  
#' 
#' @param train.xgb.dat The data to be used for training the XGBoost model
#' @param booster for now just "gbtree" (default) 
#' @param objective one of "survival:cox" (default), "binary:logistic" or "reg:squarederror"
#' @param eval_metric one of "cox-nloglik" (default), "auc" or "rmse",
#' @param minimize whether the eval_metric is to be minimized or maximized
#' @param seed a seed for set.seed() to assure one can get the same results twice.  If NULL 
#' the program will generate a random seed.  Whether specified or NULL, the seed is stored in the output
#' object for future reference.  
#' @param folds an optional list where each element is a vector of indeces for a 
#' test fold.  Default is NULL.  If specified then nfold is ignored a la xgb.cv().
#' @param doxgb  A list specifying how the program is to do the xgb tune and 
#' fit.  The list can have elements $nfold, $nrounds,
#' and $early_stopping_rounds, each numerical values of length 1, $folds, a list as 
#' used by xgb.cv() do identify folds for cross validation, and $eta, $gamma, $max_depth, 
#' $min_child_weight, $colsample_bytree, $lambda, $alpha and $subsample, each a numeric 
#' of length 2 giving the lower and upper values for the respective tuning 
#' parameter.  The meaning of these terms is as in 'xgboost' xgb.train().  If 
#' not provided defaults will be used.  Defaults
#' can be seen from the output object$doxgb element, again a list. In case not NULL, 
#' the seed and folds option values override the $seed and $folds values. 
#' @param track 0 (default) to not track progress, 2 to track progress.   
#' 
#' @return a tuned XGBoost model fit 
#' 
#' @seealso
#   \code{\link{predict_nested_xgb}} , 
#'   \code{\link{xgb.simple}} , \code{\link{rederive_xgb}} , \code{\link{nested.glmnetr}} 
#' 
#' @author Walter K Kremers with contributions from Nicholas B Larson
#'
#' @importFrom smoof makeSingleObjectiveFunction 
#' @importFrom ParamHelpers makeParamSet makeNumericParam makeIntegerParam 
#' @importFrom mlrMBO makeMBOControl setMBOControlTermination mbo 
#' @importFrom xgboost xgb.cv xgb.train xgb.DMatrix getinfo 
#'
#' @export
#' 
xgb.tuned = function(train.xgb.dat, 
                            booster     = "gbtree",
                            objective   = "survival:cox",
                            eval_metric = NULL,
                            minimize = NULL,
                            seed = NULL, 
                            folds = NULL,  
                            doxgb = NULL, 
                            track = 0 ) { 
  
  if (is.null(doxgb)) { doxgb = list() }
  
  if  (!is.null(seed)) { doxgb$seed = seed   
  } else if (is.null( doxgb$seed )) { seed = round(runif(1)*1e9)  ; doxgb$seed = seed }
  if (seed == 0) { seed = round(runif(1)*1e9) ; doxgb$seed = seed }
  
  if ( (is.null(folds)) & (!is.null(doxgb$folds)) ) { folds = doxgb$folds } 
  
  if (!is.null(folds)) { doxgb$folds = folds ; doxgb$nfold = length(doxgb$folds) } 
  
  if ((is.null(doxgb$folds)) & (is.null(doxgb$nfold)))  { doxgb$nfold = 5 } 
  
  if ((is.null(doxgb$folds)) & (!is.null(doxgb$nfold)))  { doxgb$nfold = max(doxgb$nfold, 3) } 
  
  if (is.null( doxgb$nrounds )) { doxgb$nrounds = 1000 } 
  if (is.null( doxgb$early_stopping_rounds )) { doxgb$early_stopping_rounds = min(100,max(10,doxgb$nrounds/5)) } 

  if (is.null(doxgb$eta) ) { doxgb$eta[1] = 0.01 } 
  if ( is.na(doxgb$eta[2]) ) { doxgb$eta[2] = max( 0.3, doxgb$eta[1]) } 
  
  if (is.null(doxgb$gamma)) { doxgb$gamma[1] = -7 } 
  if (is.na(doxgb$gamma[2])) { doxgb$gamma[2] = max( 6, doxgb$gamma[1]) } 
  
  if (is.null(doxgb$max_depth)) { doxgb$max_depth[1] = 2 } 
  if (is.na(doxgb$max_depth[2])) { doxgb$max_depth[2] = max( 20, doxgb$max_depth[1]) } 
  
  if (is.null(doxgb$min_child_weight)) { doxgb$min_child_weight[1] = 1 } 
  if (is.na(doxgb$min_child_weight[2])) { doxgb$min_child_weight[2] = max( 10, doxgb$min_child_weight[1]) } 
  
  if (is.null(doxgb$colsample_bytree)) { doxgb$colsample_bytree[1] = 0.5 } 
  if (is.na(doxgb$colsample_bytree[2])) { doxgb$colsample_bytree[2] = min(max( 1, doxgb$colsample_bytree[1]),1) } 
  
  if (is.null(doxgb$lambda)) { doxgb$lambda[1] = -10 } 
  if (is.na(doxgb$lambda[2])) { doxgb$lambda[2] = max( 10, doxgb$lambda[1]) } 
  
  if (is.null(doxgb$alpha)) { doxgb$alpha[1] = -10 } 
  if (is.na(doxgb$alpha[2])) { doxgb$alpha[2] = max( 10, doxgb$alpha[1]) } 
  
  if (is.null(doxgb$subsample)) { doxgb$subsample[1] = 0.5 } 
  if (is.na(doxgb$subsample[2])) { doxgb$subsample[2] = min(max( 1, doxgb$subsample[1]), 1) } 

#  early_stopping_rounds = doxgb$early_stopping_rounds

  if (objective == "survival:cox") { family = "cox" 
  } else if (objective == "binary:logistic") { family = "binomial" 
  } else { family = "gaussian" } 
  
  if (track >= 3) { print("XGB Tuned - about to get folds")}
  
  stratified = 1 
  if (is.null(folds)) { 
    if  (is.null(seed))   { seed = round(runif(1)*1e9) 
    } else if (seed == 0) { seed = round(runif(1)*1e9) 
    }
    set.seed(seed) 
    label = getinfo(train.xgb.dat,'label')
    if (objective == "survival:cox") {
      y_ = abs( label )
      event = ( label > 0) * 1 
    } else {
      y_ = label  
      event = NULL
    }
    folds_v = get.foldid(y_, event, family, doxgb$nfold, stratified) 
    folds = list()
    for (i1 in c(1:doxgb$nfold)) { folds[[i1]] = c(1:length(y_))[folds_v == i1] } 
  }
  
  if (track >= 3) { print("XGB Tuned - about to get eval_metric")}
  
  if        (objective == "survival:cox"    ) { minimize = TRUE ; if (is.null(eval_metric)) { eval_metric = "cox-nloglik" }
  } else if (objective == "binary:logistic" ) { 
    if (is.null(eval_metric)) { eval_metric = "auc" ; minimize = FALSE }
    if (eval_metric =="mlogloss") { minimize=TRUE }
  } else if (objective == "reg:squarederror") { minimize = TRUE ;  if (is.null(eval_metric)) { eval_metric = "rmse" }
  } else if (objective == "count:poisson"   ) { minimize = TRUE ;  if (is.null(eval_metric)) { eval_metric = "poisson-nloglik" }
  } else if (is.null(minimize)) { minimize = TRUE }
  
  ## DEFINE objective function #####
  
  if (track >= 5) {
  verbose = TRUE 
  show.info = TRUE  
  } else if (track >= 3) {
    verbose = FALSE 
    show.info = TRUE  
  } else {
    verbose = FALSE 
    show.info = FALSE 
  }
    
  if (track >= 3) { print("XGB Tuned - get early_stopping_rounds")}
  
#  early_stopping_rounds = max(10,nrounds/5)
#  early_stopping_rounds = min(100,early_stopping_rounds)
  
  if (track >= 3) { print("XGB Tuned - defind fn function")}

  fn = function(x) {
    cv <- xgb.cv(params = list(
      booster          = booster,
      objective        = objective,
      eval_metric      = eval_metric,
      nthread          = 16,
      eta              = x["eta"],
      max_depth        = x["max_depth"],
      min_child_weight = x["min_child_weight"],
      gamma            = x["gamma"],
      subsample        = x["subsample"],
      colsample_bytree = x["colsample_bytree"],
      lambda           = x["lambda"],
      alpha            = x["alpha"]),
      data = train.xgb.dat,
      nrounds = doxgb$nrounds,
      prediction = FALSE,
      showsd = TRUE,
      folds = folds , 
      early_stopping_rounds = doxgb$early_stopping_rounds,
      verbose = verbose 
    )
     
    if (track >= 3) { 
      print("XGB Tuned - get whch ")
      if (!is.null((cv$evaluation_log))) {
        print( names( cv$evaluation_log ) )
        print( cv$evaluation_log[1,] )
        } 
    }
    
    mydf = cv$evaluation_log
    if        (eval_metric == "cox-nloglik"    ) { whch = as.numeric( which(names(mydf)=="test_cox_nloglik_mean") )  
    } else if (eval_metric == "auc"            ) { whch = as.numeric( which(names(mydf)=="test_auc_mean") ) 
    } else if (eval_metric == "mlogloss"       ) { whch = as.numeric( which(names(mydf)=="test_mlogloss_mean") ) 
    } else if (eval_metric == "rmse"           ) { whch = as.numeric( which(names(mydf)=="test_rmse_mean") )        
    } else if (eval_metric == "poisson-nloglik") { whch = as.numeric( which(names(mydf)=="test_poisson_nloglik_mean") ) 
    }
##  whch = 4 # always ??     
    
    if (track >= 3) { cat(" whch = ", whch,"\n") ; print("XGB Tuned - get value from whch ") }
    
    if (minimize) {
      value = min( as.matrix(mydf)[, whch]) 
      whch.max = which.min( as.matrix(mydf)[, whch])
    } else  {
      value = max( as.matrix(mydf)[, whch]) 
      whch.max = which.max( as.matrix(mydf)[, whch])
    }
    
    if (track >= 2) { cat(paste0(" ", whch.max,",")) }

    if (track >= 4) { cat("value = \n") ; print(value) }
    
#   Though it works in a stand alone function, devtools::check marks (or did mark) these lines as bad code 
#    if        (objective == "survival:cox"    ) { value = cv$evaluation_log[, min(test_cox_nloglik_mean)] 
#    } else if (objective == "binary:logistic" ) { value = cv$evaluation_log[, max(test_auc_mean)] 
#    } else if (objective == "reg:squarederror") { value = cv$evaluation_log[, min(test_rmse_mean)] 
#    } else if (objective == "count:poisson"   ) { value = cv$evaluation_log[, min(test_poisson_nloglik_mean)] 
#    }
#    whch = as.numeric( which(names(mydf)=="test_cox_nloglik_mean") ) 
#    print( min( as.matrix(cv$evaluation_log)[, whch]) )
#    return( min( cv$evaluation_log[, 4]) )
#    print( min( as.matrix(mydf)[, whch]) )
    
    return( value )
  }
  
  if (track >= 3) { print("XGB Tuned - set par.set ") }
  
  par.set = makeParamSet(
    makeNumericParam("eta",              lower = doxgb$eta[1],               upper = doxgb$eta[2] ),
    makeNumericParam("gamma",            lower = doxgb$gamma[1] ,            upper = doxgb$gamma[2] , trafo = function(x) 2^x),
    makeIntegerParam("max_depth",        lower = doxgb$max_depth[1] ,        upper = doxgb$max_depth[2] ),
    makeIntegerParam("min_child_weight", lower = doxgb$min_child_weight[1] , upper = doxgb$min_child_weight[2] ),
    makeNumericParam("colsample_bytree", lower = doxgb$colsample_bytree[1] , upper = doxgb$colsample_bytree[2] ),
    makeNumericParam("lambda",           lower = doxgb$lambda[1],            upper = doxgb$lambda[2], trafo = function(x) 2^x),
    makeNumericParam("alpha",            lower = doxgb$alpha[1],             upper = doxgb$alpha[2], trafo = function(x) 2^x),
    makeNumericParam("subsample",        lower = doxgb$subsample[1],         upper = doxgb$subsample[2] ) 
  )

  if (track >= 3) { print("XGB Tuned - set obg.fun") }
  
  obj.fun  <- smoof::makeSingleObjectiveFunction(
    name = "xgb_cv_bayes",
    fn = fn,
    par.set = par.set,
    minimize = minimize 
  )

  ## END objective function #####
  
  if (track >= 3) { print("XGB Tuned - set MBO control") }
  
  # Set MBO Control parameters
  control = makeMBOControl()
  control = setMBOControlTermination(control, iters = 20)
  
  # show.info = TRUE to give more progress from MBO ...
  
  if (track >= 3) { print("XGB Tuned - set mbo.run") }
  
  set.seed(seed)                                                                ## 240218 
  
  mbo.run = mbo(fun = obj.fun,
                design = NULL,
                control = control,
                show.info = show.info )
  
  # Extract full results
  
#  perf.df1 <- mbo.run$opt.path$env$path %>% arrange(y)
#  mydf = mbo.run$opt.path$env$path
#  perf.df = mydf[order(mydf[,which(nms=='y')]),]

  ## Extract optimal parameters and fit XGBoost model ############################
  
  if (track >= 3) { print("XGB Tuned - set param.final") }
  
  param.final <- list(
    booster          = booster,
    objective        = objective,
    eval_metric      = eval_metric,
    nthread          = 16,
    eta              = mbo.run$x$eta,
    max_depth        = mbo.run$x$max_depth,
    min_child_weight = mbo.run$x$min_child_weight,
    gamma            = 2^(mbo.run$x$gamma),
    subsample        = mbo.run$x$subsample,
    colsample_bytree = mbo.run$x$colsample_bytree,
    lambda           = 2^(mbo.run$x$lambda),
    alpha            = 2^(mbo.run$x$alpha)
  )

  if (track >= 3) { print("XGB Tuned - fit xgb.cv.fit") }
  
  set.seed(seed)                                                                ## 240218 
  
  xgb.cv.fit = xgb.cv(param.final,
                      data = train.xgb.dat,
                      nrounds = doxgb$nrounds,
                      showsd = TRUE,
                      folds = folds, 
                      early_stopping_rounds = doxgb$early_stopping_rounds,
                      verbose = 0)
  
  if (track >= 3) { print("XGB Tuned - set up return object") }
  
  set.seed(seed)                                                                ## 240218 
  
  xgb.tuned.fit = xgb.train(params = param.final, data = train.xgb.dat, nrounds = xgb.cv.fit$best_iteration)
  doxgb$nrounds.final = xgb.cv.fit$best_iteration
  xgb.tuned.fit$doxgb = doxgb   
  xgb.tuned.fit$param.final = param.final
#  names( xgb.tuned.fit ) 
#  xgb.tuned.fit$niter
#  xgb.tuned.fit$params
#  xgb.tuned.fit$param.final
#  xgb.tuned.fit$param.final = param.final   ## use xgb.tuned.fit$params 
#  xgb.tuned.fit$nrounds.final = xgb.cv.fit$best_iteration ## ## use xgb.tuned.fit$niter 
##  class(xgb.tuned.fit) = "xgb.tuned" ## class already assigned as "xgb.Booster"
  return(xgb.tuned.fit)
  if (track >= 3) { print("XGB Tuned - return") }
}

################################################################################
################################################################################
