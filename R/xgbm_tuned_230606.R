##### xgbm_tuned_230322.R #######################################################################################
################################################################################################################
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
#' @param folds_xgb number of folds used in xgb.cv() call 
#' @param minimize whether the eval_metric is to be minimized or maximized
#' @param nrounds max number of iterations 
#' 
#' @author Walter K Kremers with contributions from Nicholas B Larson
#'
#' @importFrom xgboost xgb.cv xgb.train xgb.DMatrix
#' 
#' @return an XGBoost model fit 
#' 
#' @export
#'
#' @examples 
#' \donttest{
#' # Simulate some data for a Cox model 
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
#' Surv.xgb = ifelse( sim.data$event==1, sim.data$yt, -sim.data$yt )
#' data.full <- xgboost::xgb.DMatrix(data = sim.data$xs, label = Surv.xgb)
#' # for this example we use a small number for folds_n and nrounds to shorten run time 
#' xgbfit = xgb.simple( data.full, objective = "survival:cox", folds_xgb=5, nrounds=20)
#' preds = predict(xgbfit, sim.data$xs)
#' summary( preds ) 
#' preds[1:8]
#' }
#' 
xgb.simple = function(train.xgb.dat, 
                      booster     = "gbtree",
                      objective   = "survival:cox",
                      eval_metric = NULL,
                      folds_xgb = 5,
                      minimize = NULL,
                      nrounds=1000 ) {
  
  if        (objective == "survival:cox"    ) { minimize = TRUE ; if (is.null(eval_metric)) { eval_metric = "cox-nloglik" }
  } else if (objective == "binary:logistic" ) { if (is.null(minimize)) { minimize = FALSE } ; if (is.null(eval_metric)) { eval_metric = "auc" }
  } else if (objective == "reg:squarederror") { minimize = TRUE ;  if (is.null(eval_metric)) { eval_metric = "rmse" }
  } else if (objective == "count:poisson"   ) { minimize = TRUE ;  if (is.null(eval_metric)) { eval_metric = "poisson-nloglik" }
  } else { minimize = TRUE }
  
  param.final = list(booster = booster, objective = objective, eval_metric = eval_metric,
                     nthread = 2, eta = 1, max_depth = 2, nrounds = 2 )
  param.final = list(booster = booster, objective = objective, eval_metric = eval_metric,
                     nthread = 16, eta = .2, max_depth = 10, gamma=0.1, min_child_weight=5, 
                     colsample_bytree=0.75, alpha=0.1, subsample=0.6 )
  
  folds_xgb = max(folds_xgb, 5)
  
  early_stopping_rounds = max(10,nrounds/5)
  early_stopping_rounds = min(100,early_stopping_rounds)

  xgb.cv.fit = xgb.cv(param.final,
                      data = train.xgb.dat,
                      nrounds = nrounds,
                      showsd = TRUE,
                      nfold = folds_xgb,
                      early_stopping_rounds = early_stopping_rounds,
                      verbose = 0)
  
  xgb.simple = xgb.train(params = param.final, data = train.xgb.dat, nrounds = xgb.cv.fit$best_iteration)
  return(xgb.simple)
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
#' these.  Still, for xgb.tuned() to run it seems that one should install the 
#' DiceKriging and rgenoud packages.  
#' 
#' @param train.xgb.dat The data to be used for training the XGBoost model
#' @param booster for now just "gbtree" (default) 
#' @param objective one of "survival:cox" (default), "binary:logistic" or "reg:squarederror"
#' @param eval_metric one of "cox-nloglik" (default), "auc" or "rmse",
#' @param folds_xgb number of folds used in xgb.cv() call 
#' @param minimize whether the eval_metiric is to be minimized or maximized
#' @param nrounds max number of iterations 
#' 
#' @author Walter K Kremers with contributions from Nicholas B Larson
#'
#' @importFrom smoof makeSingleObjectiveFunction 
#' @importFrom ParamHelpers makeParamSet makeNumericParam makeIntegerParam 
#' @importFrom mlrMBO makeMBOControl setMBOControlTermination mbo 
#' @importFrom xgboost xgb.cv xgb.train xgb.DMatrix  
## @importFrom DiceKriging
## @importFrom rgenoud
#'
#' @return an XGBoost model fit 
#' @export
#'
#' @examples 
#' \donttest{
#' # Simulate some data for a Cox model 
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
#' Surv.xgb = ifelse( sim.data$event==1, sim.data$yt, -sim.data$yt )
#' data.full <- xgboost::xgb.DMatrix(data = sim.data$xs, label = Surv.xgb)
#' # for this example we use a small number for folds_n and nrounds to shorten 
#' # run time.  This may still take a minute or so.  
#' # xgbfit=xgb.tuned(data.full,objective="survival:cox",folds_xgb=5,nrounds=20)
#' # preds = predict(xgbfit, sim.data$xs)
#' # summary( preds ) 
#' }
#' 
xgb.tuned = function(train.xgb.dat,   
                     booster     = "gbtree", 
                     objective   = "survival:cox",
                     eval_metric = NULL,
                     folds_xgb=5, 
                     minimize = NULL, 
                     nrounds=1000) {
  
#   if (is.null(seed)) { seed = round(runif(1)*1000000000) }
#   set.seed(seed) 
  
  folds_xgb = max(folds_xgb, 3)
  
  if        (objective == "survival:cox"    ) { minimize = TRUE ; if (is.null(eval_metric)) { eval_metric = "cox-nloglik" }
  } else if (objective == "binary:logistic" ) { 
    if (is.null(eval_metric)) { eval_metric = "auc" ; minimize = FALSE }
    if (eval_metric =="mlogloss") { minimize=TRUE }
  } else if (objective == "reg:squarederror") { minimize = TRUE ;  if (is.null(eval_metric)) { eval_metric = "rmse" }
  } else if (objective == "count:poisson"   ) { minimize = TRUE ;  if (is.null(eval_metric)) { eval_metric = "poisson-nloglik" }
  } else if (is.null(minimize)) { minimize = TRUE }
  
  ## DEFINE objective function #####
  
  verbose = FALSE 
  show.info = FALSE 
  
  early_stopping_rounds = max(10,nrounds/5)
  early_stopping_rounds = min(100,early_stopping_rounds)

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
      nrounds = nrounds,
      prediction = FALSE,
      showsd = TRUE,
      nfold = folds_xgb ,
      early_stopping_rounds = early_stopping_rounds,
      verbose = verbose 
    )
#    print( class(cv$evaluation_log) ) 
    
    mydf = cv$evaluation_log
    if        (eval_metric == "cox-nloglik"    ) { whch = as.numeric( which(names(mydf)=="test_cox_nloglik_mean") )  
    } else if (eval_metric == "auc"            ) { whch = as.numeric( which(names(mydf)=="test_auc_mean") ) 
    } else if (eval_metric == "mlogloss"       ) { whch = as.numeric( which(names(mydf)=="test_mlogloss_mean") ) 
    } else if (eval_metric == "rmse"           ) { whch = as.numeric( which(names(mydf)=="test_rmse_mean") )        
    } else if (eval_metric == "poisson-nloglik") { whch = as.numeric( which(names(mydf)=="test_poisson_nloglik_mean") ) 
    }
##  whch = 4 # always ??     
    
    if (minimize) { value = min( as.matrix(mydf)[, whch]) 
    } else  { value = max( as.matrix(mydf)[, whch]) 
    }
    
#   Though it works in a stand alone function, devtools::check marks these lines as bad code 
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
  
  par.set = makeParamSet(
    makeNumericParam("eta", lower = 0.01, upper = 0.3),
    makeNumericParam("gamma", lower = -7, upper = 6, trafo = function(x) 2^x),
    makeIntegerParam("max_depth", lower = 2, upper = 20),
    makeIntegerParam("min_child_weight", lower= 1, upper = 10),
    makeNumericParam("colsample_bytree", lower = 0.5, upper = 1),
    makeNumericParam("lambda", lower = -10, upper = 10, trafo = function(x) 2^x),
    makeNumericParam("alpha", lower = -10, upper = 10, trafo = function(x) 2^x),
    makeNumericParam("subsample", lower = 0.5, upper = 1)
  )
  
  obj.fun  <- smoof::makeSingleObjectiveFunction(
    name = "xgb_cv_bayes",
    fn = fn,
    par.set = par.set,
    minimize = minimize 
  )

  ## END objective function #####
  
  # Set MBO Control parameters
  control = makeMBOControl()
  control = setMBOControlTermination(control, iters = 20)
  
  # show.info = TRUE to give more progress ...
  
  mbo.run = mbo(fun = obj.fun,
                design = NULL,
                control = control,
                show.info = show.info )
  
  # Extract full results
  
#  perf.df1 <- mbo.run$opt.path$env$path %>% arrange(y)
#  mydf = mbo.run$opt.path$env$path
#  perf.df = mydf[order(mydf[,which(nms=='y')]),]

  ## Extract optimal parameters and fit XGBoost model ############################
  
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

  xgb.cv.fit = xgb.cv(param.final,
                      data = train.xgb.dat,
                      nrounds = nrounds,
                      showsd = TRUE,
                      nfold = folds_xgb,
                      early_stopping_rounds = early_stopping_rounds,
                      verbose = 0)
  
  xgb.tuned = xgb.train(params = param.final, data = train.xgb.dat, nrounds = xgb.cv.fit$best_iteration)
  return(xgb.tuned)
}
  
## END TUNED XGBoost fit =======================================================

################################################################################
################################################################################
