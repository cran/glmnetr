##### nested.glmnetr_230426 ###################################################################################################################
## store folds and use same CV folds throughout for methods that allow, store these in a list 
## modify fits 4 adn 5 for ANN, lasso terms plus lasso + 1 more hidden layer level, lasso offset + 2 more layer levels 
## check factor.foldid() program for case where there are is only one level present  
## random select n obs for use in case control like analysis for ann Cox model 
## for using offset also rediervie the lasso model within the ANN routine as this becomes a hyperparamter
###############################################################################################################################################
#' Using nested cross validation, describe and compare fits of various cross validation informed machine learning models.
#' 
#' @description Performs a nested cross validation for cross validation informed 
#' relaxed lasso, (artificial) Neural Network (ANN) with two hidden layers, Gradient 
#' Boosting Machine (GBM), Recursive Partitioning (RPART) and step wise regression.  That is
#' hyper parameters for all these models are informed by cross validation (CV), and a second 
#' layer of CV is used to evaluate the performance of these CV informed model fits.  For
#' step wise regression CV is used to inform either ap value for entry or degress of 
#' freedom (df) for the final model shoice.  For input 
#' we require predictors (features) to be in numeric matrix format with no missing 
#' values.  This is similar to how the glmnet package expects predictors.  For 
#' survival data we allow input of start time as an option, and require stop time, 
#' and an event indicator, 1 for event and 0 for censoring, as separate terms. This
#' may seem unorthodox as it might seem simpler to accept a Surv() object as input.  However
#' the XGBoost routines require a different data format with only a "stop time" 
#' variable, taking a positive value to indicate being associated with an event, and 
#' the negative of the time when associated with a censoring.  Further with XGBoost 
#' we combine both the predictors and outcome (or label) in an xgb.DMatrix() object.  Further
#' teh evaluation of the loss functinwhen fitting a neural network model also requires
#' the individual vectors for time to event and the event indicator and not a Surv()
#' object.  Instead of translating between these formats we take as input the individual 
#' data elements and construct whatever object is needed for the particular model.  
#' 
#' @param xs     predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
#' the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
#' and not a data frame. 
#' @param start  optional start times in case of a Cox model.  - class numeric of length same as number of patients (n)
#' @param y_     dependent variable as a vector: time, or stop time for Cox model, Y_ 0 or 1 for binomal (logistic), numeric for gaussian. 
#' Must be a vector of length same as number of sample size. 
#' @param event   event indicator, 1 for event, 0 for census, Cox model only.
#' Must be a numeric vector of length same as sample size.
#' @param family  model family, "cox", "binomial" or "gaussian" (default) 
#' @param steps_n number of steps done in stepwise regression fitting 
#' @param folds_n generally the number of folds for each level of cross validation.  If 
#' a vector of length 2 the first element indicates the number of folds for cross 
#' validation, and the second entry can take the value 1 to indicate
#' that no nested cross validation is to be performed.  In this case the function 
#' will only derive the models based upon the full data set.  This is useful for the 
#' "lasso informed" Artificial Neural Network (ANN) models which are based upon 
#' both a lasso fit and an ANN fit. See the predict_ann_tab() function regarding
#' getting predicteds for this "lasso informed" models from a nested.glmnetr() output.  
#' @param dolasso fit and do cross validation for lasso model, 0 or 1
#' @param doxgb 1 to fit and evaluate cross validation informed XGBoost 
#' model, and the number of folds used in model derivation will be the same as 
#' folds_n.  0 (default) for no XGBoost model fits.  An integer > 1 to specify 
#' number of folds used in model fit other than folds_n. If length(doxgb) > 1 then the
#' maximum number of rounds (nrounds) when training the xgb model is set to doxgb[2], 
#' where the default is 1000.  
#' @param doann fit and evaluate a cross validation informed Artificial Neural Network 
#' (ANN) model with two hidden levels.  1 for yes, 0 for no (default). To take some 
#' control over the ANN fit specify a list can be passed as the arguement for doann.  The
#' list can have elements $epochs, $epochs2, $myler, $myler2, $eppr, $eppr2, $lenv1, $lenz2, 
#' $actv, $drpot, $wd, wd2, l1, l12, $fold_n, $minloss and $gotoend.  These different 
#' values are then passed to the ann_tab_cv() function, with the meanings described 
#' in the help for that function, with some exception.  When there are two similar 
#' values like $epoch and $epoch2epoch the first applies to the ANN fits without transfer 
#' learning and the second to the models fit with transfer learning from the lasso 
#' model.  Elements of this list left unspecified will take default values.
#' @param dorpart fit and do a nested cross validation for an RPART model.  As rpart() does its
#' own approximation for cross validation there is no new functions for RPART itself. 
#' @param doaic   fit and do cross validation for AIC fit, 0 or 1.   
#' This is provided primarily as a reference.   
#' @param ensemble fit a sort of ensemble model by fitting a lasso model and then using its 
#' predictions as a) an additional term in the feature set or b) an offset 
#' (base_margin in XGBoost) to the XGBoost model.  This is a vector of length 3, with
#' position 1 indicating if the "standard" XGBoost model is to be fit, position 2 indicating 
#' whether the lasso predictor is be included as an additional feature, and position 3 indicating whether 
#' the lasso predictor should be included as an offset. 1 for yes, 0 (default) for no.
#' @param dostep  fit and do cross validation for stepwise regression fit, 0 or 1, 
#' as discussed in James, Witten, Hastie and Tibshirani, 2nd edition.  
#' @param method  method for choosing model in stepwise procedure, "loglik" or "concordance".
#' Other procedures use the "loglik". 
#' @param lambda  lambda vector for the lasso fit
#' @param gamma   gamma vector for the relaxed lasso fit, default is c(0,0.25,0.5,0.75,1)
#' @param relax   fit the relaxed lasso model when fitting a lasso model  
#' @param limit   limit the small values for lambda after the initial fit.  This 
#' will have minimal impact on the cross validation.  Default is 2 for moderate 
#' limitation, 1 for less limitation, 0 for none.   
#' @param fine    use a finer step in determining lambda.  Of little value unless one 
#' repeats the cross valiaiton many times to more finely tune the hyper paramters.  
#' See the 'glmnet' package documentation  
#' @param track   track progress by printing to console elapsed and split times
#' @param seed    a seed for set.seed() to assure one can get the same results twice.  If NULL 
#' the program will generate a random seed.  Whether specified or NULL, the seed is stored in the output
#' object for future reference.  
#' @param foldid  a vector of integers to associate each record to a fold.  Should be integers between 1 and folds_n.
#' These will only be used in the outer folds. 
#' @param ties method for handling ties in Cox model for relaxed model component.  Default 
#' is "efron", optionally "breslow".  For penalized fits "breslow" is 
#' always used as in the 'glmnet' package.
#' @param stratified 1 to generate fold IDs stratified on outcome or event indicators for the binomial or Cox model.
#' @param time    track progress by printing to console elapsed and split times.  Suggested to use
#' track option instead as time option will be discontinued.  
#' @param ... additional arguments that can be passed to glmnet() 
#' 
#' @return - Cross validation informed LASSO, GBM, RPART or STEPWISE model fits, 
#' together with estimates of model performance derived using nested cross validation.
#' 
#' @seealso
#'   \code{\link{glmnetr}} , \code{\link{cv.glmnetr}}  , \code{\link{glmnetr.simdata}} , \code{\link{summary.nested.glmnetr}} , \code{\link{glmnetr.compcv}} , \code{\link{plot.nested.glmnetr}} , \code{\link{predict_ann_tab}} 
#' 
#' @export
#' 
#' @importFrom utils installed.packages  
#' @importFrom stats runif logLik predict cov cor 
#' @importFrom survival Surv coxph coxph.control concordance concordancefit 
#' @importFrom glmnet cv.glmnet 
#' @importFrom Matrix rankMatrix 
#' @importFrom rpart rpart prune 
#' @importFrom xgboost xgb.DMatrix
#' @importFrom torch torch_manual_seed
#'
#' @examples
#' \donttest{
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$y_ 
#' # for this example we use a small number for folds_n to shorten run time 
#' nested.glmnetr.fit = nested.glmnetr( xs, NULL, y_, NULL, family="gaussian", folds_n=3)
#' plot(nested.glmnetr.fit) 
#' plot(nested.glmnetr.fit, coefs=TRUE) 
#' summary(nested.glmnetr.fit) 
#' summary(nested.glmnetr.fit, cvfit=TRUE) 
#' }
#' 
nested.glmnetr = function(xs, start=NULL, y_, event=NULL, family="gaussian", steps_n=0, folds_n=10,  
                          dolasso=1, doxgb=0, doann=0, dorpart=0, dostep=0, doaic=0, ensemble=0, 
                          method="loglik", lambda=NULL, gamma=NULL, relax=TRUE, limit=1, fine=0, 
                          track=0, seed=NULL, foldid=NULL, ties="efron", stratified=1, time=NULL, ... ) {

  if ( is.null(xs) ) { print(" xs cannot be NULL, program will stop ") ; stop ; }
  if (sum(is.na(xs))> 0) { print(" xs cannot have missing values, program will stop ") ; stop ; }
  if ( is.null(y_) ) { print(" y_ cannot be missing or NULL, program will stop ") ; stop ; }
  if (sum(is.na(y_))> 0) { print(" y_ cannot have missing values, program will stop ") ; stop ; }
  if (family=="cox") {
    if        (is.null(event))       { print(" event cannot be NULL for Cox model, program will stop ") ; stop ; 
    } else if (sum(is.na(event))> 0) { print(" event cannot have missing values, program will stop ") ; stop ; 
    } else if (min(y_) < 0)          { print(" survival times cannot be negative for Cox model, program will stop ") ; stop ; 
    }
  }
  
  if (is.null(folds_n)) { folds_n = 10 }
  if (length(folds_n) > 1) { 
    do_ncv = folds_n[2]  
    folds_n = folds_n[1]
  } else { 
    do_ncv = 1 
  }
  folds_n = max(folds_n, 3)

  doxgb_ = 0 
  folds_xgb = folds_n
  nrounds = 1000 
  if (is.null(doxgb)) { 
    doxgb = 0 
    }
  if (doxgb[1] > 0) {  
    doxgb_ = 1 
    if (doxgb[1] > 2) { folds_xgb = doxgb[1] } 
    if (length(doxgb) == 1) {
      nrounds = 1000 
    } else {
      nrounds = doxgb[2]
    }
  }
  
#  if (track >= 3) { cat(paste(" \n   doxgb_=", doxgb_, "  nrounds=", nrounds, "\n\n")) }
  
  doann_ = 0 
  if (is.null(doann)) { doann = 0 }
  if ((!is.list(doann)) & ( doann[1] >= 1 )) { doann_ = 1 ;  doann = as.list(doann) }
  if (is.list(doann)) { doann_ = 1 } 
  
  if (doann_ == 1) {
    temp_ = list(epochs=200, epochs2=200, mylr=0.005, mylr2=0.001, eppr=-2,
                 eppr2=-2, lenz1=16, lenz2=8, actv=1, drpot=0,
                 wd=0, wd2=0, l1=0, l12=0, fold_n=folds_n, 
                 minloss=0, gotoend=0, bestof=1) 
    if (is.null(doann$epochs )) { epochs = temp_$epochs  ; doann$epochs  = epochs  } else { epochs = doann$epochs  }
    if (is.null(doann$epochs2)) { epochs2= temp_$epochs2 ; doann$epochs2 = epochs2 } else { epochs2 = doann$epochs2 }  
    if (is.null(doann$mylr   )) { mylr   = temp_$mylr    ; doann$mylr    = mylr    } else { mylr    = doann$mylr    }
    if (is.null(doann$mylr2  )) { mylr2  = temp_$mylr2   ; doann$mylr2   = mylr2   } else { mylr2   = doann$mylr2   }
    if (is.null(doann$eppr   )) { eppr   = temp_$eppr    ; doann$eppr    = eppr    } else { eppr    = doann$eppr    }
    if (is.null(doann$eppr2  )) { eppr2  = temp_$eppr2   ; doann$eppr2   = eppr2   } else { eppr2   = doann$eppr2   }
    if (is.null(doann$lenz1  )) { lenz1  = temp_$lenz1   ; doann$lenz1   = lenz1   } else { lenz1   = doann$lenz1   }
    if (is.null(doann$lenz2  )) { lenz2  = temp_$lenz2   ; doann$lenz2   = lenz2   } else { lenz2   = doann$lenz2   }
    if (is.null(doann$actv   )) { actv   = temp_$actv    ; doann$actv    = actv    } else { actv    = doann$actv    }
    if (is.null(doann$drpot  )) { drpot  = temp_$drpot   ; doann$drpot   = drpot   } else { drpot   = doann$drpot   }
    if (is.null(doann$wd     )) { wd     = temp_$wd      ; doann$wd      = wd      } else { wd      = doann$wd      }
    if (is.null(doann$wd2    )) { wd2    = temp_$wd2     ; doann$wd2     = wd2     } else { wd2     = doann$wd2     }
    if (is.null(doann$l1     )) { l1     = temp_$l1      ; doann$l1      = l1      } else { l1      = doann$l1      }
    if (is.null(doann$l12    )) { l12    = temp_$l12     ; doann$l1      = l12     } else { l12     = doann$l12     }
    if (is.null(doann$fold_n )) { fold_n = temp_$fold_n  ; doann$fold_n  = fold_n  } else { fold_n  = doann$fold_n  }
    if (is.null(doann$minloss)) { minloss= temp_$minloss ; doann$minloss = minloss } else { minloss = doann$minloss }
    if (is.null(doann$gotoend)) { gotoend= temp_$gotoend ; doann$gotoend = gotoend } else { gotoend = doann$gotoend }
    if (is.null(doann$bestof )) { bestof = temp_$bestof  ; doann$bestof  = bestof  } else { bestof  = doann$bestof  }
    if (fold_n == 0) { fold_n = folds_n ; doann$fold_n = fold_n } 
    if (fold_n == 1) { 
      fold_n = 3
      doann$fold_n = fold_n
      cat(paste(" --- fold_n cannot be 1 and is so set to 3 ---"))
    }
    if (track >= 0) {
      cat(paste0(" epochs=", epochs, " epochs2=", epochs2, " mylr=", mylr, " mylr2=", mylr2, " eppr=", eppr, "\n" )) 
      cat(paste0(" eppr2=", eppr2, " lenz1=", lenz1, " lenz2=", lenz2, " actv=", actv, " drpot=", drpot, "\n" )) 
      cat(paste0(" wd=", wd, " wd2=", wd2, " l1=", l1, " l12=", l12, "fold_n=", fold_n, "\n minloss=", minloss, " gotoend=", gotoend, "\n" )) 
    }
  } 

  if (!("torch" %in% installed.packages()[,1])) {
    doann_ = 0 
    print("\n  ***** torch is NOT installed and so neural networks models will not be fit *****\n")
  }
  
  if (!is.null(time)) {
    track=time
    cat(" Please use track option instead of time option.  time option will be dropped in future versions.")
  }
  if (track  < 0) { eppr = -3 ; eppr2 = -3 ; }
  if (track == 0) { eppr = -2 ; eppr2 = -2 ; }
  if (track >= 1) { cat(paste0("\n", " ##############################################################################################" , "\n")) }
  if (track >= 1) {
    time_start = diff_time() 
    time_last = NULL 
  }
  
  if (ties != "breslow") { ties="efron" }
  
  if (is.null(ensemble)  ) { ensemble = c(1,0,0) } 
  if (length(ensemble)==1) { ensemble = c(rep(ensemble,3)) }
  if (sum(ensemble)   ==0) { ensemble = c(1,0,0) } 
  temp_ = ensemble
  lng = min(length(temp_), 8)
  ensemble = rep(0,8) ; 
  ensemble[1:lng] = temp_ 

  if (((doxgb_==1) | (doann_==1) | (dorpart==1)) & (sum(ensemble[c(2:8)]) >= 1)) { dolasso = 1 } 
  
  tol_ = 10^(-5) 
  
  xs_ncol = dim(xs)[2]                                                          ## number of candidate predictor variables 
  nobs    = dim(xs)[1] 
  one     = c( rep(1, nobs) )
  
#  nvars = dim(xs)[2]
#  if (is.null(dfmax)) { dfmax = nvars + 1 }
#  if (is.null(penalty.factor)) {penalty.factor = rep(1,nvars)}
  
  if (steps_n==0) { steps_n=xs_ncol }
  steps_n = min(steps_n, xs_ncol) 
  
  if ( is.null(family)) {family = NA}
  if ( is.na(family) & (!is.null(start)) ) { family = "cox" }
  if ( is.na(family) & (!is.null(event)) & is.null(y_)) { 
    y_ = event 
    event = NULL 
    family = "binomial" 
  } else if ( is.na(family) & (!is.null(event     )) ) { 
    family = "cox" 
  }
  if ( is.na(family) & ( length(table(y_)) == 2 ) ) { family = "binomial" }
  if ( is.na(family) & ( length(table(y_)) >  2 ) ) { family = "gaussian" }
  if ( family == "poison" ) { family = "poisson" }
  if ( family %in% c("gausian", "gause", "normal") ) { family = "gaussian" }
  if ( family %in% c("COX", "cox", "coxph") ) { family = "cox" }
  
  if (is.null(foldid)) { 
    if (is.null(seed)) { seed = round(runif(1)*1e9) }
    set.seed(seed) 
    if (stratified >= 1) {
      if   (family == "binomial") {  foldid = factor.foldid(y_   , folds_n)  ;  table(y_   , foldid) ;
      } else if (family == "cox") {  foldid = factor.foldid(event, folds_n)  ;  table(event, foldid) ;  
      } else {foldid = sample( rep( 1:folds_n , ceiling(nobs/folds_n) )[ 1:nobs ] , nobs ) ;  
      }
      table(foldid)  ; 
    } else {
      foldid = sample( rep( 1:folds_n , ceiling(nobs/folds_n) )[ 1:nobs ] , nobs ) ; 
    }
  }
  
  if (doann_ == 1) { 
    if (length(seed) == 1) { 
      tseed = round(runif(1)*1e9) 
      seed = c(seed, tseed) 
    } else {
      tseed = seed[2] 
    }
    torch_manual_seed( tseed ) 
  } else { tseed = 0 }

  ## log likelihoods & concordances from LASSO and relaxed lasso by Cross Validation 
  lasso.devian.cv = matrix ( data=rep(0,(7*folds_n)), nrow = folds_n, ncol = 7 ) 
  lasso.cal.devian.cv = matrix ( data=rep(0,(7*folds_n)), nrow = folds_n, ncol = 7 ) 
  lasso.agree.cv  = matrix ( data=rep(0,(7*folds_n)), nrow = folds_n, ncol = 7 )
  lasso.intcal.cv = matrix ( data=rep(0,(7*folds_n)), nrow = folds_n, ncol = 7 )
  lasso.lincal.cv = matrix ( data=rep(0,(7*folds_n)), nrow = folds_n, ncol = 7 )
  ## number of non-zero coefficients in LASSO models 
  lasso.nzero.cv =  matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6)
  lassogammacv = matrix ( data=rep(0,(2*folds_n)), nrow = folds_n, ncol = 2)
  nms = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0", "ridge" )
  nms = c("1se", "min", "1seR", "minR", "1seR.G0", "minR.G0", "ridge" )
  colnames(lasso.devian.cv) = nms 
  colnames(lasso.cal.devian.cv) = nms 
  colnames(lasso.agree.cv ) = nms
  colnames(lasso.intcal.cv) = nms
  colnames(lasso.lincal.cv) = nms
  colnames(lasso.nzero.cv)  = nms[c(1:6)]
  colnames(lassogammacv)    = c("1se", "min")
  
  ## log likelihoods & concordances from XGBoost by Cross Validation 
  xgb.devian.cv = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 ) 
  xgb.lincal.cv = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 )
  xgb.agree.cv  = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 )
  colnames(xgb.devian.cv) = c("Simple", "Feature", "Offset", "Tuned", "Feature", "Offset") 
  colnames(xgb.lincal.cv) = c("Simple", "Feature", "Offset", "Tuned", "Feature", "Offset") 
  colnames(xgb.agree.cv ) = c("Simple", "Feature", "Offset", "Tuned", "Feature", "Offset") 
  
  ## log likelihoods & concordances from RPART by Cross Validation   
  rpnz = function(obj) { 
    frame = obj$frame ;  leaves = frame$var == "<leaf>" ;  used <- unique(frame$var[!leaves]) ; 
    return( length(used) )
  }
  rpart.nzero.cv  = matrix ( data=rep(0,(9*folds_n)), nrow = folds_n, ncol = 9 ) 
  rpart.devian.cv = matrix ( data=rep(0,(9*folds_n)), nrow = folds_n, ncol = 9 ) 
  rpart.lincal.cv = matrix ( data=rep(0,(9*folds_n)), nrow = folds_n, ncol = 9 )
  rpart.agree.cv  = matrix ( data=rep(0,(9*folds_n)), nrow = folds_n, ncol = 9 )
  nms = c("cp=0.00", "cp=0.01", "cp=0.02", "F cp=0.00", "F cp=0.01", "F cp=0.02", 
          "O cp=0.00", "O cp=0.01", "O cp=0.02" )
  colnames(rpart.nzero.cv)  = nms 
  colnames(rpart.devian.cv) = nms 
  colnames(rpart.lincal.cv) = nms 
  colnames(rpart.agree.cv)  = nms 
  
  ## log likelihoods & concordances from Neural Networks by Cross Validation   
  ann.devian.cv   = matrix ( data=rep(1,(6*folds_n)), nrow = folds_n, ncol = 6 ) 
  ann.cal.devian.cv=matrix ( data=rep(1,(6*folds_n)), nrow = folds_n, ncol = 6 ) 
  ann.agree.cv    = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 )
  ann.intcal.cv   = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 )
  ann.lincal.cv   = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 )
  nms = c("Uninformed", "lasso terms", "lasso feat.", "init w's", "update w's", "offset") 
  colnames(ann.devian.cv) = nms 
  colnames(ann.cal.devian.cv) = nms 
  colnames(ann.intcal.cv) = nms 
  colnames(ann.lincal.cv) = nms 
  colnames(ann.agree.cv)  = nms 

  ## log likelihoods & conncordances for STEPWISE by Cross Validation 
  ## and coxph models based upon LASSO terms by Cross Validation
  step.devian.cv = matrix ( data=rep(0,(4*folds_n)), nrow = folds_n, ncol = 4 )
  step.lincal.cv = matrix ( data=rep(0,(3*folds_n)), nrow = folds_n, ncol = 3 )
  step.agree.cv  = matrix ( data=rep(0,(3*folds_n)), nrow = folds_n, ncol = 3 )
  step_df_cv     = matrix ( data=rep(0,(3*folds_n)), nrow = folds_n, ncol = 3 )
  step_p_cv      = matrix ( data=rep(0,(folds_n  )), nrow = folds_n, ncol = 1 )
  colnames(step.devian.cv) = c("df", "p", "AIC", "null")
  colnames(step.lincal.cv) = c("df", "p", "AIC")
  colnames(step.agree.cv ) = c("df", "p", "AIC")
  colnames( step_df_cv  )  = c("df", "p", "AIC")
  colnames( step_p_cv   )  = c("p stepwise")

  #########################################################################################################
  #########################################################################################################
  foldid_ = NULL 
  
  ##### LASSO fit #########################################################################################
  if (dolasso == 1) {
    if (track >= 1) { cat(paste0("\n", " ########## Initial (CV) lasso fit of all data ########################################" , "\n")) }
    if (relax) { cv_glmnet_fit_f = cv.glmnetr( xs, start, y_, event, family=family, lambda=lambda, 
                                               gamma=gamma, folds_n=folds_n, limit=limit, fine=fine, track=track, ties=ties, ... ) 
    } else {
      if (family=="cox") {
        if ( is.null(start) ) { cv_glmnet_fit_f = cv.glmnet( xs, Surv(y_, event)       , family="cox", lambda=lambda, relax=FALSE, ... ) 
        } else                { cv_glmnet_fit_f = cv.glmnet( xs, Surv(start, y_, event), family="cox", lambda=lambda, relax=FALSE, ... ) 
        } 
      } else {
       cv_glmnet_fit_f = cv.glmnet( xs, y_, family=family, lambda=lambda, relax=FALSE, ... ) 
      }
    }
    
    if ((family == "cox") & (is.null(start))) {
      cv_ridge_fit_f = cv.glmnet( xs, Surv(y_, event), family=family, alpha=0, ... ) 
    } else if ((family == "cox") & (!is.null(start))) {
      cv_ridge_fit_f = cv.glmnet( xs, Surv(start, y_, event), family=family, alpha=0, ... ) 
    } else {
      cv_ridge_fit_f = cv.glmnet( xs, y_, family=family, alpha=0, ... ) 
    }
    
    #------------------------------------------
    
    if (family=="cox") {    
      if ( is.null(start) ) { SURV = Surv(y_, event) 
      } else {                SURV = Surv(start, y_, event) }
      lasso_perf_cox = function(SURV, pred) {
        fit0 = coxph( SURV ~ pred, init=c(1)) 
        retvec = c(dev1 = -2*fit0$loglik[1] / fit0$nevent ,
                   dev2 = -2*fit0$loglik[2] / fit0$nevent , 
                   agree = fit0$concordance[6] ,
                   intcal = 0 , 
                   lincal = fit0$coefficients ) 
        return( retvec )
      }
    }
    
    if (family == "binomial") {
      lasso_perf_bin = function(yy, pred) {
        fit0 = glm( yy ~ pred , family=family)
        p_ = 1/(1+exp(-pred)) ;
        retvec = c( dev1 = -2*( t(log(p_))%*%yy + t(log(1-p_))%*%(1-yy) ) / length(yy) , 
          dev2 = fit0$deviance / length(yy) ,
          agree  = concordance(fit0)[[1]] ,          
          intcal = fit0$coefficients[1] , 
          lincal = fit0$coefficients[2] ) 
        return( retvec )
      }
    }
    
    if (family=="gaussian") {
      lasso_perf_gau = function(yy, pred) {
        fit0 = glm( yy ~ pred , family=family)
        returnvec = c(  dev1 = sum((yy - pred)^2) / length(yy) ,
                 dev2 = fit0$deviance / length(yy) ,
                 agree = ifelse( var(pred)>0 , cor(x=yy, y=pred)^2 , 0 ) ,
                 intcal = fit0$coefficients[1] , 
                 lincal = fit0$coefficients[2] )  
      }
    }
    
    lasso_perf = function(yy, pred, family) {
      if        (family == "cox")      { lasso_perf_cox(yy, pred) 
      } else if (family == "binomial") { lasso_perf_bin(yy, pred) 
      } else if (family == "gaussian") { lasso_perf_gau(yy, pred) } 
    }
    
    pred1se   = predict(cv_glmnet_fit_f, xs, lam=cv_glmnet_fit_f$lambda.1se, gam=1)
    predmin   = predict(cv_glmnet_fit_f, xs, lam=cv_glmnet_fit_f$lambda.min, gam=1)
    pred1seR  = predict(cv_glmnet_fit_f, xs, lam="lambda.1se" , gam="gamma.1se" )
    predminR  = predict(cv_glmnet_fit_f, xs, lam="lambda.min" , gam="gamma.min" )
    pred1seR0 = predict(cv_glmnet_fit_f, xs, lam=cv_glmnet_fit_f$lambda.1se, gam=0)
    predminR0 = predict(cv_glmnet_fit_f, xs, lam=cv_glmnet_fit_f$lambda.min, gam=0)
    predridge = predict(cv_ridge_fit_f , xs)
      
    #      print( cor(cbind(pred1se, predmin, pred1seR, predminR, pred1seR0, predminR0)) ) 
    #      print( cor(cbind(predmin, predminR, predminR0)) ) 
      
    if (family == "cox") {    
      if ( is.null(start) ) { yy = Surv(y_, event) 
      } else {                yy = Surv(start, y_, event) }
    } else { yy = y_ }
      
    perfm1 = lasso_perf( yy , pred1se   , family )
    perfm2 = lasso_perf( yy , predmin   , family )
    perfm3 = lasso_perf( yy , pred1seR  , family )
    perfm4 = lasso_perf( yy , predminR  , family )
    perfm5 = lasso_perf( yy , pred1seR0 , family )
    perfm6 = lasso_perf( yy , predminR0 , family )
    perfm7 = lasso_perf( yy , predridge , family )
    lasso.devian.naive  = c( perfm1[1] , perfm2[1] , perfm3[1] , perfm4[1] , perfm5[1] , perfm6[1] , perfm7[1] )
    lasso.cal.devian.naive=c(perfm1[2] , perfm2[2] , perfm3[2] , perfm4[2] , perfm5[2] , perfm6[2] , perfm7[2] )    
    lasso.agree.naive   = c( perfm1[3] , perfm2[3] , perfm3[3] , perfm4[3] , perfm5[3] , perfm6[3] , perfm7[3] )
    lasso.intcal.naive  = c( perfm1[4] , perfm2[4] , perfm3[4] , perfm4[4] , perfm5[4] , perfm6[4] , perfm7[4] )
    lasso.lincal.naive  = c( perfm1[5] , perfm2[5] , perfm3[5] , perfm4[5] , perfm5[5] , perfm6[5] , perfm7[5] )
    nms = c("naive_agree_1se" , "naive_agree_min"  , "naive_agree_1seR" , 
            "naive_agree_minR", "naive_agree_1seR0", "naive_agree_minR0", "naive_agree_ridge" ) 
    nms = c("1se", "min", "1seR", "minR", "1seR.G0", "minR.G0", "ridge" )
    names(lasso.devian.naive) = nms 
    names(lasso.cal.devian.naive) = nms 
    names(lasso.agree.naive)  = nms 
    names(lasso.intcal.naive) = nms 
    names(lasso.lincal.naive) = nms 

    if (track >= 1) { cat(paste0(" length(lambda) = " , length(cv_glmnet_fit_f$lambda), "\n" )) } 
    
    if (track >= 1) { time_last = diff_time(time_start, time_last)  }
  }
  
  ##### XGBoost fit ##########################################################################################
  if ( (doxgb_==1) & (family=="cox") & !is.null(start) ) { 
    doxgb_ = 0  
    cat("\n xgboost() does not fit Cox model with (start, stop) time data.  doxgb set to 0. \n\n" ) 
  } 
  
  if (doxgb_==1) { 
    if (track >= 1) { cat(paste0("\n", " ########## Initial XGBoost fit on all data ###########################################" , "\n")) }
    
    if (family=="cox") { Surv.xgb = ifelse( event == 1, y_, -y_) }
    
    if (ensemble[1] >=1) {
      if (family=="cox") { full.xgb.dat <- xgb.DMatrix(data = xs, label = Surv.xgb)
      } else {             full.xgb.dat <- xgb.DMatrix(data = xs, label = y_)      }
    }
    
    if ((ensemble[2] >= 1) | (ensemble[3] >= 1)) {
      ofst = predict(cv_glmnet_fit_f, xs, lam="lambda.min", gam="gamma.min" )  ###### min RELAXED lasso as an offset #########
    }
    
    if (ensemble[2] >= 1) {
      if (family=="cox") { full.xgb.datF <- xgb.DMatrix(data = cbind(xs,ofst), label = Surv.xgb)
      } else {             full.xgb.datF <- xgb.DMatrix(data = cbind(xs,ofst), label = y_)      }
    }
    
    if (ensemble[3] >= 1) {
      if (family=="cox") { full.xgb.datO <- xgb.DMatrix(data = xs, label = Surv.xgb, base_margin=ofst)
      } else {             full.xgb.datO <- xgb.DMatrix(data = xs, label = y_, base_margin=ofst)      }
    }
    
    xgb.devian.naive = c(0,0,0, 0,0,0)
    xgb.lincal.naive = c(0,0,0, 0,0,0)
    xgb.agree.naive  = c(0,0,0, 0,0,0) 
    names(xgb.devian.naive) = c("Simple", "Feature", "Offset", "Tuned", "Feature", "Offset") 
    names(xgb.lincal.naive) = c("Simple", "Feature", "Offset", "Tuned", "Feature", "Offset") 
    names(xgb.agree.naive)  = c("Simple", "Feature", "Offset", "Tuned", "Feature", "Offset") 

    xgb_cox_perf = function(hrhat, y=y_, evnt=event, tol=tol_) {
      hrhat[(hrhat < tol)] = tol
      hrhat[(hrhat >  (1/tol))] =  (1/tol)
      XB = log( hrhat )
      fit1 = coxph( Surv(y, evnt) ~ XB, init=c(1), control=coxph.control(iter.max=0)) 
      fit2 = coxph( Surv(y, evnt) ~ XB) 
      returnvec = c(-2*fit1$loglik[1] / fit1$nevent, 
                    fit2$coefficients[1], fit2$concordance[6], -2*fit2$loglik[2] / fit2$nevent)
      return(returnvec) 
    }

    xgb_bin_perf = function(phat, y=y_, tol=tol_) {
      phat[(phat < tol)] = tol
      phat[(phat > 1 - tol)] = 1 - tol
      XB = log( phat / (1-phat) )
      p_ = phat
      fit1   = glm( y ~ XB, family=family ) 
      returnvec = c(-2*( t(log(p_))%*%y + t(log(1-p_))%*%(1-y) ) / length(y) , 
                    fit1$coefficients[2] , 
                    concordance(y ~ phat)$concordance  )
      return(returnvec) 
    }

    xgb_gau_perf = function(xbhat, y__=y_) {
      if (var(xbhat) > 0) { 
        beta = cov(y__ ,xbhat )/var( xbhat ) 
        rsquare = cor(x=y__ , y=xbhat )^2 
      } else {
        beta = 0
        rsquare = 0 
      }
      returnvec = c( sum((y__  - xbhat)^2) / length(y__ ) , beta, rsquare )  
      return(returnvec)
    }
    
    xgb_perform = function(xgb_model, xs_, y=y_, evnt=event, family, tol=tol_) {
      lxb = predict(xgb_model, xs_)
      if (family == "cox") {
        xgb_cox_perf(lxb, y, evnt, tol) 
      } else if (family == "binomial") {
        xgb_bin_perf(lxb, y, tol) 
      } else if (family == "gaussian") {
        xgb_gau_perf(lxb, y) 
      }
    }

    if (family %in% c("cox","binomial", "gaussian")) { 
      
      if        (family == "cox"     ) { objective = "survival:cox" 
      } else if (family == "binomial") { objective = "binary:logistic" 
      } else                           { objective = "reg:squarederror" }
      
      ## SIMPLE XGB model full data ################################################
      if (ensemble[1] == 1) {
        xgb.simple.fit = xgb.simple(full.xgb.dat, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
        xgbperf = xgb_perform(xgb_model=xgb.simple.fit, xs_=full.xgb.dat, y=y_, evnt=event, family=family )
        xgb.devian.naive[1] = xgbperf[1] ;  xgb.lincal.naive[1] = xgbperf[2] ;  xgb.agree.naive [1] = xgbperf[3] 
      }
      if (ensemble[2] == 1) {
        xgb.simple.fitF = xgb.simple(full.xgb.datF, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
        xgbperf = xgb_perform(xgb.simple.fitF, full.xgb.datF, y=y_, evnt=event, family=family )
        xgb.devian.naive[2] = xgbperf[1] ;  xgb.lincal.naive[2] = xgbperf[1] ;  xgb.agree.naive [2] = xgbperf[2]
      }
      if (ensemble[3] == 1) {
        xgb.simple.fitO = xgb.simple(full.xgb.datO, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
        xgbperf = xgb_perform(xgb.simple.fitO, full.xgb.datO, y=y_, evnt=event, family=family )
        xgb.devian.naive[3] = xgbperf[1] ;  xgb.lincal.naive[3] = xgbperf[2] ;  xgb.agree.naive [3] = xgbperf[3]        
      }
      ## TUNED XGB fit full data ###################################################
      if (ensemble[1] == 1) {
        xgb.tuned.fit = xgb.tuned(full.xgb.dat, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
        xgbperf = xgb_perform(xgb.tuned.fit, full.xgb.dat, y=y_, evnt=event, family=family )
        xgb.devian.naive[4] = xgbperf[1] ;  xgb.lincal.naive[4] = xgbperf[2] ;  xgb.agree.naive [4] = xgbperf[3] 
      }
      if (ensemble[2] == 1) {
        xgb.tuned.fitF = xgb.tuned(full.xgb.datF, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
        xgbperf = xgb_perform(xgb.tuned.fitF, full.xgb.datF, y=y_, evnt=event, family=family )
        xgb.devian.naive[5] = xgbperf[1] ;  xgb.lincal.naive[5] = xgbperf[2] ;  xgb.agree.naive [5] = xgbperf[3] 
      }
      if (ensemble[3] == 1) {
        xgb.tuned.fitO = xgb.tuned(full.xgb.datO, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
        xgbperf = xgb_perform(xgb.tuned.fitO, full.xgb.datO, y=y_, evnt=event, family=family )
        xgb.devian.naive[6] = xgbperf[1] ;  xgb.lincal.naive[6] = xgbperf[2] ;  xgb.agree.naive [6] = xgbperf[3] 
      }
      ##########################################################################
    }
    if (track >= 1) { time_last = diff_time(time_start, time_last)  }
  }
  
  ##### Neural Network fit ######################################################################################
  if (doann_ == 1) { 
    if (track >= 1) { cat(paste0("\n", " ########## Initial Neural Network fit on all data #############################################" , "\n")) }
    ## 1 - uninformed NN                                 ## 1
    ## 2 - limit to lasso variables                      ## 4
    ## 3 - add lasso X*Beta as feature                   ## 2
    ## 4 - include lasso structure in initial weights    ## 5
    ## 5 - update weights each epoch                     ## 3

    ann_cox_perf = function(object, newdat=xs_z, newy, start, event, tol=tol_) {
      xbhat  = as.numeric( object$model(newdat) )  
#      hrhat[(hrhat<tol)] = tol
#      hrhat[(hrhat>(1/tol))] = 1/tol 
#      xbhat = log(hrhat)  
      if (is.null(start)) {
        fit1 = coxph( Surv(newy, event) ~ xbhat, init=c(1)) 
      } else {
        fit1 = coxph( Surv(start, newy, event) ~ xbhat, init=c(1)) 
      }
      return( c(dev1 = -2*fit1$loglik[1]/fit1$nevent, 
                dev2 = -2*fit1$loglik[2]/fit1$nevent, 
                agree= fit1$concordance[6] ,
                intcal = 0 , 
                lincal= fit1$coefficients[1] ) )
    }   

    ann_bin_perf = function(object=ann_fit_1_f, newdat=xs_z, newy=y_, start, event, tol=tol_) {
      phat_nn  = as.numeric( object$model(newdat) )  
      summary(phat_nn)
      phat_nnt = phat_nn 
      phat_nnt[(phat_nn < tol)] = tol
      phat_nnt[(phat_nn > (1-tol))] = 1 - tol
      xbhat  = log(phat_nnt/(1-phat_nnt))
      fit1 = glm( newy ~ xbhat , start=c(0,1), family=family) 
      p_ = 1/(1+exp(-xbhat)) ; 
      returnvec = c( dev1 = -2*( t(log(p_))%*%newy + t(log(1-p_))%*%(1-newy) ) / length(newy) ,  
                     dev2 = fit1$deviance / length(yy) ,
                     agree = concordance(newy ~ xbhat)[[1]] , 
                     intcal = fit1$coefficients[1] ,
                     lincal = fit1$coefficients[2] ) 
      return(returnvec) 
    }
    
    ann_gau_perf = function(object=ann_fit_1_f, newdat=xs_z, newy=y_, tol=tol_) {
      xbhat  = as.numeric( object$model(newdat) )  
      fit1 = glm( newy ~ xbhat , start=c(0,1), family=family) 
      returnvec = c( dev1 = mean((xbhat - newy)^2),
                     dev2 = fit1$deviance / length(yy) ,
                     agree = ifelse( var(xbhat)>0 , cor(xbhat, newy)^2 , 0 ), 
                     lincal = fit1$coefficients[1] , 
                     lincal = fit1$coefficients[2]  )   
      return(returnvec) 
    }
    
    ann_perform = function(object=ann_fit_1_f, newdat=xs_z, newy=y_, family="binomial", start, event, tol=tol_) {
      if (family=="binomial") { 
        ann_bin_perf(object, newdat, newy, tol) 
      } else if (family=="gaussian") { 
        ann_gau_perf(object, newdat, newy, tol) 
      } else if (family=="cox") { 
        ann_cox_perf(object, newdat, newy, start, event, tol) 
      }
    }
    
    ann_fit_1_f = NULL ;  ann_fit_2_f = NULL ;  ann_fit_3_f = NULL ;  ann_fit_4_f = NULL ; ann_fit_5_f = NULL ; ann_fit_6_f = NULL ; 
    
#    print(ensemble)
    xs_means = colMeans(xs)
    xs_c     = sweep(xs, 2, xs_means, FUN = "-") 
    xs_sds   = sqrt( colSums(xs_c^2) / (dim(xs)[1]-1) ) 
    xs_sds_  = xs_sds 
    xs_sds_[(xs_sds <= 1e-8)] = 1e-8 
    xs_z     = sweep(xs_c, 2, xs_sds_, FUN = "/") 
    xs_z[,(xs_sds <= 1e-8)] = 0 
#    table(round(diag(cov(xs_z)),digits=8))
    xs_z0 = xs_z[,(xs_sds_ > 1e-8)]                                              ## xs_sds or xs_sds_ ?? 
    ann_zb = list(xs_means=xs_means, xs_sds_=xs_sds_, tol=1e-8)
      
    if (sum(ensemble[2:6]) > 0) {
      lassopred = predict(cv_glmnet_fit_f, xs, lam="lambda.min" , gam="gamma.min" )
      lassobeta = predict(cv_glmnet_fit_f,     lam="lambda.min" , gam="gamma.min" )
      if (family=="cox") {
        fit0 = coxph(Surv(y_, event) ~ lassopred) 
        calbeta = fit0$coefficients
        lassopred = calbeta * lassopred 
      } else if (family %in% c("binomial", "gaussian")) {
        fit0 = glm(y_ ~ lassopred, family=family) 
        calbeta = fit0$coefficients
        lassopred = calbeta[1] + calbeta[2] * lassopred 
      }
      dim(xs_z)
      if (family %in% c("binomial", "gaussian")) {
        xs_z1 = xs_z[,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])]]          ## pick up non zero features, remove intercept from this list 
      } else if (family == "cox") {
        xs_z1 = xs_z[,(lassobeta[[1]] != 0)]
      }
      dim(xs_z1)
      xs_z2 = cbind(lasso=lassopred,xs_z1)                                      ## add lasso prediction as first column 
      dim(xs_z2)
#      table(diag(cov(xs_z0)))
#      table(diag(cov(xs_z1)))
#      table(diag(cov(xs_z2)))
      ann_zb$lassobeta=lassobeta ; ann_zb$calbeta=calbeta ; 
    }
    
    getwd = 0 ; getwd2 = 0 ; 
    if (wd  < 0) { getwd  = 1 ; wd  = abs(wd ) * cv_ridge_fit_f$lambda.min } 
    if (wd2 < 0) { getwd2 = 1 ; wd2 = abs(wd2) * cv_ridge_fit_f$lambda.min }
    
    getl1 = 0 ; getl12 = 0 ; 
    if (l1  < 0) { getl1  = 1 ; l1  = abs(l1 ) * cv_glmnet_fit_f$relaxed$lambda.min.g0 } 
    if (l12 < 0) { getl12 = 1 ; l12 = abs(l12) * cv_glmnet_fit_f$relaxed$lambda.min.g0 }
      
    if (family %in% c("cox", "binomial", "gaussian")) {
#     cat(paste(fold_n, family, epochs, eppr, wd, lenz1, lenz2, mylr) )
      if (ensemble[1] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN uninformed **\n")) }
        ## xs_z0 - standardized data 
        ann_fit_1_f = ann_tab_cv(myxs=xs_z0, mystart=start, myy=y_, myevent=event, fold_n=fold_n, family=family, 
                                 epochs=epochs, eppr=eppr, lenz1=lenz1, lenz2=lenz2, mylr=mylr, actv=actv, drpot=drpot, wd=wd, l1=l1, minloss=minloss, gotoend=gotoend) 
        ann_perf1 = ann_perform(ann_fit_1_f, xs_z0, y_, family, start, event)
      } else { ann_perf1 = c(10,10,0,0,0) }
      if (ensemble[4] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN informed on relaxed lasso terms **\n")) }
        ann_fit_2_f = ann_tab_cv(myxs=xs_z1, mystart=start, myy=y_, myevent=event, fold_n=fold_n, family=family, 
                                 epochs=epochs, eppr=eppr, lenz1=lenz1, lenz2=lenz2, mylr=mylr, actv=actv, drpot=drpot, wd=wd, l1=l1, minloss=minloss, gotoend=gotoend) 
        perfm2 = ann_perform(ann_fit_2_f, xs_z1, y_, family, start, event)
      } else { perfm2 = c(10,10,0,0,0) }
      if (ensemble[2] == 1) {
        if (eppr >= -2)  { cat(paste0("  ** fitting ANN with feature X*Beta and relaxed lasso terms **\n")) } 
        lenz1_ = lenz1 + 1 ; lenz2_ = lenz2 + 1 ; 
        ann_fit_3_f = ann_tab_cv(myxs=xs_z2, mystart=start, myy=y_, myevent=event , fold_n=fold_n, family=family, 
                                 epochs=epochs, eppr=eppr, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr, actv=actv, drpot=drpot, wd=wd, l1=l1, minloss=minloss, gotoend=gotoend) 
        ann_perf = ann_perform(ann_fit_3_f, xs_z2, y_, family, start, event)
      } else { perfm3 = c(10,10,0,0,0) }
      if (ensemble[3] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN including X*Beta from relaxed lasso in initial weight structure  **\n")) } 
        lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
        ann_fit_4_f = ann_tab_cv(myxs=xs_z2, mystart=start, myy=y_, myevent=event, fold_n=fold_n, family=family, 
                                 epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, drpot=drpot, wd=wd2, l1=l12,lasso=1, adjustw=0, minloss=minloss, gotoend=gotoend) 
        perfm4 = ann_perform(ann_fit_4_f, xs_z2, y_, family, start, event)
      } else { perfm4 = c(10,10,0,0,0) }
      if (ensemble[5] == 1) {
        if (eppr >= -2) { cat(paste0("  ** fitting ANN including X*Beta from relaxed lasso in weight structure each epoch **\n")) } 
        lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
        ann_fit_5_f = ann_tab_cv(myxs=xs_z2, mystart=start, myy=y_, myevent=event, fold_n=fold_n, family=family, 
                                 epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, drpot=drpot, wd=wd2, l1=l12, lasso=1, adjustw=1, minloss=minloss, gotoend=gotoend) 
        perfm5 = ann_perform(ann_fit_5_f, xs_z2, y_, family, start, event)
      } else { perfm5 = c(10,10,0,0,0) }
      if (ensemble[6] == 1.41456) {
        if (eppr >= -2) { cat(paste0("  ** fitting ANN including X*Beta as offset - DOESN'T WORK properly **\n")) } 
        ann_fit_6_f = ann_tab_cv(myxs=xs_z1, mystart=start, myy=y_, myevent=event, fold_n=fold_n, family=family, 
                                 epochs=epochs2, eppr=eppr2, lenz1=lenz1, lenz2=lenz2, mylr=mylr2, actv=actv, drpot=drpot, wd=wd, l1=l1, minloss=minloss, gotoend=gotoend) 
        perfm6 = ann_perform(ann_fit_6_f, xs_z1, y_, family, start, event)   ## offset=lassopred, 
      } else { perfm6 = c(10,10,0,0,0)  }
      ann.devian.naive  = c( perfm1[1] , perfm2[1] , perfm3[1] , perfm4[1] , perfm5[1] , perfm6[1] )
      ann.cal.devian.naive=c(perfm1[2] , perfm2[2] , perfm3[2] , perfm4[2] , perfm5[2] , perfm6[2] )    
      ann.agree.naive   = c( perfm1[3] , perfm2[3] , perfm3[3] , perfm4[3] , perfm5[3] , perfm6[3] )
      ann.lincal.naive  = c( perfm1[4] , perfm2[4] , perfm3[4] , perfm4[4] , perfm5[4] , perfm6[4] )
      ann.intcal.naive  = c( perfm1[5] , perfm2[5] , perfm3[5] , perfm4[5] , perfm5[5] , perfm6[5] )
      nms = c("Uninformed", "lasso terms", "lasso feat.", "init w's", "update w's", "offset") 
      names(ann.devian.naive)  = nms 
      names(ann.cal.devian.naive) = nms
      names(ann.intcal.naive)  = nms
      names(ann.lincal.naive)  = nms
      names(ann.agree.naive)   = nms
    } 
    if (track >= 1) { time_last = diff_time(time_start, time_last)  }
  }
  
  ##### RPART fit ##########################################################################################
  if (dorpart==1) { 
    if (track >= 1) { cat(paste0("\n", " ########## Initial RPART fit on all data #############################################" , "\n")) }

    rpart.train.00 = NULL ; rpart.trainO.00 = NULL ; rpart.trainF.00 = NULL ;
    
    rpart.devian.naive = c(0,0,0, 0,0,0, 0,0,0)
    rpart.nzero.naive  = c(0,0,0, 0,0,0, 0,0,0)
    rpart.lincal.naive = c(0,0,0, 0,0,0, 0,0,0)
    rpart.agree.naive  = c(0,0,0, 0,0,0, 0,0,0) 
    
    nms = c("cp=0.00", "cp=0.01", "cp=0.02", "F cp=0.00", "F cp=0.01", "F cp=0.02", "O cp=0.00", "O cp=0.01", "O cp=0.02" ) 
    names(rpart.devian.naive)  = nms 
    names(rpart.nzero.naive )  = nms
    names(rpart.lincal.naive)  = nms
    names(rpart.agree.naive )   = nms

    #<--------------------------------------------------------------------------
    trainxs = xs 
    testxs = xs 
    trainstart = start 
    teststart = start 
    trainy_ = y_ 
    testy_  = y_ 
    trainevent = event 
    testevent  = event 
    #<--------------------------------------------------------------------------
    
    df.trainxs = as.data.frame(trainxs) 
    df.testxs  = as.data.frame(testxs) 
    
    if ((ensemble[2] == 1) | (ensemble[3] == 1)) { 
#     trainofst = predict(cv_glmnet_fit  , trainxs, lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
#     testofst  = predict(cv_glmnet_fit  , testxs , lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
      trainofst = predict(cv_glmnet_fit_f, trainxs, lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
      length(trainofst) ; 
      testofst  = predict(cv_glmnet_fit_f, testxs , lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
      length(testofst) ; 
      trainxsO  = cbind(trainxs,ofst=trainofst) 
      dim(xs) ; dim(trainxs) ; dim(trainxsO) ; length(trainofst) ; 
      testxsO   = cbind(testxs ,ofst=testofst) 
      df.trainxsO = as.data.frame(trainxsO)
      df.testxsO  = as.data.frame(testxsO )
#      dim(df.testxsO) ; dim(df.trainxsO)
    }
    
    rpart_cox_perf = function(rpart_model, dframe, y__, tol=1e-5) {
      hrhat  = predict(rpart_model, dframe, type="vector")  
      hrhat[(hrhat<tol)] = tol
      hrhat[(hrhat>(1/tol))] = 1/tol 
      xbhat = log(hrhat) 
      fit1 = coxph( y__ ~ xbhat, init=c(1), control=coxph.control(iter.max=0)) 
      fit2 = coxph( y__ ~ xbhat)   
      return( c(devian= -2*fit1$loglik[1]/fit1$nevent, lincal= fit2$coefficients[1], 
                agree= fit2$concordance[6]) )
    }
    
    rpart_bin_perf = function(rpart_model, dframe, y__, tol=1e-5) {
      phat_  = predict(rpart_model, dframe, type="prob")[,2]  
      phat_t = phat_ 
      phat_t[(phat_ < tol)] = tol
      phat_t[(phat_ > (1-tol))] = 1 - tol
      xbhat  = log(phat_t/(1-phat_t))
      fit1 = glm( y__ ~ xbhat , start=c(0,1), family=family) 
      p_ = 1/(1+exp(-xbhat)) ; 
      returnvec = c( -2*( t(log(p_))%*%y__ + t(log(1-p_))%*%(1-y__) ) / length(y__) ,  
                     fit1$coefficients[2] , 
                     concordance(y__ ~ xbhat  )[[1]] )    ## fit1$loglik[1] ## ??
      return(returnvec) 
    }
    
    rpart_gau_perf = function(rpart_model, dframe, y__) {
      xbhat = predict(rpart_model, dframe, type="vector") 
      if (var( xbhat ) > 0) { 
        beta = cov(y__ ,xbhat )/var( xbhat ) 
        rsquare = cor(x=y__ , y=xbhat )^2 
      } else {
        beta = 0
        rsquare = 0 
      }
      returnvec = c( sum((y__  - xbhat)^2) / length(y__ ) , beta, rsquare )  
      return(returnvec)
    }
    
    rpart_perform = function(rpart_model, dframe, y__, family, tol=1e-5) {
      if (family == "cox") {
        returnvec = rpart_cox_perf(rpart_model, dframe, y__, tol=tol_)
      } else if (family == "binomial") {
        returnvec = rpart_bin_perf(rpart_model, dframe, y__, tol=tol_)
      }else if (family == "gaussian") {
        returnvec = rpart_gau_perf(rpart_model, dframe, y__ ) 
      }
      returnvec = c(returnvec, rpnz(rpart_model))
    }
    
    if (family=="cox") {
      if (is.null(start)) { 
        trainy__ = Surv(trainy_, trainevent)
        testy__  = Surv(testy_ , testevent)
        if ((ensemble[2] ==1) | (ensemble[3] ==1)) { df.trainO   = as.data.frame(cbind(trainy_, trainevent, trainxsO) ) }
      } else { 
        trainy__ = Surv(trainstart, trainy_, trainevent) 
        testy__  = Surv(teststart , testy_ , testevent ) 
        if ((ensemble[2] ==1) | (ensemble[3] ==1)) { df.trainO   = as.data.frame(cbind(trainstart, trainy_, trainevent, trainxsO) ) }
      }
    } else {
      trainy__ = trainy_
      testy__  = testy_ 
    }
      
    if (family == "cox") {
      rpmethod = "exp" 
    } else if (family == "binomial") {
      rpmethod = "class" 
    } else if (family == "gaussian") {
      rpmethod = "anova"
    }
      
    if (ensemble[1] == 1) {
      #        colnames(df.train)[1] = c("trainy_") 
      rpart.train.00 <- rpart( trainy__ ~ ., data=df.trainxs, method=rpmethod, cp = 0 )
      rpart.train.01 <- prune(rpart.train.00, cp = 0.01)
      rpart.train.02 <- prune(rpart.train.00, cp = 0.02)
      rp_perf = rpart_perform(rpart.train.00 , df.testxs, testy__, family )
      rpart.devian.naive[1] = rp_perf[1] ; rpart.lincal.naive[1] = rp_perf[2] ; rpart.agree.naive[1] = rp_perf[3] ; rpart.nzero.naive[1] = rp_perf[4]
      rp_perf = rpart_perform(rpart.train.01 , df.testxs, testy__, family )
      rpart.devian.naive[2] = rp_perf[2] ; rpart.lincal.naive[2] = rp_perf[2] ; rpart.agree.naive[2] = rp_perf[3] ; rpart.nzero.naive[2] = rp_perf[4]
      rp_perf = rpart_perform(rpart.train.02 , df.testxs, testy__, family )
      rpart.devian.naive[3] = rp_perf[3] ; rpart.lincal.naive[3] = rp_perf[2] ; rpart.agree.naive[3] = rp_perf[3] ; rpart.nzero.naive[3] = rp_perf[4]
    }
      
    if (ensemble[2] == 1) { 
      rpart.trainF.00 <- rpart( trainy__ ~ ., data=df.trainxsO, method=rpmethod, cp = 0 )
      rpart.trainF.01 <- prune(rpart.trainF.00, cp = 0.01)
      rpart.trainF.02 <- prune(rpart.trainF.00, cp = 0.02)
      rp_perf = rpart_perform(rpart.trainF.00 , df.testxsO, testy__, family )
      rpart.devian.naive[4] = rp_perf[1] ; rpart.lincal.naive[4] = rp_perf[2] ; rpart.agree.naive[4] = rp_perf[3] ; rpart.nzero.naive[4] = rp_perf[4]
      rp_perf = rpart_perform(rpart.trainF.01 , df.testxsO, testy__, family )
      rpart.devian.naive[5] = rp_perf[2] ; rpart.lincal.naive[5] = rp_perf[2] ; rpart.agree.naive[5] = rp_perf[3] ; rpart.nzero.naive[5] = rp_perf[4]
      rp_perf = rpart_perform(rpart.trainF.02 , df.testxsO, testy__, family )
      rpart.devian.naive[6] = rp_perf[3] ; rpart.lincal.naive[6] = rp_perf[2] ; rpart.agree.naive[6] = rp_perf[3] ; rpart.nzero.naive[6] = rp_perf[4]
    } 
      
    if ((ensemble[3] ==1) & (!(family %in% c("binomial"))))  {
      xsnames = names(df.trainxsO)
      form1 = formula( paste( "trainy__ ~ ", paste(xsnames, collapse = " + " ), " + offset(ofst)" ) ) 
      rpart.trainO.00 <- rpart( form1 , data=df.trainxsO, method=rpmethod, cp = 0 )
      rpart.trainO.01 <- prune(rpart.trainO.00, cp = 0.01)
      rpart.trainO.02 <- prune(rpart.trainO.00, cp = 0.02)
      rp_perf = rpart_perform(rpart.trainO.00 , df.testxsO, testy__, family )
      rpart.devian.naive[7] = rp_perf[1] ; rpart.lincal.naive[7] = rp_perf[2] ; rpart.agree.naive[7] = rp_perf[3] ; rpart.nzero.naive[7] = rp_perf[4]
      rp_perf = rpart_perform(rpart.trainO.01 , df.testxsO, testy__, family )
      rpart.devian.naive[8] = rp_perf[2] ; rpart.lincal.naive[8] = rp_perf[2] ; rpart.agree.naive[8] = rp_perf[3] ; rpart.nzero.naive[8] = rp_perf[4]
      rp_perf = rpart_perform(rpart.trainO.02 , df.testxsO, testy__, family )
      rpart.devian.naive[9] = rp_perf[3] ; rpart.lincal.naive[9] = rp_perf[2] ; rpart.agree.naive[9] = rp_perf[3] ; rpart.nzero.naive[9] = rp_perf[4]
    }
    rpart.fit.00 = rpart.train.00 ; rpart.fitO.00 = rpart.trainO.00 ; rpart.fitF.00 = rpart.trainF.00 ;
    rpart.nzero = rpart.nzero.naive 
    if (track >= 1) { time_last = diff_time(time_start, time_last) } 
  }
  
  ##### STEPWISE fit ##########################################################################################    
  if (dostep == 1) { 
    if (track >= 1) { cat(paste0("\n", " ########## Initial stepwise model on all data ########################################" , "\n")) }
    cv.stepreg.fit.all  = cv.stepreg(xs, start, y_, event, steps_n, folds_n, method=method, family=family,foldid=foldid_,track=track)
    if (track >= 1) { time_last = diff_time(time_start, time_last) }
  }
  
  ##### AIC fit ##########################################################################################    
  if (doaic ==1) { 
    if (track >= 1) { cat(paste0("\n", " ########## Initial AIC model on all data #############################################" , "\n")) }
    if (dostep==1) { 
      func.fit.aic = aicreg(xs, start, y_, event, family=family, object=cv.stepreg.fit.all, track=track) 
    } else {
      func.fit.aic = aicreg(xs, start, y_, event, steps_n=steps_n, family=family, object=NULL, track=track) 
    }
    if (track >= 1) { time_last = diff_time(time_start, time_last) } 
  }
  
  ##### track time before outer nested loop ###############################################################
  if (track >= 2) { 
#    time_last = diff_time(time_start, time_last) 
    cat(paste0("\n", " ########## full data set analysis completed #####################################" , "\n")) 
  } 
  
  ##### FOLDS ##############################################################################################
  ##### FOLDS ##############################################################################################
  ## i_ = 1 
  if (do_ncv != 0) {
  for (i_ in 1:folds_n) {  
    ##=== begin folds ==========================================================
    if (track >= 1) { cat(paste0("\n", " ########## Entering Nested Cross Validation outer fold  ", i_, "  of  " , folds_n , "  ############################" , "\n")) }
  
    ##### set up train and test data sets in matrix form for glmnet & stepreg #####
    trainxs = xs[(foldid!=i_),]   
    testxs  = xs[(foldid==i_),]  
    test1   = c(rep(1,dim(testxs)[1])) 
    testxs1 = cbind(test1, testxs)                                              ## for models with Intercept 
    dim(trainxs) 
    dim(testxs) 
    
    if ( is.null(start) ) {
      trainstart = NULL
      teststart  = NULL
    } else {
      trainstart = start[(foldid!=i_)]    
      teststart  = start[(foldid==i_)]   
    }
  
    trainy_ = y_[(foldid!=i_)]    
    testy_  = y_[(foldid==i_)]   
    
    if (family == "cox") {
      trainevent = event[(foldid!=i_)] ;  
      testevent  = event[(foldid==i_)] ; 
    } else {
      trainevent = NULL
      testevent  = NULL 
    }    
    
    if (family=="cox") {
      if ( is.null(start) ) {
        train_data = data.frame(trainy_, trainevent, trainxs)
        test_data  = data.frame(testy_ , testevent , testxs )
      } else {
        train_data = data.frame(trainstart, trainy_, trainevent, trainxs)
        test_data  = data.frame(teststart , testy_ , testevent , testxs )
      }
    } else {
      train_data = data.frame(trainy_, trainxs)
      test_data  = data.frame(testy_ , testxs )
    }
    
    foldid_ = NULL
    ##### LASSO fit ##########################################################################################
    if (dolasso==1) { 
      if (track >= 1) { cat( " ## Starting LASSO model fits ##", "\n" ) } 

      if (relax) { 
     #  xs=trainxs; start=trainstart; y_=trainy_; event=trainevent; lambda=lambda; gamma=gamma; folds_n=folds_n; limit=limit; fine=fine; track=1; family=family ;         
        cv_glmnet_fit = cv.glmnetr( trainxs, trainstart, trainy_, trainevent, family=family, 
                                    lambda=lambda, gamma=gamma, folds_n=folds_n, limit=limit, fine=fine, 
                                    track=0, ties=ties, ... ) 
      } else {
        if (family=="cox") {   
          if ( is.null(start) ) { cv_glmnet_fit = cv.glmnet( trainxs, Surv(trainy_, trainevent)            , family="cox", lambda=lambda, relax=FALSE, ... ) 
          } else                { cv_glmnet_fit = cv.glmnet( trainxs, Surv(trainstart, trainy_, trainevent), family="cox", lambda=lambda, relax=FALSE, ... ) 
          }
        } else {
          cv_glmnet_fit = cv.glmnet( trainxs, trainy_, family=family, lambda=lambda, relax=FALSE, ... ) 
        }
      }
      
      if ((family == "cox") & (is.null(trainstart))) {
        cv_ridge_fit = cv.glmnet( trainxs, Surv(trainy_, trainevent), family=family, alpha=0, ... ) 
      } else if ((family == "cox") & (!is.null(trainstart))) {
        cv_ridge_fit = cv.glmnet( trainxs, Surv(trainstart, trainy_, trainevent), family=family, alpha=0, ... ) 
      } else {
        cv_ridge_fit = cv.glmnet( trainxs, trainy_, family=family, alpha=0, ... ) 
      }
      
      foldid_ = cv_glmnet_fit$foldid 
      
      lasso.nzero.cv[i_,1] = cv_glmnet_fit$nzero [ cv_glmnet_fit$index[2] ]
      lasso.nzero.cv[i_,2] = cv_glmnet_fit$nzero [ cv_glmnet_fit$index[1] ]     
      if (relax) {
        lasso.nzero.cv[i_,3] = cv_glmnet_fit$relaxed$nzero.1se 
        lasso.nzero.cv[i_,4] = cv_glmnet_fit$relaxed$nzero.min 
        lasso.nzero.cv[i_,5] = cv_glmnet_fit$nzero [ cv_glmnet_fit$relaxed$index.g0[2] ]
        lasso.nzero.cv[i_,6] = cv_glmnet_fit$nzero [ cv_glmnet_fit$relaxed$index.g0[1] ]
      }
      lassogammacv[i_,1]= cv_glmnet_fit$relaxed$gamma.1se 
      lassogammacv[i_,2]= cv_glmnet_fit$relaxed$gamma.min 
      
      #### LASSO predictions on TEST data ######################################
      # object=cv_glmnet_fit ; xs_new=testxs ; lam="lambda.min" ; gam="gamma.min" ;  comment=TRUE
      pred1se   = predict(cv_glmnet_fit, testxs, lam=cv_glmnet_fit$lambda.1se, gam=1)  ##### 1se lambda model predictions ################
      predmin   = predict(cv_glmnet_fit, testxs, lam=cv_glmnet_fit$lambda.min, gam=1)  ##### minimizing lambda model predictions #########
      pred1seR  = predict(cv_glmnet_fit, testxs, lam="lambda.1se" , gam="gamma.1se" )  ###### 1se RELAXED lasso ##########################
      predminR  = predict(cv_glmnet_fit, testxs, lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso ##########################
      pred1seR0 = predict(cv_glmnet_fit, testxs, lam=cv_glmnet_fit$relaxed$lambda.1se.g0, gam=0)  ###### 1se gamma = 0 RELAXED lasso #####
      predminR0 = predict(cv_glmnet_fit, testxs, lam=cv_glmnet_fit$relaxed$lambda.min.g0, gam=0)  ###### min gamma = 0 RELAXED lasso #####
      predridge = predict(cv_ridge_fit , testxs)  ###### min gamma = 0 RELAXED lasso #####
      
      if (family == "cox") {    
        if ( is.null(start) ) { yy = Surv(testy_, testevent) 
        } else {                yy = Surv(teststart, testy_, testevent) }
      } else { yy = testy_ }
      perfm1 = lasso_perf( yy , pred1se  , family )
      perfm2 = lasso_perf( yy , predmin  , family )
      perfm3 = lasso_perf( yy , pred1seR , family )
      perfm4 = lasso_perf( yy , predminR , family )
      perfm5 = lasso_perf( yy , pred1seR0, family )
      perfm6 = lasso_perf( yy , predminR0, family )
      perfm7 = lasso_perf( yy , predridge, family )
      lasso.devian.cv[i_,]  = c( perfm1[1] , perfm2[1] , perfm3[1] , perfm4[1] , perfm5[1] , perfm6[1] , perfm7[1] )
      lasso.cal.devian.cv[i_,]=c(perfm1[2] , perfm2[2] , perfm3[2] , perfm4[2] , perfm5[2] , perfm6[2] , perfm7[2] )    
      lasso.agree.cv[i_,]   = c( perfm1[3] , perfm2[3] , perfm3[3] , perfm4[3] , perfm5[3] , perfm6[3] , perfm7[3] )
      lasso.intcal.cv[i_,]  = c( perfm1[4] , perfm2[4] , perfm3[4] , perfm4[4] , perfm5[4] , perfm6[4] , perfm7[4] )
      lasso.lincal.cv[i_,]  = c( perfm1[5] , perfm2[5] , perfm3[5] , perfm4[5] , perfm5[5] , perfm6[5] , perfm7[5] )
  
      if (track >= 1) { time_last = diff_time(time_start, time_last) } 
    } 
    
    ##### XGBoost ############################################################################################
    if (doxgb_==1) { 
      if (track >= 1) { cat(paste0(" ## Starting XGBoost model fits ##" , "\n")) }

#      train.xgb.dat <- xgb.DMatrix(data = trainxs, label = trainy_)       

      if (family=="cox") {
        trainSurv.xgb = ifelse( trainevent == 1, trainy_, -trainy_)
        testSurv.xgb  = ifelse( testevent  == 1, testy_ , -testy_ )
      }
      
      if (ensemble[1]==1) {
        if (family=="cox") {
          train.xgb.dat <- xgb.DMatrix(data = trainxs, label = trainSurv.xgb)
          test.xgb.dat  <- xgb.DMatrix(data = testxs, label = testSurv.xgb)
        } else {
          train.xgb.dat <- xgb.DMatrix(data = trainxs, label = trainy_) 
          test.xgb.dat  <- xgb.DMatrix(data = testxs , label = testy_ ) 
        }
      }
      
      if ((ensemble[2] == 1)|(ensemble[3] == 1)) { 
        trainofst = predict(cv_glmnet_fit, trainxs, lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
        testofst  = predict(cv_glmnet_fit, testxs , lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
      }
      
      if (ensemble[2] == 1) { 
        if (family=="cox") {
          train.xgb.datF <- xgb.DMatrix(data = cbind(trainxs, ofst=trainofst), label = trainSurv.xgb)
          test.xgb.datF  <- xgb.DMatrix(data = cbind(testxs , ofst=testofst ), label = testSurv.xgb )
        } else {
          train.xgb.datF <- xgb.DMatrix(data = cbind(trainxs, ofst=trainofst), label = trainy_ ) 
          test.xgb.datF  <- xgb.DMatrix(data = cbind(testxs , ofst=testofst ), label = testy_  ) 
        }
      }
      
      if (ensemble[3] == 1) { 
        if (family=="cox") {
          train.xgb.datO <- xgb.DMatrix(data = trainxs, label = trainSurv.xgb, base_margin=trainofst)
          test.xgb.datO  <- xgb.DMatrix(data = testxs , label = testSurv.xgb , base_margin=testofst )
        } else {
          train.xgb.datO <- xgb.DMatrix(data = trainxs, label = trainy_ , base_margin=trainofst) 
          test.xgb.datO  <- xgb.DMatrix(data = testxs , label = testy_  , base_margin=testofst ) 
        }
      }

      if (family %in% c("cox", "binomial", "gaussian")) {
        
        if        (family == "cox"     ) { objective = "survival:cox" 
        } else if (family == "binomial") { objective = "binary:logistic" 
        } else                           { objective = "reg:squarederror" }
        
        ## SIMPLE model fit train ######################################################
        if (ensemble[1] == 1) {
          if (track >= 2) { print("XGB ensemble simple") }
          xgb.simple.train = xgb.simple(train.xgb.dat, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
          xgbperf = xgb_perform(xgb.simple.train, test.xgb.dat, y=testy_, evnt=testevent, family=family )
          xgb.devian.cv[i_,1] = xgbperf[1] ;  xgb.lincal.cv[i_,1] = xgbperf[2] ;  xgb.agree.cv [i_,1] = xgbperf[3] 
        }
        if (ensemble[2] == 1) {
          if (track >= 2) { print("XGB ensemble simple factor") 
              if (track >= 2) { 
                print( summary(trainofst) ) 
                print( summary(testofst ) ) }
              } 
          xgb.simple.trainF = xgb.simple(train.xgb.datF, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
          xgbperf = xgb_perform(xgb.simple.trainF, test.xgb.datF, y=testy_, evnt=testevent, family=family) 
          xgb.devian.cv[i_,2] = xgbperf[1] ;  xgb.lincal.cv[i_,2] = xgbperf[2] ;  xgb.agree.cv [i_,2] = xgbperf[3] 
        }
        if (ensemble[3] == 1) {
          if (track >= 2) { print("XGB ensemble simple offset") }
          xgb.simple.trainO = xgb.simple(train.xgb.datO, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
          xgbperf = xgb_perform(xgb.simple.trainO, test.xgb.datO, y=testy_, evnt=testevent, family=family) 
          xgb.devian.cv[i_,3] = xgbperf[1] ;  xgb.lincal.cv[i_,3] = xgbperf[2] ;  xgb.agree.cv [i_,3] = xgbperf[3] 
        }
        ## TUNED model fit train #######################################################
        if (ensemble[1] == 1) {
          if (track >= 2) { print("XGB ensemble tuned") }
          xgb.tuned.train = xgb.tuned(train.xgb.dat, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
          xgbperf = xgb_perform(xgb.tuned.train, test.xgb.dat, y=testy_, evnt=testevent, family=family) 
          xgb.devian.cv[i_,4] = xgbperf[1] ;  xgb.lincal.cv[i_,4] = xgbperf[2] ;  xgb.agree.cv [i_,4] = xgbperf[3] 
        }
        if (ensemble[2] == 1) {
          if (track >= 2) { print("XGB ensemble tuned factor") }
          xgb.tuned.trainF = xgb.tuned(train.xgb.datF, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
          xgbperf = xgb_perform(xgb.tuned.trainF, test.xgb.datF, y=testy_, evnt=testevent, family=family) 
          xgb.devian.cv[i_,5] = xgbperf[1] ;  xgb.lincal.cv[i_,5] = xgbperf[2] ;  xgb.agree.cv [i_,5] = xgbperf[3] 
        }
        if (ensemble[3] == 1) {
          if (track >= 2) { print("XGB ensemble tuned offset") }
          xgb.tuned.trainO = xgb.tuned(train.xgb.datO, objective = objective, folds_xgb=folds_xgb, nrounds=nrounds)
          xgbperf = xgb_perform(xgb.tuned.trainO, test.xgb.datO, y=testy_, evnt=testevent, family=family) 
          xgb.devian.cv[i_,6] = xgbperf[1] ;  xgb.lincal.cv[i_,6] = xgbperf[2] ;  xgb.agree.cv [i_,6] = xgbperf[3] 
        }
        ################################################################################
      }
      if (track >= 1) { time_last = diff_time(time_start, time_last) } 
      if (track >= 1) { cat(paste0(" (Continuing Nested Cross Validation outer fold  ", i_, "  of  " , folds_n , " )" , "\n")) }
    }

    ##### Neural Network ##########################################################################################
    if (doann_==1) { 
      
      if (track >= 1) { cat(" ## Starting Neural Network model fits ##" , "\n") }
      
      if (track == 1) { eppr_t = eppr ; eppr = min(eppr, -1) ;  eppr2_t = eppr2 ; eppr2 = min(eppr2, -1) ; }
      
        ## center about train mean 
        trainxs_means = colMeans(trainxs)
        trainxs_c     = sweep(trainxs, 2, trainxs_means, FUN = "-") 
        testxs_c      = sweep(testxs , 2, trainxs_means, FUN = "-") 
        ## get SD = 1 
        trainxs_sds   = sqrt( colSums(trainxs_c^2) / (dim(trainxs)[1]-1) ) 
        trainxs_sds_  = trainxs_sds 
        trainxs_sds_[(trainxs_sds < 1e-8)] = 1e-8 
        trainxs_z     = sweep(trainxs_c, 2, trainxs_sds_, FUN = "/") 
        trainxs_z[,(trainxs_sds < 1e-8)] = 0 
        testxs_z      = sweep(testxs_c , 2, trainxs_sds_, FUN = "/") 
        testxs_z[,(trainxs_sds < 1e-8)] = 0 
        trainxs_z0   = trainxs_z[, (trainxs_sds_ > 0) ]                         ## add lasso prediction as first column 
        testxs_z0    = testxs_z[, (trainxs_sds_ > 0) ]                           ## add lasso prediction as first column 
        
        # table(round(diag(cov(trainxs_z)),digits=8))
        if (sum(ensemble[c(2:6)]) > 0.5) {        
#          predminR      = predict(cv_glmnet_fit, testxs , lam="lambda.min" , gam="gamma.min" )  
          lassopredtrain = predict(cv_glmnet_fit, trainxs, lam="lambda.min" , gam="gamma.min" )
          lassopredtest  = predict(cv_glmnet_fit, testxs , lam="lambda.min" , gam="gamma.min" )
          lassobeta = predict(cv_glmnet_fit, lam="lambda.min" , gam="gamma.min" )
          if (family=="cox") {
            # trainy_, myevent=trainevent
            fit0 = coxph(Surv(trainy_, trainevent) ~ lassopredtrain) 
            fit0$loglik
            calbeta = fit0$coefficients
            lassopredtrain = calbeta * lassopredtrain 
            lassopredtest  = calbeta * lassopredtest 
          } else if (family %in% c("binomial", "gaussian")) {
            fit0 = glm(trainy_ ~ lassopredtrain, family=family) 
            beta = fit0$coefficients
            lassopredtrain = calbeta[1] + calbeta[2] * lassopredtrain 
            lassopredtest  = calbeta[1] + calbeta[2] * lassopredtest  
          }
          
##          concordance(y_ ~ lassopred)
##          summary(1/(1+exp(-lassopred)))
##          (lassobeta[[1]] == 0)
          trainxs_z1 = trainxs_z[,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])]]           ## pick up non zero features, remove intercept from this list 
          trainxs_z2 = cbind(lasso=lassopredtrain,trainxs_z1)                                ## add lasso prediction as first column 
          testxs_z1  = testxs_z[,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])]]            ## pick up non zero features, remove intercept from this list 
          testxs_z2  = cbind(lasso=lassopredtest,testxs_z1)                                  ## add lasso prediction as first column 
          # table(round(diag(cov(trainxs_z2)),digits=8))
        }

      if (getwd  == 1) { wd  = cv_ridge_fit$lambda.min } 
      if (getwd2 == 1) { wd2 = cv_ridge_fit$lambda.min }

      if (getl1  == 1) { l1  = cv_glmnet_fit$relaxed$lambda.min.g0 } 
      if (getl12 == 1) { l12 = cv_glmnet_fit$relaxed$lambda.min.g0 }
        
      if (family %in% c("cox", "binomial", "gaussian")) {
        if (ensemble[1] == 1) { 
          if (eppr >= -2) { cat(paste0("\n  ** fitting uninformed ANN **\n")) }
          ann_fit_1 = ann_tab_cv(myxs=trainxs_z0, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=fold_n, family=family, 
                                 epochs=epochs, eppr=eppr, lenz1=lenz1, lenz2=lenz2, mylr=mylr, actv=actv, drpot=drpot, wd=wd, l1=l1, minloss=minloss, gotoend=gotoend) 
          perfm1 = ann_perform(ann_fit_1, testxs_z0, testy_, family, teststart, testevent)
        } else { perfm1 = c(10,10,0,0,0) }
        if (ensemble[4] == 1) { 
          if (eppr >= -2) { cat(paste0("  ** fitting ANN informed on relaxed lasso terms **\n")) }
          ## limit to lasso variables 
          ann_fit_2 = ann_tab_cv(myxs=trainxs_z1, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=fold_n, family=family,
                                 epochs=epochs, eppr=eppr, lenz1=lenz1, lenz2=lenz2, mylr=mylr, actv=actv, drpot=drpot, wd=wd, l1=l1, minloss=minloss, gotoend=gotoend) 
          perfm2 = ann_perform(ann_fit_2, testxs_z1, testy_, family, teststart, testevent)
        } else { perfm2 = c(10,10,0,0,0) }
        if (ensemble[2] == 1) {
          if (eppr >= -2) { cat(paste0("  ** fitting ANN with feature X*Beta and relaxed lasso terms **\n")) } 
          lenz1_ = lenz1 + 1 ; lenz2_ = lenz2 + 1 ; 
          ann_fit_3 = ann_tab_cv(myxs=trainxs_z2, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=fold_n, family=family, 
                                 epochs=epochs, eppr=eppr, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr, actv=actv, drpot=drpot, wd=wd, l1=l1, minloss=minloss, gotoend=gotoend) 
          perfm3 = ann_perform(ann_fit_3, testxs_z2, testy_, family, teststart, testevent)
        } else { perfm3 = c(10,10,0,0,0) }
        if (ensemble[3] == 1) { 
          if (eppr >= -2) { cat(paste0("  ** fitting ANN including X*Beta from relaxed lasso in initial weight structure  **\n")) } 
          lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
          ann_fit_4 = ann_tab_cv(myxs=trainxs_z2, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=fold_n, family=family, 
                                 epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, drpot=drpot, wd=wd, l1=l1, lasso=1, adjustw=0, minloss=minloss, gotoend=gotoend) 
          perfm4 = ann_perform(ann_fit_4, testxs_z2, testy_, family, teststart, testevent)
        } else { perfm4 = c(10,10,0,0,0) }
        if (ensemble[5] == 1) { 
          if (eppr >= -2) { cat(paste0("  ** fitting ANN including X*Beta from relaxed lasso in weight structure each epoch **\n")) } 
          lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
          ann_fit_5 = ann_tab_cv(myxs=trainxs_z2, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=fold_n, family=family, 
                                 epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, drpot=drpot, wd=wd2, l1=l12, lasso=1, adjustw=1, minloss=minloss, gotoend=gotoend) 
          perfm5 = ann_perform(ann_fit_5, testxs_z2, testy_, family, teststart, testevent)
        } else { perfm5 = c(10,10,0,0,0) }
        if (ensemble[6] == 9) { 
          if (eppr >= -2) { cat(paste0("  ** fitting ANN including X*Beta from relaxed lasso in weight structure each epoch **\n")) } 
          ann_fit_6 = ann_tab_cv(myxs=trainxs_z1, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=fold_n, family=family, 
                                 epochs=epochs2, eppr=eppr2, lenz1=lenz1, lenz2=lenz2, mylr=mylr2, drpot=drpot, wd=wd2, l1=l12, minloss=minloss, gotoend=gotoend) 
          perfm6 = ann_perform(ann_fit_6, testxs_z1, testy_, family, teststart, testevent)  ## offset=lassopred, 
        } else { perfm6 = c(10,10,0,0,0) }
        ann.devian.cv[i_,]  = c( perfm1[1] , perfm2[1] , perfm3[1] , perfm4[1] , perfm5[1] , perfm6[1] )
        ann.cal.devian.cv[i_,]=c(perfm1[2] , perfm2[2] , perfm3[2] , perfm4[2] , perfm5[2] , perfm6[2] )    
        ann.agree.cv [i_,]  = c( perfm1[3] , perfm2[3] , perfm3[3] , perfm4[3] , perfm5[3] , perfm6[3] )
        ann.intcal.cv[i_,]  = c( perfm1[4] , perfm2[4] , perfm3[4] , perfm4[4] , perfm5[4] , perfm6[4] )
        ann.lincal.cv[i_,]  = c( perfm1[5] , perfm2[5] , perfm3[5] , perfm4[5] , perfm5[5] , perfm6[5] )
      } 
      if (track == 1) { eppr = eppr_t ; eppr2 = eppr2_t }
      if (track >= 1) { time_last = diff_time(time_start, time_last) } 
    } 
    
    ##### RPART fit ##########################################################################################
    if (dorpart==1) { 
      if (track >= 1) { cat(" ## Starting RPART model fits ##" , "\n") }
    
#<--------------------------------------------------------------------------        
#      rpart.devian.naive = c(0,0,0, 0,0,0, 0,0,0)
#      rpart.nzero.naive  = c(0,0,0, 0,0,0, 0,0,0)
#      rpart.lincal.naive = c(0,0,0, 0,0,0, 0,0,0)
#      rpart.agree.naive  = c(0,0,0, 0,0,0, 0,0,0) 
#      nms = c("cp=0.00", "cp=0.01", "cp=0.02", "F cp=0.00", "F cp=0.01", "F cp=0.02", "O cp=0.00", "O cp=0.01", "O cp=0.02" ) 
#      names(rpart.devian.naive)  = nms 
#      names(rpart.nzero.naive )  = nms
#      names(rpart.lincal.naive)  = nms
#      names(rpart.agree.naive )  = nms
#     
#      trainxs = xs 
#      testxs = xs 
#      trainstart = start 
#      teststart = start 
#      trainy_ = y_ 
#      testy_  = y_ 
#      trainevent = event 
#      testevent  = event 
#<--------------------------------------------------------------------------

      df.trainxs = as.data.frame(trainxs) 
      df.testxs  = as.data.frame(testxs) 
      
      if ((ensemble[2] == 1) | (ensemble[3] == 1)) { 
         trainofst = predict(cv_glmnet_fit  , trainxs, lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
        #trainofst = predict(cv_glmnet_fit_f, trainxs, lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
        length(trainofst) ; 
         testofst  = predict(cv_glmnet_fit  , testxs , lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
        #testofst  = predict(cv_glmnet_fit_f, testxs , lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
        length(testofst) ; 
        trainxsO  = cbind(trainxs,ofst=trainofst) 
        dim(xs) ; dim(trainxs) ; dim(trainxsO) ; length(trainofst) ; 
        testxsO   = cbind(testxs ,ofst=testofst) 
        df.trainxsO = as.data.frame(trainxsO)
        df.testxsO  = as.data.frame(testxsO )
        #      dim(df.testxsO) ; dim(df.trainxsO)
      }
      
      if (family=="cox") {
        if (is.null(start)) { 
          trainy__ = Surv(trainy_, trainevent)
          testy__  = Surv(testy_ , testevent)
          if ((ensemble[2] ==1) | (ensemble[3] ==1)) { df.trainO   = as.data.frame(cbind(trainy_, trainevent, trainxsO) ) }
        } else { 
          trainy__ = Surv(trainstart, trainy_, trainevent) 
          testy__  = Surv(teststart , testy_ , testevent ) 
          if ((ensemble[2] ==1) | (ensemble[3] ==1)) { df.trainO   = as.data.frame(cbind(trainstart, trainy_, trainevent, trainxsO) ) }
        }
      } else {
        trainy__ = trainy_
        testy__  = testy_ 
      }
      
      if (family == "cox") {
        rpmethod = "exp" 
      } else if (family == "binomial") {
        rpmethod = "class" 
      } else if (family == "gaussian") {
        rpmethod = "anova"
      }
      
      if (ensemble[1] == 1) {
        #        colnames(df.train)[1] = c("trainy_") 
        rpart.train.00 <- rpart( trainy__ ~ ., data=df.trainxs, method=rpmethod, cp = 0 )
        rpart.train.01 <- prune(rpart.train.00, cp = 0.01)
        rpart.train.02 <- prune(rpart.train.00, cp = 0.02)
        rp_perf = rpart_perform(rpart.train.00 , df.testxs, testy__, family )
        rpart.devian.cv[i_,1] = rp_perf[1] ; rpart.lincal.cv[i_,1] = rp_perf[2] ; rpart.agree.cv[i_,1] = rp_perf[3] ; rpart.nzero.cv[i_,1] = rp_perf[4]
        rp_perf = rpart_perform(rpart.train.01 , df.testxs, testy__, family )
        rpart.devian.cv[i_,2] = rp_perf[2] ; rpart.lincal.cv[i_,2] = rp_perf[2] ; rpart.agree.cv[i_,2] = rp_perf[3] ; rpart.nzero.cv[i_,2] = rp_perf[4]
        rp_perf = rpart_perform(rpart.train.02 , df.testxs, testy__, family )
        rpart.devian.cv[i_,3] = rp_perf[3] ; rpart.lincal.cv[i_,3] = rp_perf[2] ; rpart.agree.cv[i_,3] = rp_perf[3] ; rpart.nzero.cv[i_,3] = rp_perf[4]
      }
      
      if (ensemble[2] == 1) { 
        rpart.trainF.00 <- rpart( trainy__ ~ ., data=df.trainxsO, method=rpmethod, cp = 0 )
        rpart.trainF.01 <- prune(rpart.trainF.00, cp = 0.01)
        rpart.trainF.02 <- prune(rpart.trainF.00, cp = 0.02)
        rp_perf = rpart_perform(rpart.trainF.00 , df.testxsO, testy__, family )
        rpart.devian.cv[i_,4] = rp_perf[1] ; rpart.lincal.cv[i_,4] = rp_perf[2] ; rpart.agree.cv[i_,4] = rp_perf[3] ; rpart.nzero.cv[i_,4] = rp_perf[4]
        rp_perf = rpart_perform(rpart.trainF.01 , df.testxsO, testy__, family )
        rpart.devian.cv[i_,5] = rp_perf[2] ; rpart.lincal.cv[i_,5] = rp_perf[2] ; rpart.agree.cv[i_,5] = rp_perf[3] ; rpart.nzero.cv[i_,5] = rp_perf[4]
        rp_perf = rpart_perform(rpart.trainF.02 , df.testxsO, testy__, family )
        rpart.devian.cv[i_,6] = rp_perf[3] ; rpart.lincal.cv[i_,6] = rp_perf[2] ; rpart.agree.cv[i_,6] = rp_perf[3] ; rpart.nzero.cv[i_,6] = rp_perf[4]
      } 
      
      if ((ensemble[3] ==1) & (!(family %in% c("binomial"))))  {
        xsnames = names(df.trainxsO)
        form1 = formula( paste( "trainy__ ~ ", paste(xsnames, collapse = " + " ), " + offset(ofst)" ) ) 
        rpart.trainO.00 <- rpart( form1 , data=df.trainxsO, method=rpmethod, cp = 0 )
        rpart.trainO.01 <- prune(rpart.trainO.00, cp = 0.01)
        rpart.trainO.02 <- prune(rpart.trainO.00, cp = 0.02)
        rp_perf = rpart_perform(rpart.trainO.00 , df.testxsO, testy__, family )
        rpart.devian.cv[i_,7] = rp_perf[1] ; rpart.lincal.cv[i_,7] = rp_perf[2] ; rpart.agree.cv[i_,7] = rp_perf[3] ; rpart.nzero.cv[i_,7] = rp_perf[4]
        rp_perf = rpart_perform(rpart.trainO.01 , df.testxsO, testy__, family )
        rpart.devian.cv[i_,8] = rp_perf[2] ; rpart.lincal.cv[i_,8] = rp_perf[2] ; rpart.agree.cv[i_,8] = rp_perf[3] ; rpart.nzero.cv[i_,8] = rp_perf[4]
        rp_perf = rpart_perform(rpart.trainO.02 , df.testxsO, testy__, family )
        rpart.devian.cv[i_,9] = rp_perf[3] ; rpart.lincal.cv[i_,9] = rp_perf[2] ; rpart.agree.cv[i_,9] = rp_perf[3] ; rpart.nzero.cv[i_,9] = rp_perf[4]
      }
      rpart.nzero = rpart.nzero.naive 
      if (track >= 1) { time_last = diff_time(time_start, time_last) } 
    }
    
    ###### STEPWISE fits #####################################################################################
    if (dostep == 1) { 
      if (track >= 1) { cat(paste0(" ## Starting STEPWISE model fits ## ", "\n")) }
      ## xs_cv=trainxs ; start_cv=trainstart ; y_cv=trainy_ ; event_cv=trainevent ; steps_n=steps_n ; folds_n_cv=folds_n ; method=method ; family=family ; 
      cv.stepreg.fit = cv.stepreg(trainxs, trainstart, trainy_, trainevent, steps_n, folds_n, method=method, family=family, foldid=foldid_, track=0) 
      stepreg.fit.all.best = cv.stepreg.fit$stepreg.fit.all.best 
#      names(cv.stepreg.fit) ; cv.stepreg.fit$best.p ; cv.stepreg.fit$df.p ; cv.stepreg.fit$pvalue
#            names(cv.stepreg.fit) ; cv.stepreg.fit$best.p ; cv.stepreg.fit$func.fit.p      
#      str(cv.stepreg.fit)
#      cat(paste0( cv.stepreg.fit$stepreg.fit.all.best[1,1] , cv.stepreg.fit$stepreg.fit.all.best[2,2], cv.stepreg.fit$stepreg.fit.all.best[3,3], cv.stepreg.fit$stepreg.fit.all.best[4,4])) 

#     stepreg.fit.all.best[1:10, 1:7] 
      #### stepwise tuned by DF ################################################
#      df = cv.stepreg.fit$best.df
#      step.devian.cv [i_,4] = stepreg.fit.all.best [ df , 4 ]        
#      step.devian.cv [i_,1] = stepreg.fit.all.best [ df , 5 ]  
#      step.agree.cv  [i_,1] = stepreg.fit.all.best [ df , 6 ] 
#      step_df_cv    [i_,1] = df 

      func.fit.df = cv.stepreg.fit$func.fit.df                                  ##  cox model suggested by Cross Validation 
      mod_df      = length(func.fit.df$coefficients)                            ##  Number of terms in Cross Validation model 
      testxb  = predict( func.fit.df, as.data.frame( testxs) )                  ##  predicted XB scores for test data based upon train data model 
      if (family == "cox") {
        if (is.null(start)) { 
          fit0 = coxph( Surv(testy_, testevent) ~ rep(1,length(testxb)), init=c(1), control=coxph.control(iter.max=0))
          fit3 = coxph( Surv(testy_, testevent) ~ testxb, init=c(1))
        } else { 
          fit0 = coxph( Surv(testy_, teststart, testevent) ~ rep(1,length(testxb)), init=c(1), control=coxph.control(iter.max=0))
          fit3 = coxph( Surv(testy_, teststart, testevent) ~ testxb, init=c(1))  
        }
        step.devian.cv[i_,4] = -2*fit0$loglik[1] / fit0$nevent  
        step.devian.cv[i_,1] = -2*fit3$loglik[1] / fit3$nevent  
        step.lincal.cv[i_,1] = fit3$coefficients[1]
        step.agree.cv [i_,1] = fit3$concordance[6]   
        step_df_cv   [i_,1] = mod_df 
      } else if (family == "binomial") {      
        fit1 = glm( testy_ ~ testxb, family=family )  
        p_ = sum(testy_)/length(testy_) 
        step.devian.cv[i_,4] = -2 * sum( log(p_)*testy_ + log(1-p_)*(1-testy_) ) / length(testy_)
        p_ = 1/(1+exp(-testxb)) 
        step.devian.cv[i_,1] = -2 * ( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_)
        step.lincal.cv[i_,1] = fit1$coefficients[2] 
        step.agree.cv [i_,1] = concordance(testy_ ~ testxb)[[1]]
        step_df_cv   [i_,1] = mod_df - 1 
      }  else if (family == "gaussian") {  
        step.devian.cv[i_,4] = (var(testy_) * (length(testy_)-1)) / length(testy_) 
        step.devian.cv[i_,1] = sum((testy_-testxb)^2)  / length(testxb) 
        step.lincal.cv[i_,1] = cov(testy_, testxb)/var(testxb) 
        step.agree.cv [i_,1] = cor(x=testy_,y=testxb)^2         
        step_df_cv    [i_,1] = mod_df - 1 
      }
      
      #### stepwise tuned by p (critical value) ################################
      func.fit.p = cv.stepreg.fit$func.fit.p                                  ##  cox model suggested by Cross Validation 
      mod_df     = length(func.fit.p$coefficients)                            ##  Number of terms in Cross Validation model 
      testxb  = predict( func.fit.p, as.data.frame( testxs) )              ##  predicted scores for test data based upon train data model 
      if (family == "cox") {
        if (is.null(start)) { 
          fit3 = coxph( Surv(testy_, testevent) ~ testxb, init=c(1))
        } else { 
          fit3 = coxph( Surv(testy_, teststart, testevent) ~ testxb, init=c(1))  
        }
        step.devian.cv[i_,2] = -2*fit3$loglik[1] / fit3$nevent  
        step.lincal.cv [i_,2] = fit3$coefficients[1] 
        step.agree.cv [i_,2] = fit3$concordance[6]   
        step_df_cv   [i_,2] = mod_df 
      } else if (family == "binomial") {   
        fit1 = glm( testy_ ~ testxb, family=family )  
        p_ = 1/(1+exp(-testxb)) 
        step.devian.cv[i_,2] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_)
        step.lincal.cv[i_,2] = fit1$coefficients[2] 
        step.agree.cv [i_,2] = concordance(testy_ ~ testxb)[[1]]
        step_df_cv   [i_,2] = mod_df - 1 
      }  else if (family == "gaussian") {  
        step.devian.cv[i_,2] = sum((testy_ - testxb)^2) / length(testxb) 
        step.lincal.cv[i_,2] = cov(testy_, testxb) / var(testxb) 
        step.agree.cv [i_,2] = cor(testy_, testxb)^2         
        step_df_cv   [i_,2] = mod_df - 1 
      }
      step_p_cv    [i_,1] = cv.stepreg.fit$best.p 
      if (track >= 1) { time_last = diff_time(time_start, time_last) } 
    }
    
    ###### AIC fits ##########################################################################################
    if (doaic==1){
      if (track >= 1) { cat(" ## Starting AIC model fits ##" , "\n") }
      if (dostep==1) {
        stepreg.fit.best = as.matrix( cv.stepreg.fit$stepreg.fit.all.best ) 
#        class(stepreg.fit) 
      } else {
        stepreg.fit = stepreg(trainxs, trainstart, trainy_, trainevent, steps_n=steps_n)  
#        class(stepreg.fit) 
        class(stepreg.fit) = "data.frame" 
        stepreg.fit.best = stepreg.fit[(stepreg.fit$best==1) , ] 
#        class( stepreg.fit.best )
      }
#      stepreg.fit.best[1:10,1:7]
      aic = 2*stepreg.fit.best[,1] - 2*stepreg.fit.best[,5]                     
      mod_df  = which.min( aic )                                                    ## may be differnt from first p > 0.15 XX
   
      if (family=="cox") {
        beta = stepreg.fit.best [ mod_df , (xs_ncol+8):(2*xs_ncol+7) ]              ## get beta, the regression estimate
#        class(beta)
#        beta = matrix(beta, nrow=xs_ncol, ncol=1) 
      } else {
        beta = stepreg.fit.best [ mod_df , (xs_ncol+8):(2*xs_ncol+8) ]              ## get beta, the regression estimate
#        class(beta)
#        beta = matrix(beta, nrow=(xs_ncol+1), ncol=1) 
      }
      
      beta = as.numeric(beta)
      if (family=="cox") { testxb  = testxs  %*% beta
      } else             { testxb  = testxs1 %*% beta }                      ## get predicteds 
      
      if (family=="cox") {
        if (is.null(start)) { SURV = Surv(           testy_, testevent)
        } else {              SURV = Surv(teststart, testy_, testevent) }
        fit0  = coxph( SURV ~ rep(1,length(testxb)), control=coxph.control(iter.max=0) )  
        fit1  = coxph( SURV ~ testxb, init=c(1) )  
        step.devian.cv[i_,4] = -2*fit0$loglik[1] / fit0$nevent  
        
        step.devian.cv[i_,3] = -2*fit1$loglik[1] / fit1$nevent
        step.lincal.cv[i_,3] = fit1$coefficients[1] 
        step.agree.cv [i_,3] = fit1$concordance[6]   
        step_df_cv   [i_,3] = mod_df
      } else if (family=="binomial") {
        fit1 = glm( testy_ ~ testxb, family=family )  
        if (dostep!=1) {
          p_ = sum(testy_)/length(testy_) 
          step.devian.cv[i_,4] = -2*sum( log(p_)*testy_ + log(1-p_)*(1-testy_) )
        }
        p_ = 1/(1+exp(-testxb)) 
        step.devian.cv[i_,3] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testxb)  
        step.lincal.cv[i_,3] = fit1$coefficients[2] 
        step.agree.cv [i_,3] = concordance(testy_ ~ testxb)[[1]]   
        step_df_cv   [i_,3] = mod_df - 1 
      } else if (family=="gaussian") {
        if (dostep!=1) {
          step.devian.cv[i_,4] = (var(testy_) * (length(testy_)-1)) / length(testy_)
        }
        step.devian.cv[i_,3] = sum((testy_ - testxb)^2)  / length(testxb) 
        step.lincal.cv[i_,3] = cov(testy_, testxb) / var(testxb) 
        step.agree.cv [i_,3] = cor(testy_, testxb)^2 
        step_df_cv   [i_,3] = mod_df - 1 
      } 
      if (track >= 1) { time_last = diff_time(time_start, time_last) } 
    }
    ## end within fold AIC FIT #################################################
  }
  }
  ##### END FOLDS ##########################################################################################
  ##### END FOLDS ##########################################################################################

  if (track >= 1) { cat(paste0("\n", " ################################################################################################" , "\n")) } 

#  row.names(c.naive) = NULL
  if        (family=="cox"     ) { nevents=sum(event) 
                                   temp = coxph(Surv(y_, event) ~ 1, ties=ties)
                                   null.deviance = -2*(temp$loglik[1])/nevents 
  } else if (family=="binomial") { nevents=sum(y_) 
                                   p = nevents/length(y_) 
                                   null.deviance = -2 * ( p*log(p) + (1-p)*log(1-p) ) 
  } else                         { nevents=NA          
                                   null.deviance = (length(y_)-1)*var(y_)/length(y_) 
  } 
  
  nestedcv = list( sample=c(family=family, n=dim(xs)[1], nevents=nevents, xs.columns=dim(xs)[2], xs.df=rankMatrix(xs)[1], null.dev.n=null.deviance), 
                   tuning=c(steps_n=steps_n,folds_n=folds_n,method=method,dolasso=dolasso,doann_=doann_,doxgb_=doxgb_,dorpart=dorpart,dostep=dostep,doaic=doaic,
                   seed=seed[1], tseed=seed[2], ties=ties, limit=limit),
                   ensemble = ensemble, 
                   doann = doann, 
                   doxgb = doxgb, 
                   do_ncv=do_ncv, 
                   foldid=foldid ) 

  if (doann_ == 1) {
                   nestedcv$ann.devian.cv = ann.devian.cv
                   nestedcv$ann.intcal.cv = ann.intcal.cv
                   nestedcv$ann.lincal.cv = ann.lincal.cv
                   nestedcv$ann.agree.cv  = ann.agree.cv
                   nestedcv$ann.devian.naive = ann.devian.naive
                   nestedcv$ann.agree.naive  = ann.agree.naive 
                   nestedcv$ann.intcal.naive = ann.intcal.naive
                   nestedcv$ann.lincal.naive = ann.lincal.naive
                   nestedcv$ann_zb         = ann_zb 
                   nestedcv$ann_fit_1   = ann_fit_1_f 
                   nestedcv$ann_fit_2   = ann_fit_2_f 
                   nestedcv$ann_fit_3   = ann_fit_3_f 
                   nestedcv$ann_fit_4   = ann_fit_4_f
                   nestedcv$ann_fit_5   = ann_fit_5_f
                   nestedcv$ann_fit_6   = ann_fit_6_f 
                   } 
  
  if (dolasso == 1) { nestedcv$lasso.nzero.cv  = lasso.nzero.cv 
                      nestedcv$lasso.devian.cv = lasso.devian.cv
                      nestedcv$lasso.cal.devian.cv = lasso.cal.devian.cv
                      nestedcv$lasso.intcal.cv = lasso.intcal.cv
                      nestedcv$lasso.lincal.cv = lasso.lincal.cv
                      nestedcv$lasso.agree.cv  = lasso.agree.cv
                      nestedcv$lasso.agree.naive=lasso.agree.naive 
                      nestedcv$cv_glmnet_fit   = cv_glmnet_fit_f 
                      nestedcv$cv_ridge_fit    = cv_ridge_fit_f } 
  
  if (doxgb_ == 1)   { nestedcv$xgb.devian.cv    = xgb.devian.cv
                      nestedcv$xgb.lincal.cv    = xgb.lincal.cv
                      nestedcv$xgb.agree.cv     = xgb.agree.cv
                      nestedcv$xgb.devian.naive = xgb.devian.naive
                      nestedcv$xgb.lincal.naive = xgb.lincal.naive 
                      nestedcv$xgb.agree.naive  = xgb.agree.naive 
                      nestedcv$xgb.simple.fit   = xgb.simple.fit
                      if (ensemble[2] == 1) { nestedcv$xgb.simple.fitF = xgb.simple.fitF }
                      if (ensemble[3] == 1) { nestedcv$xgb.simple.fitO = xgb.simple.fitO }
                      nestedcv$xgb.tuned.fit    = xgb.tuned.fit 
                      if (ensemble[2] == 1) { nestedcv$xgb.tuned.fitF = xgb.tuned.fitF }
                      if (ensemble[3] == 1) { nestedcv$xgb.tuned.fitO = xgb.tuned.fitO }
                    }   

  if (dorpart == 1) { nestedcv$rpart.nzero.cv  = rpart.nzero.cv 
                      nestedcv$rpart.devian.cv = rpart.devian.cv
                      nestedcv$rpart.lincal.cv = rpart.lincal.cv
                      nestedcv$rpart.agree.cv  = rpart.agree.cv 
                      nestedcv$rpart.nzero     = rpart.nzero
                      nestedcv$rpart.agree.naive = rpart.agree.naive
                      nestedcv$rpart.fit.00    = rpart.fit.00
                      if (ensemble[2] == 1) { 
                        nestedcv$rpart.fitF.00 = rpart.fitF.00 
                      }
                      if ((ensemble[3] == 1) & (family %in% c("cox"))) { 
                        nestedcv$rpart.fitO.00 = rpart.fitO.00 
                      }
                    } 

  if ((dostep == 1) | (doaic==1)) { nestedcv$step.devian.cv = step.devian.cv 
                                    nestedcv$step.lincal.cv = step.lincal.cv 
                                    nestedcv$step.agree.cv  = step.agree.cv
                                    nestedcv$step_df_cv     = step_df_cv 
                                    nestedcv$step_p_cv      = step_p_cv }
  
  if (dostep  == 1) { nestedcv$cv.stepreg.fit = cv.stepreg.fit.all }
  
  if (doaic   == 1) { nestedcv$func.fit.aic = func.fit.aic } 

  class(nestedcv) <- c("nested.glmnetr")

  names (nestedcv)

  return( nestedcv )
  cat("\n Program End \n")
  if (track >= 2) { time_last = diff_time(time_start, time_last) } 
}

###############################################################################################################################################
###############################################################################################################################################
