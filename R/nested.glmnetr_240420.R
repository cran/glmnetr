################################################################################
##### nested.glmnetr_yymmdd ####################################################
################################################################################
#' Using (nested) cross validation, describe and compare some machine learning model performances
#' 
#' @description Performs a nested cross validation for cross validation informed 
#' relaxed lasso, Gradient Boosting Machine (GBM), Random Forest (RF), (artificial) 
#' Neural Network (ANN) with two hidden layers, Recursive Partitioning (RPART) 
#' and step wise regression.  That is
#' hyper parameters for all these models are informed by cross validation (CV) (or in the 
#' case of RF by out-of-bag calculations), and a second  layer of CV 
#' (or analogously for the RF) is 
#' used to evaluate the performance of these CV informed model fits.  For
#' step wise regression CV is used to inform either a p-value for entry or degrees of 
#' freedom (df) for the final model choice.  For input 
#' we require predictors (features) to be in numeric matrix format with no missing 
#' values.  This is similar to how the glmnet package expects predictors.  For 
#' survival data we allow input of start time as an option, and require stop time, 
#' and an event indicator, 1 for event and 0 for censoring, as separate terms. This
#' may seem unorthodox as it might seem simpler to accept a Surv() object as input.  However,
#' multiple packages we use for model fitting models require data
#' in various formats and this choice was the most straight forward for constructing the 
#' data formats required.  As an example, 
#' the XGBoost routines require a data format specific to the XGBoost 
#' package, not a matrix, not a data frame.  Note, for XGBoost and survival models, 
#' only a "stop time" variable, taking a positive value to indicate 
#' being associated with an event, and the negative of the time when 
#' associated with a censoring, is passed to the input data object for 
#' analysis.  
#' 
#' @param xs     predictor input - an n by p matrix, where n (rows) is sample 
#' size, and p (columns) the number of predictors.  Must be in matrix form for 
#' complete data, no NA's, no Inf's, etc., and not a data frame. 
#' @param start  optional start times in case of a Cox model.  A numeric (vector) 
#' of length same as number of patients (n).  Optionally start may be specified as 
#' a column matrix in which case the colname value is used when outputing summaries.
#' @param y_     dependent variable as a vector: time, or stop time for Cox 
#' model, Y_ 0 or 1 for binomal (logistic), numeric for gaussian. Must be a 
#' vector of length same as number of sample size. Optionally y_ may be specified as 
#' a column matrix in which case the colname value is used when outputing summaries.
#' @param event   event indicator, 1 for event, 0 for census, Cox model only.
#' Must be a numeric vector of length same as sample size.  Optionally event may be specified as 
#' a column matrix in which case the colname value is used when outputing summaries.
#' @param family  model family, "cox", "binomial" or "gaussian" (default) 
#' @param folds_n the number of folds for the outer loop of the nested cross 
#' validation, and if not overridden by the individual model specifications, also 
#' the number of folds for the inner loop of the nested cross validation, i.e. 
#' the number of folds used in model derivation.   
#' @param stratified 1 to generate fold IDs stratified on outcome or event 
#' indicators for the binomial or Cox model.
#' @param do_ncv 1 by default to do the Nested Cross Validation, or 0
#' to only fit the various models without doing the Nested part. In this case the 
#' nested.glmnetr() function will only derive the models based upon the full 
#' data set.  This may be useful when exploring various models without having to 
#' the Nested Cross Validation for assessing model performance. for example when
#' wanting to examine a "lasso informed" extreme gradient boosting 
#' models (GBM) or Artificial Neural Network (ANN) models which are based upon 
#' both a lasso fit and a GBM or ANN fit. See the predict_ann_tab() function 
#' regarding getting predicteds for this "lasso informed" models from a 
#' nested.glmnetr() output.  
#' @param dolasso fit and do cross validation for lasso model, 0 or 1
#' @param doxgb fit and evaluate a cross validation informed XGBoost (GBM) 
#' model.  1 for yes, 0 for no (default).  By default the number of folds used when 
#' training the GBM model will be the same as the number of folds used in the outer
#' loop of the nested cross validation, and the maximum number of rounds when 
#' training the GBM model is set to 1000.  To control these values one 
#' may specify a list for the doxgb argument.  The list can have 
#' elements  $nfold, $nrounds,
#' and $early_stopping_rounds, each numerical values of length 1, $folds, a list as 
#' used by xgb.cv() do identify folds for cross validation, and $eta, $gamma, $max_depth, 
#' $min_child_seight, $colsample_bytree, $lambda, $alpha and $subsample, each a numeric 
#' of length 2 giving the lower and upper values for the respective tuning 
#' parameter.  Here we deviate from nomenclature used elsewhere in the package to
#' be able to use terms those used in the 'xgboost' (and mlrMBO) package, in particular as used
#' in xgb.train(), e.g. nfold instead of folds_n and folds instead of foldid.  If 
#' not provided defaults will be used.  Defaults
#' can be seen from the output object$doxgb element, again a list. In case not NULL, 
#' the seed and folds option values override the $seed and $folds values.  
#' 
#' If to shorten run time the user sets nfold to a value 
#' other than folds_n we recommend that nfold = folds_n/2 or folds_n/3.  Then the 
#' folds will be formed by collapsing the folds_n folds allowing a better comparisons of 
#' model performances between the different machine learning models. Typically 
#' one would want to keep the full data model but the GBM models can cause the 
#' output object to require large amounts of storage space so optionally one can 
#' choose to not keep the final model when the goal is basically only to assess 
#' model performance for the GBM.  In that case the tuning papramgers for the final
#' tuned modle ae retained faciliting reccalcualtion of the final model, this will
#' also require the original training data.  
#' @param dorf fit and evaluate a random forest (RF) 
#' model.  1 for yes, 0 for no (default).  Also, if dorf is specified by a list, 
#' then RF models will be fit.  The randomForestSRC package is used.  This list can have 
#' three elements.  One is the vector mtryc, and contains values for mtry.  The program 
#' searches over the different values to find a better fir for the final model.  If 
#' not specified mtryc is set to 
#' round( sqrt(dim(xs)[2]) * c(0.67 , 1, 1.5, 2.25, 3.375) ).  The second list element 
#' the vector ntreec.  The first item (ntreec[1]) specifies the number of trees to 
#' fit in evaluating the models specified by the different mtry values.  The second 
#' item (ntreec[2]) specifies the number of trees to fit in the final model.  The
#' default is ntreec = c(25,250).  The third element in the list is the numeric variable keep, with  
#' the value 1 (default) to store the model fit on all data in the output object, or the value 0 
#' to not store the full data model fit.  Typically 
#' one would want to keep the full data model but the RF models can cause the 
#' output object to require large amounts of storage space so optionally one can 
#' choose to not keep the final model when the goal is basically only to assess 
#' model performance for the RF.   Random forests use the out-of-bag (OOB) data elements
#' for assessing model fit and hyperparameter tuning and so cross validation is 
#' not used for tuning.  Still, because of the number of trees in the forest random forest 
#' can take long to run.  
#' @param dorpart fit and do a nested cross validation for an RPART model.  As rpart() does its
#' own approximation for cross validation there is no new functions for cross validation. 
#' @param doann fit and evaluate a cross validation informed Artificial Neural Network 
#' (ANN) model with two hidden levels.  1 for yes, 0 for no (default). By default 
#' the number of folds used when training the ANN model will be the same as the 
#' number of folds used in the outer loop of the nested cross validation.  To override 
#' this, for example to shrtn run time, one may specify a list for the doann argument 
#' where the element $folds_ann_n gives the number of folds used when training the 
#' ANN.  To shorten run we recommend folds_ann_n = folds_n/2 or folds_n/3, and at 
#' least 3.  Then the folds will be formed by collapsing the folds_n folds using 
#' in fitting other models allowing a better comparisons of model performances 
#' between the different machine learning models.  The list can also have 
#' elements $epochs, $epochs2, $myler, $myler2, $eppr, $eppr2, $lenv1, $lenz2, 
#' $actv, $drpot, $wd, wd2, l1, l12, $lscale, $scale, $minloss and $gotoend.  These 
#' arguments are then passed to the ann_tab_cv_best() function, with the 
#' meanings described in the help for that function, with some exception.  When 
#' there are two similar values like $epoch and $epoch2 the first applies to the 
#' ANN models trained without transfer learning and the second to the models 
#' trained with transfer learning from the lasso model.  Elements of this list 
#' unspecified will take default values.  The user may also specify the element 
#' $bestof (a positive integer) to fit bestof models with different random 
#' starting weights and biases while taking the best performing of the different 
#' fits based upon CV as the final model.  The default value for bestof is 1.     
#' @param dostep  fit and do cross validation for stepwise regression fit, 0 or 1, 
#' as discussed in James, Witten, Hastie and Tibshirani, 2nd edition.  
#' @param doaic   fit and do cross validation for AIC fit, 0 or 1.   
#' This is provided primarily as a reference. 
#' @param ensemble This is a vector 8 characters long and specifies a set of ensemble 
#' like model to be fit based upon the predicteds form a relaxed lasso model fit, by 
#' either inlcuding the predicteds as an additional term (feature) in the machine 
#' learning model, or including the predicteds similar to an offset.  For XGBoost, 
#' the offset is specified in the model with the "base_margin" in the XGBoost 
#' call.  For the Artificial Neural Network models fit using the ann_tab_cv_best() function, 
#' one can initialize model weights (parameters) to account for the predicteds in 
#' prediction and either let these weights by modified each epoch or update and maintain 
#' these weights during the fitting process.  For ensemble[1] = 1 a model is fit
#' ignoring these predicteds, ensemble[2]=1 a model is fit including the predicteds 
#' as an additional feature.  For ensemble[3]=1 a model is fit using the predicteds
#' as an offset when running the xgboost model, or a model is fit including the 
#' predicteds with initial weights corresponding to an offset, but then weights are 
#' allowed to be tuned over the epochs.  For i >= 4 ensemble[i] only applies to 
#' the neural network models.  For ensemble[4]=1 a model is fit like for 
#' ensemble[3]=1 but the weights are reassigned to correspond to an offset after 
#' each epoch.  For i in (5,6,7,8) ensemble[i] is similar to ensemble[i-4] except
#' the original predictor (feature) set is replaced by the set of non-zero terms 
#' in the relaxed lasso model fit. If ensemble is specified as 0 or NULL, then ensemble 
#' is assigned c(1,0,0,0, 0,0,0,0).  If ensemble is specified as 1, then ensemble 
#' is assigned c(1,0,0,0, 0,1,0,1).
#' @param method  method for choosing model in stepwise procedure, "loglik" or "concordance".
#' Other procedures use the "loglik". 
#' @param lambda  lambda vector for the lasso fit
#' @param gamma   gamma vector for the relaxed lasso fit, default is c(0,0.25,0.5,0.75,1)
#' @param relax   fit the relaxed lasso model when fitting a lasso model  
#' @param steps_n number of steps done in stepwise regression fitting 
#' @param seed optional, either NULL, or a numerical/integer vector of length 2, for R and torch 
#' random generators, or a list with two two vectors, each of length folds_n+1, for 
#' generation of random folds of the outer cross validation loop, and the remaining 
#' folds_n terms for the random generation of the folds or the bootstrap samples for the 
#' model fits of the inner loops.  This can be used to replicate model fits.  Whether 
#' specified or NULL, the seed is stored 
#' in the output object for future reference.  The stored seed is a list with two
#' vectors seedr for the seeds used in generating the random fold splits, and seedt
#' for generating the random initial weights and biases in the torch neural network 
#' models.  The first element in each of these vectors is for the all data fits and
#' remaining elements for the folds of the inner cross validation.  The integers assigned to 
#' seed should be positive and not more than 2147483647.   
#' @param foldid  a vector of integers to associate each record to a fold.  Should 
#' be integers from 1 and folds_n.  These will only be used in the outer folds. 
#' @param limit limit the small values for lambda after the initial fit.  This 
#' will have minimal impact on the cross validation.  Default is 2 for moderate 
#' limitation, 1 for less limitation, 0 for none.   
#' @param fine  use a finer step in determining lambda.  Of little value unless one 
#' repeats the cross validation many times to more finely tune the hyper paramters.  
#' See the 'glmnet' package documentation  
#' @param ties method for handling ties in Cox model for relaxed model component.  Default 
#' is "efron", optionally "breslow".  For penalized fits "breslow" is 
#' always used as derived form to 'glmnet' package.
#' @param keepdata  0 (default) to delete the input data (xs, start, y_, event) 
#' from the output objects from the random forest fit and the glm() fit for the 
#' stepwise AIC model, 1 to keep.
#' @param keepxbetas  1 (default) to retain in the output object a copy of the 
#' functional outcome variable, i.e. y_ for "gaussian" and "binomial" data, and 
#' the Surv(y_,event) or Surv(start,y_,event) for "cox" data.  This allows 
#' calibration studies of the models, going beyond the linear calibration 
#' information calculated by the function.  The xbetas are calculated both for
#' the model derived using all data as well as for the hold out sets (1/k of the 
#' data each) for the models derived within the cross validation ((k-1)/k of the 
#' data for each fit).
#' @param track   1 (default) to track progress by printing to console elapsed and 
#' split times, 0 to not track
#' @param ... additional arguments that can be passed to glmnet() 
#'  
#' @return - Model fit performance for LASSO, GBM, Random 
#' Forest, RPART, artificial neural network (ANN) or STEPWISE models are 
#' estimated using k-cross validation.  Full data model fits for these models 
#' are also calculated independently (prior to) the performance evaluation, 
#' often using a second layer of k-cross validation. 
#' 
#' @seealso
#'   \code{\link{glmnetr.simdata}} , \code{\link{summary.nested.glmnetr}} , \code{\link{glmnetr.compcv}} , \code{\link{plot.nested.glmnetr}} , 
#'   \code{\link{predict.nested.glmnetr}} , \code{\link{predict_nested_xgb}} , \code{\link{predict_nested_rf}} , \code{\link{predict_ann_tab}}, 
#'   \code{\link{cv.glmnetr}} , \code{\link{xgb.tuned}} , \code{\link{rf_tune}} , \code{\link{ann_tab_cv}} , \code{\link{cv.stepreg}} 
#' 
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @export
#' 
#' @importFrom utils installed.packages  
#' @importFrom stats runif logLik predict cov cor 
#' @importFrom survival Surv coxph coxph.control concordance concordancefit 
#' @importFrom glmnet cv.glmnet 
#' @importFrom Matrix rankMatrix 
#' @importFrom xgboost xgb.DMatrix
## #' @importFrom randomForestSRC predict.rfsrc
#' @importFrom rpart rpart prune 
#' @importFrom torch torch_manual_seed
#' 
#' @examples
#' \donttest{
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$y_ 
#' # for this example we use a small number for folds_n to shorten run time 
#' nested.glmnetr.fit = nested.glmnetr( xs, NULL, y_, NULL, family="gaussian", folds_n=3)
#' plot(nested.glmnetr.fit, type="devrat", ylim=c(0.7,1)) 
#' plot(nested.glmnetr.fit, type="lincal", ylim=c(0.9,1.1)) 
#' plot(nested.glmnetr.fit, type="lasso") 
#' plot(nested.glmnetr.fit, type="coef") 
#' summary(nested.glmnetr.fit) 
#' glmnetr.compcv(nested.glmnetr.fit) 
#' summary(nested.glmnetr.fit, cvfit=TRUE) 
#' }
#' 
nested.glmnetr = function(xs, start=NULL, y_, event=NULL, family="gaussian", do_ncv=1, folds_n=10, stratified=1,  
                          dolasso=1, doxgb=0, dorf=0, doann=0, dorpart=0, dostep=0, doaic=0, 
                          ensemble=0, method="loglik", lambda=NULL, gamma=NULL, relax=TRUE, steps_n=0,
                          seed=NULL, foldid=NULL, limit=1, fine=0, ties="efron", keepdata=0, keepxbetas=1,  
                          track=0, ... ) {
  
  pver = "glmnetr version 0.4-6 (2024-04-21)" 
#  pver = "0.4-6 dev 240421" 
  
  if (is.null(keepdata)) { keepdata = 0 }
  if (is.null(keepxbetas)) { keepxbetas = 0 }
  
  if (is.matrix(start)) { start_name = colnames(start) ; start = as.numeric(start) } else { start_name = "NULL" }
  if (is.matrix(y_   )) { y_name     = colnames(y_)    ; y_    = as.numeric(y_   ) } else { y_name     = "NULL" }
  if (is.matrix(event)) { event_name = colnames(event) ; event = as.numeric(event) } else { event_name = "NULL" }
  
  stp = 0 
  if ( is.null(xs) ) { cat("\n xs cannot be NULL, program will stop\n") ; stp = 1 ; }
  if (sum(is.na(xs))> 0) { cat("\n xs cannot have missing values, program will stop\n") ; stp = 1 ; }
  if ( is.null(y_) ) { cat("\n y_ cannot be missing or NULL, program will stop\n") ; stp = 1 ; }
  if (sum(is.na(y_))> 0) { cat("\n y_ cannot have missing values, program will stop\n") ; stp = 1 ; }
  if (family=="cox") {
    if        (is.null(event))       { cat("\n event cannot be NULL for Cox model, program will stop\n") ; stp = 1 ; 
    } else if (sum(is.na(event))> 0) { cat("\n event cannot have missing values, program will stop\n") ; stp = 1 ; 
    } else if (min(y_) < 0)          { cat("\n survival times cannot be negative for Cox model, program will stop\n") ; stp = 1 ; 
    }
  }
  if (stp == 1) { stop } 
  
  if ( is.logical(y_) & (family == "binomial")) { 
    cat("\n  The y_ variable is of class logical and is converted to numeric\n")
    y_ = y_ * 1 
  }
  
  if (is.null(folds_n)) { folds_n = 10 }
  folds_n = max(folds_n, 2)
  
  ## set up doxgb ==============================================================  
  doxgb_ = 0 
  folds_xgb_n = folds_n 
  if (is.null(doxgb)) { doxgb = 0 }
  if (!is.list(doxgb)) {
    if ( ( doxgb[1] == 0 )) { 
      doxgb_ = 0 
      folds_xgb_n = 0 
    } else if ( doxgb[1] == 1 ) { 
      doxgb_ = 1 
    } else if ( doxgb[1] > 1 ) { 
      doxgb_ = 1 
      folds_xgb_n = max( doxgb[1] , 3) 
    } else {
      doxgb_ = 0 
      folds_xgb_n = 0 
    }
    if ( length(doxgb) >= 2 ) { 
      nrounds = doxgb[2] 
    } else ( nrounds = NULL )
    xgbkeep = 1 
    doxgb = list( nfold=folds_xgb_n, nrounds=nrounds, keep=xgbkeep )
  } else {
    xgbkeep = 1 
    if (!is.null(doxgb$keep)) {
      if (!is.na(doxgb$keep)) {
        xgbkeep = doxgb$keep 
      }
    }
    if (is.null(doxgb$nfold)) { 
      doxgb_ = 1
      doxgb$nfold = folds_n 
      folds_xgb_n = folds_n 
    } else if ( doxgb$nfold < 1) {
      doxgb_ = 0 
      folds_xgb_n = 0
    } else {
      doxgb_ = 1 
      if (doxgb$nfold < 3) { doxgb$nfold = 3 }
      folds_xgb_n = doxgb$nfold 
    }
  }
  
  if ( (doxgb_ == 1) & (family == "cox") & (!is.null(start))) { 
    doxgb_ = 0
    cat( "  The gradient boosting machine fitting routine does not fit Cox model with (start,stop) time data\n",
         " gradient boosting machine models will not be fit\n\n")
  }
  
#  cat(paste('  xgbkeep = ', xgbkeep)) 

  if (is.list(doxgb)) {
    if ( is.null(doxgb$keep) ) { 
      xgbkeep = 1 
      doxgb$keep = xgbkeep 
    } else {
      xgbkeep = doxgb$keep
    } 
    if (is.null( doxgb$nrounds )) { doxgb$nrounds  = 1000 } 
  }

  ## set up dorf ===============================================================  
  ## https://www.randomforestsrc.org/articles/survival.html#ensemble-chf-and-survival-function
  dorf_ = 0 
  
  if (is.null(dorf)) { dorf = 0 }
  
  if (!is.list(dorf)) { 
    if (dorf[1] == 1) { 
      dorf_ = 1 
      dorf = list( keep = 1, nsplitc=8  )
    }
  } else {
    dorf_ = 1 
  }
  
  if (dorf_ == 1) {
    mtryc = dorf$mtryc 
    ntreec = dorf$ntreec 
    nsplitc = dorf$nsplitc 
    if ( is.null(dorf$keep) ) { dorf$keep = 1 }
    rfkeep = dorf$keep 
  }
  
  ## set up doann ==============================================================  
  doann_ = 0 
  if (is.null(doann)) { doann = 0 }
  if ((!is.list(doann)) & ( doann[1] >= 1 )) { doann_ = 1 ;  doann = as.list(doann) }
  if (is.list(doann)) { doann_ = 1 } 
  
  if ( (doann_==1) & (family == "cox") & (!is.null(start))) { 
    cat( "  The neural network fitting routine does not fit Cox model with (start,stop) time data\n",
         " neural network models will not be fit\n\n")
    doann_ = 0
  }

  if ( (dorpart == 1) & (family == "cox") & (!is.null(start))) { 
    dorpart = 0
    cat( "  The recursive partitioning fitting routine does not fit Cox model with (start,stop) time data\n",
         " recursive partitioning machine models will not be fit\n\n")
  }
  
  eppr = NULL 
  eppr2 = NULL 
  if (doann_ == 1) {
    temp_ = list(epochs=200, epochs2=200, mylr=0.005, mylr2=0.001, eppr=-3, eppr2=-3,
                 lenz1=16, lenz2=8, actv=1, drpot=0, lscale=5, scale=1, 
                 wd=0, wd2=0, l1=0, l12=0, folds_ann_n=folds_n, 
                 minloss=0, gotoend=0, bestof=1) 
    if (is.null(doann$epochs )) { epochs = temp_$epochs  ; doann$epochs  = epochs  } else { epochs  = doann$epochs  }
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
    if (is.null(doann$lscale )) { lscale = temp_$lscale  ; doann$lscale  = lscale  } else { lscale  = doann$lscale  }
    if (is.null(doann$scale  )) { scale  = temp_$scale   ; doann$scale   = scale   } else { scale   = doann$scale   }
    if (is.null(doann$folds_ann_n)) { folds_ann_n = temp_$folds_ann_n ; doann$folds_ann_n = folds_ann_n  } else { folds_ann_n  = doann$folds_ann_n  }
    if (is.null(doann$minloss)) { minloss= temp_$minloss ; doann$minloss = minloss } else { minloss = doann$minloss }
    if (is.null(doann$gotoend)) { gotoend= temp_$gotoend ; doann$gotoend = gotoend } else { gotoend = doann$gotoend }
    if (is.null(doann$bestof )) { bestof = temp_$bestof  ; doann$bestof  = bestof  } else { bestof  = doann$bestof  }
    if (folds_ann_n == 0) { doann_ = 0 } 
    if (folds_ann_n == 1) { 
      folds_ann_n = 3
      doann$folds_ann_n = folds_ann_n
      cat(paste(" --- folds_ann_n cannot be 1 so is set to 3 ---"))
    }
    if ((track >= 1) & (doann_ == 1)) {
      cat(paste0(" epochs=", epochs, " epochs2=", epochs2, " mylr=", mylr, " mylr2=", mylr2, " eppr=", eppr,
                 " eppr2=", eppr2, " lenz1=", lenz1, " lenz2=", lenz2, " actv=", actv, " drpot=", drpot, "\n" )) 
      cat(paste0(" wd=", wd, " wd2=", wd2, " l1=", l1, " l12=", l12, " lscale=", lscale, " scale=", scale, " folds_ann_n=", folds_ann_n,  
                 " minloss=", minloss, " gotoend=", gotoend, " bestof=", bestof, "\n" )) 
    }
  } else { folds_ann_n = folds_n }
  
  if (!("torch" %in% installed.packages()[,1])) {
    doann_ = 0 
    cat("\n  ***** torch is NOT installed and so neural networks models will not be fit *****\n")
  }
  
  if (track <= 0) { if (is.null(eppr)) { eppr = -3 } ;  if (is.null(eppr)) { eppr2 = -3 } }
  if (track <= 0) { eppr = min(eppr, -2) ;  eppr2 = min(eppr2,-2) }
  if (track >= 1) { cat(paste0("\n", " ##############################################################################################" , "\n")) }
  time_start = Sys.time()
  time_last = NULL 
  if (track >= 1) { time_start = diff_time() }

  ## general input check #######################################################

  if (ties != "breslow") { ties="efron" }
  
  if ( is.null(ensemble) ) { ensemble = c(1, rep(0,7)) } 
  if (sum(ensemble)   ==0) { ensemble = c(1, rep(0,7)) 
  } else if (length(ensemble)==1) { ensemble = c(1,0,0,0, 0,1,0,1) }
  temp_ = ensemble
  lng = min(length(temp_), 9)
  ensemble = rep(0,9) ; 
  ensemble[1:lng] = temp_ 
  
  if (((doxgb_==1) | (doann_==1) | (dorf_==1) | (dorpart==1)) & (sum(ensemble[c(2:8)]) >= 1)) { dolasso = 1 } 
  
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
  
  #####  SEEDS  ################################################################
  
  seed = glmnetr_seed(seed, folds_n, folds_ann_n) 

  seedr_ = seed$seedr[1]
  seedt_ = seed$seedt[1]

  if (is.null(foldid)) { 
    set.seed(seedr_) 
    foldid = get.foldid(y_, event, family, folds_n, stratified) 
  }
  table (foldid)
  
  if (doxgb_ == 1) { 
    if (!is.null(doxgb$folds)==1) {
      foldid_xgb = doxgb$folds  
      foldid_xgb_v = rep(0,nobs)
      for (i1 in c(1:folds_xgb_n)) { foldid_xgb_v[foldid_xgb[[i1]]] = i1 }
    } else if (folds_xgb_n == folds_n) {
      foldid_xgb_v = foldid 
      foldid_xgb = list() 
      for (i1 in c(1:folds_xgb_n)) { foldid_xgb[[i1]] = c(1:nobs)[foldid==i1] } 
    } else {
      if ( (folds_n / folds_xgb_n) == 2 ) {
        foldid_xgb_v = (foldid + 1) %/% 2 
      } else if ( (folds_n / folds_xgb_n) == 3 ) {
        foldid_xgb_v = (foldid + 1) %/% 3 
      } else {
        set.seed(seedr_) 
        foldid_xgb_v = get.foldid(y_, event, family, folds_xgb_n, stratified) 
      }
    } 
    table( foldid_xgb_v, foldid)
    foldid_xgb = list()
    for (i1 in c(1:folds_xgb_n)) { foldid_xgb[[i1]] = c(1:nobs)[foldid_xgb_v == i1] } 
    doxgb$folds = foldid_xgb 
  }
  
#  doann$foldid_ann = NULL 
  if (doann_ == 1) { 
    torch_manual_seed( seedt_ ) 
    if (!is.null(doann$foldid_ann)==1) {
      foldid_ann = doann$foldid_ann 
    } else if (folds_ann_n == folds_n) {
      foldid_ann = foldid
    } else {
      if ( (folds_n / folds_ann_n) == 2 ) {
        foldid_ann = (foldid + 1) %/% 2 
      } else if ( (folds_n / folds_ann_n) == 3 ) {
        foldid_ann = (foldid + 1) %/% 3 
      } else {
        set.seed(seedr_) 
        foldid_ann = get.foldid(y_, event, family, folds_ann_n, stratified)  
      }
    }
    doann$foldid_ann = foldid_ann 
    table( foldid_ann, foldid)
  } 

  ## no model statistics #######################################################
  null.m2LogLik.cv = rep(0,folds_n) 
  sat.m2LogLik.cv = rep(0,folds_n) 
  n.cv = rep(0,folds_n) 
  
  ## log likelihoods & concordances from LASSO and relaxed lasso by Cross Validation 
  zero_nby7 = matrix ( data=rep(0,(7*folds_n)), nrow = folds_n, ncol = 7 ) 
  ones_nby7 = matrix ( data=rep(1,(7*folds_n)), nrow = folds_n, ncol = 7 ) 
  lasso.devian.cv = ones_nby7
  lasso.cal.devian.cv = ones_nby7
  lasso.agree.cv  = zero_nby7
  lasso.intcal.cv = zero_nby7
  lasso.lincal.cv = zero_nby7
  ## number of non-zero coefficients in LASSO models 
  lasso.nzero.cv =  zero_nby7
  lassogammacv = matrix ( data=rep(0,(2*folds_n)), nrow = folds_n, ncol = 2)
  nms = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0", "ridge" )
  nms = c("1se", "min", "1seR", "minR", "1seR.G0", "minR.G0", "ridge" )
  colnames(lasso.devian.cv) = nms 
  colnames(lasso.cal.devian.cv) = nms 
  colnames(lasso.agree.cv ) = nms
  colnames(lasso.intcal.cv) = nms
  colnames(lasso.lincal.cv) = nms
  colnames(lasso.nzero.cv)  = nms
  colnames(lassogammacv)    = c("1se", "min")
  
  ## log likelihoods & concordances from XGBoost by Cross Validation 
  zero_nby6 = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 ) 
  ones_nby6 = matrix ( data=rep(1,(6*folds_n)), nrow = folds_n, ncol = 6 ) 
  xgb.devian.cv = ones_nby6
  xgb.cal.devian.cv = ones_nby6 
  xgb.intcal.cv = zero_nby6
  xgb.lincal.cv = zero_nby6
  xgb.agree.cv  = zero_nby6
  xgb.nzero.cv  = zero_nby6
  xgbnms = c("Simple", "Feature", "Offset", "Tuned", "Feature", "Offset") 
  colnames(xgb.devian.cv) = xgbnms
  colnames(xgb.cal.devian.cv) = xgbnms
  colnames(xgb.intcal.cv) = xgbnms
  colnames(xgb.lincal.cv) = xgbnms
  colnames(xgb.agree.cv ) = xgbnms
  colnames(xgb.nzero.cv ) = xgbnms 
  
  ## log likelihoods & concordances from Random Forest by Cross Validation   
  zero_nby3 = matrix ( data=rep(0,(3*folds_n)), nrow = folds_n, ncol = 3 ) 
  ones_nby3 = matrix ( data=rep(1,(3*folds_n)), nrow = folds_n, ncol = 3 ) 
  rf.mtry.cv     = zero_nby3
  rf.devian.cv   = ones_nby3
  rf.cal.devian.cv = ones_nby3
  rf.agree.cv    = zero_nby3
  rf.intcal.cv   = zero_nby3
  rf.lincal.cv   = zero_nby3
  nms = c("None", "Feature", "Offset")
  colnames(rf.mtry.cv)  = nms 
  colnames(rf.devian.cv) = nms 
  colnames(rf.cal.devian.cv) = nms 
  colnames(rf.agree.cv)  = nms 
  colnames(rf.intcal.cv) = nms 
  colnames(rf.lincal.cv) = nms 
  
  ## log likelihoods & concordances from RPART by Cross Validation   
  rpnz = function(obj) { 
    frame = obj$frame ;  leaves = frame$var == "<leaf>" ;  used <- unique(frame$var[!leaves]) ; 
    return( length(used) )
  }
  ones_nby9 = matrix ( data=rep(1,(9*folds_n)), nrow = folds_n, ncol = 9 )
  zero_nby9 = matrix ( data=rep(0,(9*folds_n)), nrow = folds_n, ncol = 9 )
  rpart.nzero.cv  = zero_nby9
  rpart.devian.cv = ones_nby9
  rpart.cal.devian.cv = ones_nby9
  rpart.agree.cv  = zero_nby9
  rpart.intcal.cv = zero_nby9
  rpart.lincal.cv = zero_nby9
  
  nms = c("cp=0.00", "cp=0.01", "cp=0.02", "F cp=0.00", "F cp=0.01", "F cp=0.02", 
          "O cp=0.00", "O cp=0.01", "O cp=0.02" )
  colnames(rpart.nzero.cv)  = nms 
  colnames(rpart.devian.cv) = nms 
  colnames(rpart.cal.devian.cv) = nms 
  colnames(rpart.agree.cv)  = nms 
  colnames(rpart.intcal.cv) = nms 
  colnames(rpart.lincal.cv) = nms 
  
  ones_nby8 = matrix ( data=rep(1,(8*folds_n)), nrow = folds_n, ncol = 8 )
  zero_nby8 = matrix ( data=rep(0,(8*folds_n)), nrow = folds_n, ncol = 8 )
  ## log likelihoods & concordances from Neural Networks by Cross Validation   
  ann.devian.cv   = ones_nby8
  ann.cal.devian.cv = ones_nby8
  ann.agree.cv    = ones_nby8
  ann.intcal.cv   = ones_nby8
  ann.lincal.cv   = ones_nby8
  ann.nzero.cv    = ones_nby8
  nms = c("Uninformed", "lasso feat", "lasso w's", "lasso update", "lasso terms", "l/lasso feat", "l/lasso w's", "l/lasso update") 
  colnames(ann.devian.cv) = nms 
  colnames(ann.cal.devian.cv) = nms 
  colnames(ann.intcal.cv) = nms 
  colnames(ann.lincal.cv) = nms 
  colnames(ann.agree.cv)  = nms 
  
  ## log likelihoods & conncordances for STEPWISE by Cross Validation 
  ## and coxph models based upon LASSO terms by Cross Validation
  step.devian.cv = ones_nby3
  step.cal.devian.cv = ones_nby3
  step.intcal.cv = zero_nby3
  step.lincal.cv = zero_nby3
  step.agree.cv  = zero_nby3
  step.nzero.cv  = zero_nby3
  step.p.cv      = matrix ( data=rep(0,(folds_n  )), nrow = folds_n, ncol = 1 )
  colnames(step.devian.cv) = c("df", "p", "AIC")
  colnames(step.cal.devian.cv) = c("df", "p", "AIC")
  colnames(step.intcal.cv) = c("df", "p", "AIC")
  colnames(step.lincal.cv) = c("df", "p", "AIC")
  colnames(step.agree.cv ) = c("df", "p", "AIC") 
  colnames( step.nzero.cv  )  = c("df", "p", "AIC")
  colnames( step.p.cv   )  = c("p stepwise")
  
  #########################################################################################################
  #########################################################################################################
  
  perf_gau = function(yy, xbhat) {
    fit0 = glm( yy ~ xbhat , family=family)
    returnvec = c(  dev1 = sum((yy - xbhat)^2) / length(yy) ,
                    dev2 = fit0$deviance / length(yy) ,
                    agree = ifelse( var(xbhat)>0 , cor(x=yy, y=xbhat) , 0 ) ,
                    intcal = fit0$coefficients[1] , 
                    lincal = fit0$coefficients[2] )
    return( returnvec )
  }
  
  perf_bin = function(yy, xbhat) {
    fit0 = glm( yy ~ xbhat , family=family) 
    p_ = 1/(1+exp(-xbhat)) 
    retvec = c( dev1 = -2*( t(log(p_))%*%yy + t(log(1-p_))%*%(1-yy) ) / length(yy) , 
                dev2 = fit0$deviance / length(yy) ,
                agree  = concordance(fit0)[[1]] ,          
                intcal = fit0$coefficients[1] , 
                lincal = fit0$coefficients[2] ) 
    return( retvec )
  }
  
  perf_cox = function(SURV, pred) {
    fit0 = coxph( SURV ~ pred, init=c(1)) 
    retvec = c(dev1 = -2*fit0$loglik[1] / fit0$nevent ,
               dev2 = -2*fit0$loglik[2] / fit0$nevent , 
               agree = fit0$concordance[6] ,
               intcal = 0 , 
               lincal = fit0$coefficients ) 
    return( retvec )
  }
  
  perf_gen = function(yy, pred, family) {
    if        (family == "cox")      { returnvec = perf_cox(yy, pred) 
    } else if (family == "binomial") { returnvec = perf_bin(yy, pred) 
    } else if (family == "gaussian") { returnvec = perf_gau(yy, pred) } 
    return(returnvec)
  }
  
  if (family=="cox") {
    if (is.null(start)) { 
      y__ = Surv(y_, event)
    } else { 
      y__ = Surv(start, y_, event) 
    }
  } else {
    y__ = y_
  }
  
  #########################################################################################################
  
  if (family == "cox") { xbnull = 0 
  } else if (family == "binomial") {
    p = mean(y_) 
    if (p > (1-tol_)) { p = 1/tol_
    } else if (p < tol_) { p = tol_ }
    xbnull = log(p/(1-p))
  } else {
    xbnull = mean(y_) 
  }
  
  xbetas     = matrix(rep(xbnull,1*nobs), nrow=nobs, ncol=1)
  xbetas0    = matrix(rep(xbnull,7*nobs), nrow=nobs, ncol=7)
  xbetas.cv  = matrix(rep(xbnull,1*nobs), nrow=nobs, ncol=1)
  xbetas0.cv.lasso = matrix(rep(xbnull,7*nobs), nrow=nobs, ncol=7)
  
  #########################################################################################################
  ##### LASSO fit #########################################################################################
  if (dolasso == 1) {
    
    if (track >= 1) { cat(paste0("\n", " ########## Initial (CV) lasso fit of all data ########################################" , "\n")) }
    if (relax==1) { 
      if (is.null(gamma)) { gamma = c(0, 0.25, 0.5, 0.75, 1) }
      
      cv_glmnet_fit_f = cv.glmnetr( xs, start, y_, event, family=family, 
                                    lambda=lambda, gamma=gamma, folds_n=folds_n, foldid=foldid, limit=limit, fine=fine, track=track, ties=ties, ... )  
    } else {
      if (family=="cox") {
        if ( is.null(start) ) { cv_glmnet_fit_f = cv.glmnet( xs, Surv(y_, event)       , family="cox", lambda=lambda, foldid=foldid, relax=FALSE, ... )  
        } else                { cv_glmnet_fit_f = cv.glmnet( xs, Surv(start, y_, event), family="cox", lambda=lambda, foldid=foldid, relax=FALSE, ... )  
        } 
      } else {
        cv_glmnet_fit_f = cv.glmnet( xs, y_, family=family, lambda=lambda, foldid=foldid, relax=FALSE, ... )  
      }
    }
    
    if (folds_n >= 3) {
      if ((family == "cox") & (is.null(start))) {
        cv_ridge_fit_f = cv.glmnet( xs, Surv(y_, event), family=family, alpha=0, nfolds=folds_n, foldid=foldid, ... )  
      } else if ((family == "cox") & (!is.null(start))) {
        cv_ridge_fit_f = cv.glmnet( xs, Surv(start, y_, event), family=family, alpha=0, nfolds=folds_n, foldid=foldid, ... )  
      } else {
        cv_ridge_fit_f = cv.glmnet( xs, y_, family=family, alpha=0, nfolds=folds_n, foldid=foldid, ... )  
      }
    }
    
#    cat(" table(foldid)" )
#    cat(table(foldid))
    
    #------------------------------------------
    
    if (family == "cox") { 
      if ( is.null(start) ) { SURV = Surv(y_, event) 
      } else {                SURV = Surv(start, y_, event) }
    }

    pred1se   = predict(cv_glmnet_fit_f, xs, lam=cv_glmnet_fit_f$lambda.1se, gam=1)
    predmin   = predict(cv_glmnet_fit_f, xs, lam=cv_glmnet_fit_f$lambda.min, gam=1)
    pred1seR  = predict(cv_glmnet_fit_f, xs, lam="lambda.1se" , gam="gamma.1se" )
    predminR  = predict(cv_glmnet_fit_f, xs, lam="lambda.min" , gam="gamma.min" ) 
#    predminR  = predict(cv_glmnet_fit_f, xs) #  lam="lambda.min" , gam="gamma.min" ) 
    pred1seR0 = predict(cv_glmnet_fit_f, xs, lam=cv_glmnet_fit_f$relaxed$lambda.1se.g0, gam=0)
    predminR0 = predict(cv_glmnet_fit_f, xs, lam=cv_glmnet_fit_f$relaxed$lambda.min.g0, gam=0)
#    predminR0 = predict(cv_glmnet_fit_f, xs, gam=0)
    if (folds_n >= 3) { predridge = predict(cv_ridge_fit_f , xs, s="lambda.min")
    } else {  predridge = rep(0,dim(xs)[1]) }
    
#      cat( cor(cbind(pred1se, predmin, pred1seR, predminR, pred1seR0, predminR0)) ) 
#      cat( cor(cbind(predmin, predminR, predminR0, predminRp, y_)) ) 
    
    if (family == "cox") {    
      if ( is.null(start) ) { yy = Surv(y_, event) 
      } else {                yy = Surv(start, y_, event) }
    } else { yy = y_ }
    
    perfm1 = perf_gen( yy , pred1se   , family )
    perfm2 = perf_gen( yy , predmin   , family )
    perfm3 = perf_gen( yy , pred1seR  , family )
    perfm4 = perf_gen( yy , predminR  , family )
    perfm5 = perf_gen( yy , pred1seR0 , family )
    perfm6 = perf_gen( yy , predminR0 , family )
    perfm7 = perf_gen( yy , predridge , family )
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

    lasso.nzero = rep(0,7) 
    lasso.nzero[1] = cv_glmnet_fit_f$nzero [ cv_glmnet_fit_f$index[2] ]
    lasso.nzero[2] = cv_glmnet_fit_f$nzero [ cv_glmnet_fit_f$index[1] ]     
    if (relax) {
      lasso.nzero[3] = cv_glmnet_fit_f$relaxed$nzero.1se 
      lasso.nzero[4] = cv_glmnet_fit_f$relaxed$nzero.min 
#      lasso.nzero[4] = length(predict(cv_glmnet_fit_f)$beta)
      lasso.nzero[5] = cv_glmnet_fit_f$nzero [ cv_glmnet_fit_f$relaxed$index.g0[2] ]
      lasso.nzero[6] = cv_glmnet_fit_f$nzero [ cv_glmnet_fit_f$relaxed$index.g0[1] ]
#      lasso.nzero[6] = length(predict(cv_glmnet_fit_f, gam=0)$beta)
    }
    if (folds_n >= 3) { lasso.nzero[7] = cv_ridge_fit_f$nzero[ cv_ridge_fit_f$index ][1]
    } else { lasso.nzero[7] = 0 }
    names(lasso.nzero) = nms 
    
    ## calibrate lasso ##
    if (sum(ensemble[c(2:4,6:8)]) >= 1) {
      ofst = lasso.intcal.naive[4] + lasso.lincal.naive[4] * predminR 
      if (family=="cox") { ofst = ofst - mean(ofst) }
    }
    
    xbetas0 = cbind(pred1se, predmin, pred1seR, predminR, pred1seR0, predminR0, predridge)
    colnames(xbetas0) = c("Lasso.1se", "lasso.min", "lassoR.1se", "lassoR.min", "lassoR0.1se", "lassoR0.min", "ridge" )
    colnames(xbetas0.cv.lasso) = c("Lasso.1se", "lasso.min", "lassoR.1se", "lassoR.min", "lassoR0.1se", "lassoR0.min", "ridge" )
    xbetas = cbind(xbetas, xbetas0)
#    cor(xbetas)
    
    if (track >= 1) { cat(paste0(" length(lambda) = " , length(cv_glmnet_fit_f$lambda), "\n" )) } 
    
    if (track >= 1) { time_last = diff_time(time_start, time_last)  }
  }
  
  ###############################################################################################################
  ##### XGBoost fit #############################################################################################
  if (doxgb_==1) { 
    if (track >= 1) { cat(paste0("\n", " ########## Initial XGBoost fit on all data ###########################################" , "\n")) }
    
    if (family=="cox") { Surv.xgb = ifelse( event == 1, y_, -y_) }
    
    if (sum(ensemble[c(1,5)]) >= 1) {
      if (family=="cox") { full.xgb.dat <- xgb.DMatrix(data = xs, label = Surv.xgb)
      } else {             full.xgb.dat <- xgb.DMatrix(data = xs, label = y_)      }
    }
    
    if (sum(ensemble[c(2,6)]) >= 1) {
      if (family=="cox") { full.xgb.datF <- xgb.DMatrix(data = cbind(xs,ofst), label = Surv.xgb)
      } else {             full.xgb.datF <- xgb.DMatrix(data = cbind(xs,ofst), label = y_)      }
    }
    
    if (sum(ensemble[c(3,4,7,8)]) >= 1) {
      if (family=="cox") { full.xgb.datO <- xgb.DMatrix(data = xs, label = Surv.xgb, base_margin=ofst)
      } else {             full.xgb.datO <- xgb.DMatrix(data = xs, label = y_, base_margin=ofst)      }
    }
    
    xgb_perform = function(xgb_model, xs_, y=y__, family=family, tol=tol_) {
      xbhat = xgb_xbhat(xgb_model, xs_, family, tol) 
      if (family == "cox") {
        returnvec = perf_cox( y , xbhat) 
      } else if (family == "binomial") {
        returnvec = perf_bin(y , xbhat) 
      } else if (family == "gaussian") {
        returnvec = perf_gau(y , xbhat) 
      }
      return( returnvec ) 
    }
    
    if (family %in% c("cox","binomial", "gaussian")) { 
      
      if        (family == "cox"     ) { objective = "survival:cox" 
      } else if (family == "binomial") { objective = "binary:logistic" 
      } else                           { objective = "reg:squarederror" }

      xbetas0 = matrix(rep(xbnull,6*nobs), nrow=nobs, ncol=6)
      xbetas0.cv.xgb = matrix(rep(xbnull,6*nobs), nrow=nobs, ncol=6)
      xbnames = c("xgb.simple", "xgb.sim.feat", "xgb.sim.offs", "xgb.tuned", "xgb.tun.feat", "xgb.tun.offs")
      colnames(xbetas0) = xbnames
      colnames(xbetas0.cv.xgb) = xbnames
      
      ## SIMPLE XGB model full data ################################################
      if (sum(ensemble[c(1,5)]) >= 1) {
        if (track >= 2) { cat("  XGB Simple  ") }
        xgb.simple.fit = xgb.simple(full.xgb.dat, objective = objective, seed=seedr_, folds=foldid_xgb, doxgb=doxgb, track=track )
        xgbperf1 = xgb_perform(xgb_model=xgb.simple.fit, xs_=full.xgb.dat, y=y__, family=family )
        doxgb_simple = xgb.simple.fit$doxgb
        xbetas0[,1] = xgb_xbhat(xgb.simple.fit, full.xgb.dat, family)
      } else { xgbperf1 = c(1,1,0,0,0) }
      if (sum(ensemble[c(2,6)]) >= 1) {
        if (track >= 2) { cat("\n  XGB Simple Factor  ") } 
        xgb.simple.fitF = xgb.simple(full.xgb.datF, objective = objective, seed=seedr_, folds=foldid_xgb, doxgb=doxgb, track=track )
        xgbperf2 = xgb_perform(xgb.simple.fitF, full.xgb.datF, y=y__, family=family )
        doxgb_simpleF = xgb.simple.fitF$doxgb
        xbetas0[,2] = xgb_xbhat(xgb.simple.fitF, full.xgb.datF, family)
      } else { xgbperf2 = c(1,1,0,0,0) }
      if (sum(ensemble[c(3,4,7,8)]) >= 1) {
        if (track >= 2) { cat("\n  XGB Simple Offset  ") }
        xgb.simple.fitO = xgb.simple(full.xgb.datO, objective = objective, seed=seedr_, folds=foldid_xgb, doxgb=doxgb, track=track )
        xgbperf3 = xgb_perform(xgb.simple.fitO, full.xgb.datO, y=y__, family=family )
        doxgb_simpleO = xgb.simple.fitO$doxgb
        xbetas0[,3] = xgb_xbhat(xgb.simple.fitO, full.xgb.datO, family)
      }  else { xgbperf3 = c(1,1,0,0,0) }
      if (track >= 2) { cat("\n") } 
      ## TUNED XGB fit full data ###################################################
      if (sum(ensemble[c(1,5)]) >= 1) {
        if (track >= 2) { cat("\n  XGB Tuned  ") }
        xgb.tuned.fit = xgb.tuned(full.xgb.dat, objective = objective, seed=seedr_, folds=foldid_xgb, doxgb=doxgb, track=track )
        xgbperf4   = xgb_perform(xgb.tuned.fit, full.xgb.dat, y=y__, family=family )
        xbetas0[,4] = xgb_xbhat(xgb.tuned.fit, full.xgb.dat, family)
      }  else { xgbperf4 = c(1,1,0,0,0) }
      if (sum(ensemble[c(2,6)]) >= 1) {
        if (track >= 2) { cat("\n  XGB Tuned Factor  ") }
        xgb.tuned.fitF = xgb.tuned(full.xgb.datF, objective = objective, seed=seedr_, folds=foldid_xgb, doxgb=doxgb, track=track )
        xgbperf5 = xgb_perform(xgb.tuned.fitF, full.xgb.datF, y=y__, family=family )
        xbetas0[,5] = xgb_xbhat(xgb.tuned.fitF, full.xgb.datF, family)
      }  else { xgbperf5 = c(1,1,0,0,0) }
      if (sum(ensemble[c(3,4,7,8)]) >= 1) {
        if (track >= 2) { cat("\n  XGB Tuned Offset  ") }
        xgb.tuned.fitO = xgb.tuned(full.xgb.datO, objective = objective, seed=seedr_, folds=foldid_xgb, doxgb=doxgb, track=track )
        xgbperf6 = xgb_perform(xgb.tuned.fitO, full.xgb.datO, y=y__, family=family )
        xbetas0[,6] = xgb_xbhat(xgb.tuned.fitO, full.xgb.datO, family)
      }  else { xgbperf6 = c(1,1,0,0,0) }
      if (track >= 2) { cat("\n") } 
      xgb.devian.naive     = c( xgbperf1[1] , xgbperf2[1] , xgbperf3[1] , xgbperf4[1] , xgbperf5[1] , xgbperf6[1] )
      xgb.cal.devian.naive = c( xgbperf1[2] , xgbperf2[2] , xgbperf3[2] , xgbperf4[2] , xgbperf5[2] , xgbperf6[2] )
      xgb.agree.naive      = c( xgbperf1[3] , xgbperf2[3] , xgbperf3[3] , xgbperf4[3] , xgbperf5[3] , xgbperf6[3] )
      xgb.intcal.naive     = c( xgbperf1[4] , xgbperf2[4] , xgbperf3[4] , xgbperf4[4] , xgbperf5[4] , xgbperf6[4] )
      xgb.lincal.naive     = c( xgbperf1[5] , xgbperf2[5] , xgbperf3[5] , xgbperf4[5] , xgbperf5[5] , xgbperf6[5] )
      xgb.nzero = rep(xs_ncol, 6)

      nms = c("Simple", "Feature", "Offset", "Tuned", "Feature", "Offset") 
      names(xgb.devian.naive) = nms
      names(xgb.cal.devian.naive) = nms
      names(xgb.intcal.naive) = nms
      names(xgb.lincal.naive) = nms
      names(xgb.agree.naive)  = nms
      names(xgb.nzero)  = nms
      
      ##########################################################################
    }
    
    xbetas = cbind(xbetas, xbetas0)

    if (track >= 1) { time_last = diff_time(time_start, time_last)  }
  }
  
  ###############################################################################################################
  ##### RandomForest SRC ########################################################################################
  if (dorf_==1) { 
    if (track >= 1) { cat(paste0("\n", " ########## RandomForest fit on all data ##############################################" , "\n")) }
    
    rf_tuned = NULL ; 
    
    rf.devian.naive = c(0,0,0)
    rf.cal.devian.naive = c(0,0,0)
    rf.agree.naive  = c(0,0,0) 
    rf.intcal.naive = c(0,0,0)
    rf.lincal.naive = c(0,0,0)
    rf.mtry.naive   = c(0,0,0)

    nms = c("Standard", "Feature", "Offset" ) 
    names(rf.devian.naive)  = nms 
    names(rf.cal.devian.naive) = nms 
    names(rf.agree.naive )  = nms
    names(rf.intcal.naive)  = nms
    names(rf.lincal.naive)  = nms
    names(rf.mtry.naive )   = nms
    
    rf_perform = function(rf_model, dframe, ofst=NULL, y__, family, tol=tol_) {
      xbhat = rf_xbhat(rf_model, dframe, ofst, family, tol) 
      if (family == "cox") {
        returnvec = perf_cox(y__, xbhat)
      } else if (family == "binomial") {
        returnvec = perf_bin(y__, xbhat) 
      } else if (family == "gaussian") {
        returnvec = perf_gau(y__, xbhat ) 
      }
      returnvec = c(returnvec, rf_model$mtry)
    }
    
    xbetas0 = matrix(rep(xbnull,3*nobs), nrow=nobs, ncol=3)
    xbetas0.cv.rf = matrix(rep(xbnull,3*nobs), nrow=nobs, ncol=3)
    colnames(xbetas0) = c("rf", "rf.feat", "rf.offs")
    colnames(xbetas0.cv.rf) = c("rf", "rf.feat", "rf.offs")
    
    if (sum(ensemble[c(1,5)]) >= 1) {
      if (track >= 2) { cat( "     No lasso info\n") } 
      rf_tuned = rf_tune(xs=xs, y_=y_, event=event, family=family, mtryc=mtryc, ntreec = ntreec, nsplitc=nsplitc, seed=seedr_) 
      rf_perf1  = rf_perform(rf_model=rf_tuned$rf_tuned, dframe=as.data.frame(xs), ofst=NULL, y__=y__, family=family, tol=tol_)
      xbetas0[,1] = rf_xbhat(rf_model=rf_tuned$rf_tuned, dframe=as.data.frame(xs), ofst=NULL, family=family, tol=tol_)
      if (track >= 2) { time_lasti = diff_time(time_start, time_last) } 
    }  else { rf_perf1 = rep(1,1,0,0,0,) }
    
    if (sum(ensemble[c(2,6)]) >= 1) { 
      if (track >= 2) { cat( "     Relaxed lasso predicteds as feature\n") }
      rf_tunedF = rf_tune(xs=cbind(xs, ofst), y_=y_, event=event, family=family, mtryc=mtryc, ntreec = ntreec, nsplitc=nsplitc, seed=seedr_)
      rf_perf2 = rf_perform(rf_tunedF$rf_tuned , as.data.frame(cbind(xs,ofst)), NULL, y__, family )
      xbetas0[,2] = rf_xbhat(rf_tunedF$rf_tuned , as.data.frame(cbind(xs,ofst)), NULL, family )
      if (track >= 2) { time_lasti = diff_time(time_start, time_lasti) } 
    } else { rf_perf2 = rep(1,1,0,0,0,) }
    
    if ( (sum(ensemble[c(3,4,7,8)]) >= 1) & (family == "gaussian") ) { 
      if (track >= 2) { cat( "     Relaxed lasso predicteds as offset\n") }
      yo = y_ - ofst                                                            ## 231228 
      rf_tunedO = rf_tune(xs=xs, y_=yo, event=event, family=family, mtryc=mtryc, ntreec = ntreec, nsplitc=nsplitc, seed=seedr_)
      rf_perf3   = rf_perform(rf_tunedO$rf_tuned , as.data.frame(xs), ofst, y__, family )
      xbetas0[,3] = rf_xbhat(rf_tunedO$rf_tuned , as.data.frame(xs), ofst, family )
      if (track >= 2) { time_lasti = diff_time(time_start, time_lasti) } 
    } else { rf_perf3 = rep(1,1,0,0,0,) }
    if (track >= 2) { cat("\n") }
    
    rf.devian.naive     = c( rf_perf1[1] , rf_perf2[1] , rf_perf3[1] ) 
    rf.cal.devian.naive = c( rf_perf1[2] , rf_perf2[2] , rf_perf3[2] ) 
    rf.agree.naive      = c( rf_perf1[3] , rf_perf2[3] , rf_perf3[3] ) 
    rf.intcal.naive     = c( rf_perf1[4] , rf_perf2[4] , rf_perf3[4] ) 
    rf.lincal.naive     = c( rf_perf1[5] , rf_perf2[5] , rf_perf3[5] ) 
    rf.mtry.naive       = c( rf_perf1[6] , rf_perf2[6] , rf_perf3[6] )
    
    if (sum(ensemble[c(1,5)]==1) >= 1) { 
      object1 = rf_tuned 
      temp1 = list( err.ratev = object1$err.ratev, mtryc = object1$mtryc , ntreec = object1$ntreec , nsplitc=nsplitc, 
                   mtry = object1$rf_tuned$mtry, nsplit=object1$rf_tuned$nsplit,  keep = rfkeep, seed=seedr_ ) 
      dorf_base = temp1  
    }
    if (sum(ensemble[c(2,6)]==1) >= 1) {
      object1 = rf_tunedF 
      temp1 = list( err.ratev = object1$err.ratev, mtryc = object1$mtryc , ntreec = object1$ntreec , nsplitc=nsplitc, 
                   mtry = object1$rf_tuned$mtry, nsplit=object1$rf_tuned$nsplit,  keep = rfkeep, seed=seedr_ ) 
      dorf_feat = temp1  
    }
    if ( (sum(ensemble[c(3,4,7,8)]) >= 1) & (family == "gaussian") ) {
      object1 = rf_tunedO 
      temp1 = list( err.ratev = object1$err.ratev, mtryc = object1$mtryc , ntreec = object1$ntreec , nsplitc=nsplitc, 
                   mtry = object1$rf_tuned$mtry, nsplit=object1$rf_tuned$nsplit,  keep = rfkeep, seed=seedr_ ) 
      dorf_offs = temp1 
    }
    
    xbetas = cbind(xbetas, xbetas0)

    if (track >= 1) { time_last = diff_time(time_start, time_last) } 
  }
  
  ###############################################################################################################
  ##### RPART fit ###############################################################################################
  if (dorpart==1) { 
    if (track >= 1) { cat(paste0("\n", " ########## Initial RPART fit on all data #############################################" , "\n")) }
    
    rpart.fit.00 = NULL ; rpart.fitO.00 = NULL ; rpart.fitF.00 = NULL ;

    df.xs = as.data.frame(xs)
    
    if (sum(ensemble[c(2:4,6:8)]) >= 1) { 
      xsO    = cbind(xs, ofst=ofst) 
      df.xsO = as.data.frame(xsO)
    }
    
    rpart_perform = function(rpart_model, dframe, ofst=NULL, y__, family, tol=1e-5) {
      xbhat = rpart_xbhat(rpart_model, dframe, ofst, family, tol)
      if (family == "cox") {
        returnvec = perf_cox(y__, xbhat)
      } else if (family == "binomial") {
        returnvec = perf_bin(y__, xbhat)    
      }else if (family == "gaussian") {
        returnvec = perf_gau(y__, xbhat) 
      }
      returnvec = c(returnvec, rpnz(rpart_model))
      return( returnvec )
    }
    
    if (family == "cox") {
      rpmethod = "exp" 
    } else if (family == "binomial") {
      rpmethod = "class" 
    } else if (family == "gaussian") {
      rpmethod = "anova"
    }

    xbetas0 = matrix(rep(xbnull,9*nobs), nrow=nobs, ncol=9)
    xbetas0.cv.rpart = matrix(rep(xbnull,9*nobs), nrow=nobs, ncol=9)
    xbnames = c("RPART c=0", "RPART c=0.01", "RPART c=0.02",
                "RPART fe", "RPART fe 0.01", "RPART fe 0.02",
                "RPART of", "RPART of 0.01", "RPART of 0.02")
    colnames( xbetas0 ) = xbnames 
    colnames( xbetas0.cv.rpart ) = xbnames 
    
    if (sum(ensemble[c(1,5)]) >= 1) {
      set.seed(seedr_)
      rpart.fit.00 <- rpart( y__ ~ ., data=df.xs, method=rpmethod, cp = 0 )
      rpart.fit.01 <- prune(rpart.fit.00, cp = 0.01)
      rpart.fit.02 <- prune(rpart.fit.00, cp = 0.02)
      rp_perf1 = rpart_perform(rpart.fit.00 , df.xs, NULL, y__, family )
      rp_perf2 = rpart_perform(rpart.fit.01 , df.xs, NULL, y__, family )
      rp_perf3 = rpart_perform(rpart.fit.02 , df.xs, NULL, y__, family )
      xbetas0[,1] = rpart_xbhat(rpart.fit.00 , df.xs, NULL, family )
      xbetas0[,2] = rpart_xbhat(rpart.fit.01 , df.xs, NULL, family )
      xbetas0[,3] = rpart_xbhat(rpart.fit.02 , df.xs, NULL, family )
    } else { rp_perf1 = c(1,1,0,0,0,0) ; rp_perf2 = c(1,1,0,0,0,0) ; rp_perf3 = c(1,1,0,0,0,0) } 
    
    if (sum(ensemble[c(2,6)]) >= 1) { 
      set.seed(seedr_)
      rpart.fitF.00 <- rpart( y__ ~ ., data=df.xsO, method=rpmethod, cp = 0 )
      rpart.fitF.01 <- prune(rpart.fitF.00, cp = 0.01)
      rpart.fitF.02 <- prune(rpart.fitF.00, cp = 0.02)
      rp_perf4 = rpart_perform(rpart.fitF.00 , df.xsO, NULL, y__, family )
      rp_perf5 = rpart_perform(rpart.fitF.01 , df.xsO, NULL, y__, family )
      rp_perf6 = rpart_perform(rpart.fitF.02 , df.xsO, NULL, y__, family )
      xbetas0[,4] = rpart_xbhat(rpart.fitF.00 , df.xsO, NULL, family )
      xbetas0[,5] = rpart_xbhat(rpart.fitF.01 , df.xsO, NULL, family )
      xbetas0[,6] = rpart_xbhat(rpart.fitF.02 , df.xsO, NULL, family )
    } else { rp_perf4 = c(1,1,0,0,0,0) ; rp_perf5 = c(1,1,0,0,0,0) ; rp_perf6 = c(1,1,0,0,0,0) } 
    
    if ((sum(ensemble[c(3,4,7,8)]) >= 1) & (!(family %in% c("binomial")))) {
      set.seed(seedr_)
      xsnames = names(df.xs)
      form1 = formula( paste( "y__ ~ ", paste(xsnames, collapse = " + " ), " + offset(ofst)" ) ) 
      rpart.fitO.00 <- rpart( form1 , data=df.xsO, method=rpmethod, cp = 0 )
      rpart.fitO.01 <- prune(rpart.fitO.00, cp = 0.01)
      rpart.fitO.02 <- prune(rpart.fitO.00, cp = 0.02)
      rp_perf7 = rpart_perform(rpart.fitO.00 , df.xsO, ofst, y__, family )
      rp_perf8 = rpart_perform(rpart.fitO.01 , df.xsO, ofst, y__, family )
      rp_perf9 = rpart_perform(rpart.fitO.02 , df.xsO, ofst, y__, family )
      xbetas0[,7] = rpart_xbhat(rpart.fitO.00 , df.xsO, ofst, family )
      xbetas0[,8] = rpart_xbhat(rpart.fitO.01 , df.xsO, ofst, family )
      xbetas0[,9] = rpart_xbhat(rpart.fitO.02 , df.xsO, ofst, family )
    } else { rp_perf7 = c(1,1,0,0,0,0) ; rp_perf8 = c(1,1,0,0,0,0) ; rp_perf9 = c(1,1,0,0,0,0) } 
    
    rpart.devian.naive     = c( rp_perf1[1] , rp_perf2[1] , rp_perf3[1] , rp_perf4[1] , rp_perf5[1] , rp_perf6[1] , rp_perf7[1] , rp_perf8[1] , rp_perf9[1] )
    rpart.cal.devian.naive = c( rp_perf1[2] , rp_perf2[2] , rp_perf3[2] , rp_perf4[2] , rp_perf5[2] , rp_perf6[2] , rp_perf7[2] , rp_perf8[2] , rp_perf9[2] )
    rpart.agree.naive      = c( rp_perf1[3] , rp_perf2[3] , rp_perf3[3] , rp_perf4[3] , rp_perf5[3] , rp_perf6[3] , rp_perf7[3] , rp_perf8[3] , rp_perf9[3] )
    rpart.intcal.naive     = c( rp_perf1[4] , rp_perf2[4] , rp_perf3[4] , rp_perf4[4] , rp_perf5[4] , rp_perf6[4] , rp_perf7[4] , rp_perf8[4] , rp_perf9[4] )
    rpart.lincal.naive     = c( rp_perf1[5] , rp_perf2[5] , rp_perf3[5] , rp_perf4[5] , rp_perf5[5] , rp_perf6[5] , rp_perf7[5] , rp_perf8[5] , rp_perf9[5] )
    rpart.nzero            = c( rp_perf1[6] , rp_perf2[6] , rp_perf3[6] , rp_perf4[6] , rp_perf5[6] , rp_perf6[6] , rp_perf7[6] , rp_perf8[6] , rp_perf9[6] )

    nms = c("cp=0.00", "cp=0.01", "cp=0.02", "F cp=0.00", "F cp=0.01", "F cp=0.02", "O cp=0.00", "O cp=0.01", "O cp=0.02" ) 
    names(rpart.devian.naive)  = nms 
    names(rpart.cal.devian.naive) = nms 
    names(rpart.agree.naive )  = nms
    names(rpart.intcal.naive)  = nms
    names(rpart.lincal.naive)  = nms
    names(rpart.nzero       )  = nms
    
    xbetas0.rpart = xbetas0
    
    if (track >= 1) { time_last = diff_time(time_start, time_last) } 
  }
  
  ###############################################################################################################
  ##### Neural Network fit ######################################################################################
  if ( (doann_ == 1) & (family == "cox") & (!is.null(start))) { 
    cat( "  The neural network fitting routine does not fit Cox model with (start,stop) time data\n",
         "  neural network models will not be fit\n")
    doann_ = 0
  }
  
  if (doann_ == 1) { 
    if (track >= 1) { cat(paste0("\n", " ########## Initial Neural Network fit on all data #############################################" , "\n")) }
    
    ann_perform = function(object=ann_fit_1_f, newdat=xs_z, newy=y_, family="binomial", start=NULL, event=NULL, tol=tol_) {
      if (family=="cox") { 
        SURV = Surv(newy, event)
        xbhat  = as.numeric( object$model(newdat) )  
        perf_cox(SURV, xbhat) 
      } else if (family=="binomial") { 
        phat_nn  = as.numeric( object$model(newdat) )  
        phat_nnt = phat_nn 
        phat_nnt[(phat_nn < tol)] = tol
        phat_nnt[(phat_nn > (1-tol))] = 1 - tol
        xbhat  = log(phat_nnt/(1-phat_nnt))
        perf_bin(newy, xbhat) 
      } else if (family=="gaussian") { 
        xbhat  = as.numeric( object$model(newdat) ) 
        perf_gau(newy, xbhat) 
      }
    }

    ann_predict = function(object=ann_fit_1_f, newdat=xs_z, family="binomial", tol=tol_) {
      if (family=="cox") { 
        xbhat  = as.numeric( object$model(newdat) )  
      } else if (family=="binomial") { 
        phat_nn  = as.numeric( object$model(newdat) )  
        phat_nnt = phat_nn 
        phat_nnt[(phat_nn < tol)] = tol
        phat_nnt[(phat_nn > (1-tol))] = 1 - tol
        xbhat  = log(phat_nnt/(1-phat_nnt))
      } else if (family=="gaussian") { 
        xbhat  = as.numeric( object$model(newdat) ) 
      }
      return(xbhat)
    }
        
    ann_fit_1_f = NULL ;  ann_fit_2_f = NULL ;  ann_fit_3_f = NULL ;  ann_fit_4_f = NULL ; ann_fit_5_f = NULL ; ann_fit_6_f = NULL ; ann_fit_7_f = NULL ; ann_fit_8_f = NULL ; 
    
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
    
    ## calibrate lasso ##
    if (sum(ensemble[c(2:8)]) > 0) {
      lassopred0 = predict(cv_glmnet_fit_f, xs, lam="lambda.min" , gam="gamma.min" )
      lassobeta  = predict(cv_glmnet_fit_f,     lam="lambda.min" , gam="gamma.min" )
      if (family=="cox") {
        fit0 = coxph(Surv(y_, event) ~ lassopred0) 
        lassopred = predict(fit0, as.data.frame(lassopred0))
        calbeta = c(mean(lassopred)-fit0$coefficients*mean(lassopred0), fit0$coefficients)
        #        summary( calbeta[1] + calbeta[2] * lassopred0 - lassopred)
      } else if (family %in% c("binomial", "gaussian")) {
        fit0 = glm(y_ ~ lassopred0, family=family) 
        lassopred = predict(fit0, as.data.frame(lassopred0))
        calbeta = fit0$coefficients
      }
      xs_z0L = cbind(lasso=lassopred, xs_z0)                                    ## add lasso prediction as first column 
      dim(xs_z0L)
      dim(xs_z)
      if (family == "cox") {
        xs_z1 = as.matrix( xs_z[,(lassobeta[[1]] != 0)]                           )         ## pick up non zero features 
      } else if (family %in% c("binomial", "gaussian")) {
        xs_z1 = as.matrix( xs_z[,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])]] )         ## pick up non zero features, remove intercept from this list 
      } 
      dim(xs_z1) 
      xs_z1L = cbind(lasso=lassopred,xs_z1)                                     ## add lasso prediction as first column 
      dim(xs_z1L)
      #      table(diag(cov(xs_z1L)))
      names(calbeta) = c("lassoB0", "lassoB1") 
      ann_zb$lassobeta=lassobeta ; ann_zb$calbeta=calbeta ; 
    }
    
    getwd = 0 ; getwd2 = 0 ; 
    if (folds_n >= 3) {
      if (wd  < 0) { getwd  = 1 ; wd  = abs(wd ) * cv_ridge_fit_f$lambda.min } 
      if (wd2 < 0) { getwd2 = 1 ; wd2 = abs(wd2) * cv_ridge_fit_f$lambda.min }
    } else {
      if (wd  < 0) { getwd  = 1 ; wd  = abs(wd ) } 
      if (wd2 < 0) { getwd2 = 1 ; wd2 = abs(wd2) }
    }
    
    getl1 = 0 ; getl12 = 0 ; 
    if (l1  < 0) { getl1  = 1 ; l1  = abs(l1 ) * cv_glmnet_fit_f$relaxed$lambda.min.g0 } 
    if (l12 < 0) { getl12 = 1 ; l12 = abs(l12) * cv_glmnet_fit_f$relaxed$lambda.min.g0 }
    
    ## 1 - uninformed NN       
    ## 2 - lasso feature
    ## 3 - lasso weights
    ## 4 - lasso updated weights
    ## 5 - lasso terms
    ## 6 - lasso terms & feat
    ## 7 - lasso terms & weights 
    ## 8 - lasso terms & updated weights
    
    xbetas0 = matrix(rep(xbnull,8*nobs), nrow=nobs, ncol=8)
    xbetas0.cv.ann = matrix(rep(xbnull,8*nobs), nrow=nobs, ncol=8)
    xbnames = c("ANN", "ANN.feat", "ANN.offs1", "ANN.offs2", 
                "ANN.lasso", "ANN.l.feat", "ANN.l.offs1", "ANN.l.offs2")
    colnames(xbetas0) = xbnames 
    colnames(xbetas0.cv.ann) = xbnames 

    if (family %in% c("cox", "binomial", "gaussian")) {
      #     cat(paste(fold_n, family, epochs, eppr, wd, lenz1, lenz2, mylr) )
      if (ensemble[1] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN uninformed **\n")) }
        ## xs_z0 - standardized data 
        ann_fit_1_f = ann_tab_cv_best(myxs=xs_z0, mystart=start, myy=y_, myevent=event, fold_n=folds_ann_n, family=family, 
                                      epochs=epochs, eppr=eppr, lenz1=lenz1, lenz2=lenz2, mylr=mylr, actv=actv, 
                                      drpot=drpot, wd=wd, l1=l1, scale=scale, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                      seed = c(seedr_, seedt_) ) 
        ann_perf1 = ann_perform(ann_fit_1_f, xs_z0, y_, family, start, event)
        xbetas0[,1] = ann_predict(ann_fit_1_f, xs_z0, family)
      } else { ann_perf1 = c(1,1,0,0,0) }
      
      if (ensemble[2] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso feature **\n")) }
        lenz1_ = lenz1 + 1 ; lenz2_ = lenz2 + 1 ;
        ann_fit_2_f = ann_tab_cv_best(myxs=xs_z0L, mystart=start, myy=y_, myevent=event, fold_n=folds_ann_n, family=family, 
                                      epochs=epochs, eppr=eppr, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr, actv=actv, 
                                      drpot=drpot, wd=wd, l1=l1, scale=scale, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                      seed = c(seedr_, seedt_)) 
        ann_perf2 = ann_perform(ann_fit_2_f, xs_z0L, y_, family, start, event)
        xbetas0[,2] = ann_predict(ann_fit_2_f, xs_z0L, family)
      } else { ann_perf2 = c(1,1,0,0,0) }
      
      if (ensemble[3] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso weights **\n")) }
        lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
        ann_fit_3_f = ann_tab_cv_best(myxs=xs_z0L, mystart=start, myy=y_, myevent=event, fold_n=folds_ann_n, family=family, 
                                      epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, 
                                      drpot=drpot, wd=wd2, l1=l12, lasso=1, lscale=lscale, scale=scale, resetlw=0, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                      seed = c(seedr_, seedt_)) 
        ann_perf3 = ann_perform(ann_fit_3_f, xs_z0L, y_, family, start, event)
        xbetas0[,3] = ann_predict(ann_fit_3_f, xs_z0L, family)
      } else { ann_perf3 = c(1,1,0,0,0) }
      
      if (ensemble[4] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso weights updated each epoch **\n")) }
        lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
        ann_fit_4_f = ann_tab_cv_best(myxs=xs_z0L, mystart=start, myy=y_, myevent=event, fold_n=folds_ann_n, family=family, 
                                      epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, 
                                      drpot=drpot, wd=wd2, l1=l12, lasso=1, lscale=lscale, scale=scale, resetlw=1, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                      seed = c(seedr_, seedt_)) 
        ann_perf4 = ann_perform(ann_fit_4_f, xs_z0L, y_, family, start, event)
        xbetas0[,4] = ann_predict(ann_fit_4_f, xs_z0L, family)
      } else { ann_perf4 = c(1,1,0,0,0) }
      
      if (ensemble[5] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso terms **\n")) }
        ## xs_z0 - standardized data 
        ann_fit_5_f = ann_tab_cv_best(myxs=xs_z1, mystart=start, myy=y_, myevent=event, fold_n=folds_ann_n, family=family, 
                                      epochs=epochs, eppr=eppr, lenz1=lenz1, lenz2=lenz2, mylr=mylr, actv=actv, 
                                      drpot=drpot, wd=wd, l1=l1, scale=scale, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                      seed = c(seedr_, seedt_)) 
        ann_perf5 = ann_perform(ann_fit_5_f, xs_z1, y_, family, start, event)
        xbetas0[,5] = ann_predict(ann_fit_5_f, xs_z1, family)
      } else { ann_perf5 = c(1,1,0,0,0) }
      
      if (ensemble[6] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso terms, lasso feature **\n")) }
        lenz1_ = lenz1 + 1 ; lenz2_ = lenz2 + 1 ;
        ann_fit_6_f = ann_tab_cv_best(myxs=xs_z1L, mystart=start, myy=y_, myevent=event, fold_n=folds_ann_n, family=family, 
                                      epochs=epochs, eppr=eppr, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr, actv=actv, 
                                      drpot=drpot, wd=wd, l1=l1, scale=scale, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                      seed = c(seedr_, seedt_)) 
        ann_perf6 = ann_perform(ann_fit_6_f, xs_z1L, y_, family, start, event)
        xbetas0[,6] = ann_predict(ann_fit_6_f, xs_z1L, family)
      } else { ann_perf6 = c(1,1,0,0,0) }
      
      if (ensemble[7] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso terms, lasso weights **\n")) }
        lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
        ann_fit_7_f = ann_tab_cv_best(myxs=xs_z1L, mystart=start, myy=y_, myevent=event, fold_n=folds_ann_n, family=family, 
                                      epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, 
                                      drpot=drpot, wd=wd2, l1=l12,lasso=1, lscale=lscale, scale=scale, resetlw=0, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                      seed = c(seedr_, seedt_)) 
        ann_perf7 = ann_perform(ann_fit_7_f, xs_z1L, y_, family, start, event)
        xbetas0[,7] = ann_predict(ann_fit_7_f, xs_z1L, family)
      } else { ann_perf7 = c(1,1,0,0,0) }
      
      if (ensemble[8] == 1) { 
        if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso terms, lasso weights updated each epoch **\n")) }
        lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
        ann_fit_8_f = ann_tab_cv_best(myxs=xs_z1L, mystart=start, myy=y_, myevent=event, fold_n=folds_ann_n, family=family, 
                                      epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, 
                                      drpot=drpot, wd=wd2, l1=l12,lasso=1, lscale=lscale, scale=scale, resetlw=1, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                      seed = c(seedr_, seedt_)) 
        ann_perf8 = ann_perform(ann_fit_8_f, xs_z1L, y_, family, start, event)
        xbetas0[,8] = ann_predict(ann_fit_8_f, xs_z1L, family)
      } else { ann_perf8 = c(1,1,0,0,0) }
      
      if (ensemble[9] == 1.41456) {
        if (eppr >= -2) { cat(paste0("  ** fitting ANN including X*Beta as offset - DOESN'T WORK properly **\n")) } 
        ann_fit_9_f = ann_tab_cv_best(myxs=xs_z1L, mystart=start, myy=y_, myevent=event, fold_n=folds_ann_n, family=family, 
                                      epochs=epochs2, eppr=eppr2, lenz1=lenz1, lenz2=lenz2, mylr=mylr2, actv=actv, 
                                      drpot=drpot, wd=wd, l1=l1, scale=scale, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                      seed = c(seedr_, seedt_)) 
        ann_perf9 = ann_perform(ann_fit_9_f, xs_z1, y_, family, start, event)   ## offset=lassopred, 
      } else { ann_perf9 = c(1,1,0,0,0)  }
      
      ann.devian.naive  = c( ann_perf1[1], ann_perf2[1], ann_perf3[1], ann_perf4[1], ann_perf5[1], ann_perf6[1], ann_perf7[1], ann_perf8[1] )
      ann.cal.devian.naive=c(ann_perf1[2], ann_perf2[2], ann_perf3[2], ann_perf4[2], ann_perf5[2], ann_perf6[2], ann_perf7[2], ann_perf8[2] )    
      ann.agree.naive   = c( ann_perf1[3], ann_perf2[3], ann_perf3[3], ann_perf4[3], ann_perf5[3], ann_perf6[3], ann_perf7[3], ann_perf8[3] )
      ann.intcal.naive  = c( ann_perf1[4], ann_perf2[4], ann_perf3[4], ann_perf4[4], ann_perf5[4], ann_perf6[4], ann_perf7[4], ann_perf8[4] )
      ann.lincal.naive  = c( ann_perf1[5], ann_perf2[5], ann_perf3[5], ann_perf4[5], ann_perf5[5], ann_perf6[5], ann_perf7[5], ann_perf8[5] )
      ann.nzero         = c( rep(dim(xs)[2],4), rep(lasso.nzero[4], 4) )
      
      nms = c("Uninformed", "lasso feat", "lasso w's", "lasso update", "lasso terms", "l/lasso feat", "l/lasso w's", "l/lasso update") 
      names(ann.devian.naive)  = nms 
      names(ann.cal.devian.naive) = nms
      names(ann.agree.naive)   = nms
      names(ann.intcal.naive)  = nms
      names(ann.lincal.naive)  = nms
      
      xbetas0.ann = xbetas0 
    } 
    
    if (track >= 1) { time_last = diff_time(time_start, time_last)  }
  }
  
  if (doann_  == 1) { xbetas = cbind(xbetas, xbetas0.ann) }
  if (dorpart == 1) { xbetas = cbind(xbetas, xbetas0.rpart) }
  
  ##### STEPWISE and AIC fits #################################################################################    
  
  if ((dostep == 1) | (doaic == 1)) {  
    step.devian.naive = c(1,1,1)
    step.cal.devian.naive = c(1,1,1)
    step.agree.naive   = c(0,0,0)
#    step.intcal.naive  = c(0,0,0)
#    step.lincal.naive  = c(0,0,0)
    step.nzero         = c(0,0,0)
    
    xbetas0 = matrix(rep(xbnull,3*nobs), nrow=nobs, ncol=3)
    xbetas0.cv.step = matrix(rep(xbnull,3*nobs), nrow=nobs, ncol=3)
    colnames(xbetas0)         = c("step.df", "step.p", "aic")
    colnames(xbetas0.cv.step) = c("step.df", "step.p", "aic")
  }

  ##### STEPWISE fit ##########################################################################################    
  if (dostep == 1) { 
    
    if (track >= 1) { cat(paste0("\n", " ########## Initial stepwise model on all data ########################################" , "\n")) }
    cv.stepreg.fit.all  = cv.stepreg(xs, start, y_, event, steps_n, folds_n, method=method, family=family,foldid=foldid,track=track)
    
    #### stepwise tuned by nzero (df) ##########################################
    func.fit.df = cv.stepreg.fit.all$func.fit.df                              
    nonzero      = length(func.fit.df$coefficients)                        
    if (family != "cox") {  nonzero =  nonzero - 1 }
    xb  = predict( func.fit.df, as.data.frame( xs) )              
    step_perf1 = c(perf_gen(y__, xb, family), nonzero )
    xbetas0[,1] = xb
    
    #### stepwise tuned by p (critical value) ##################################
    func.fit.p = cv.stepreg.fit.all$func.fit.p                                
    nonzero     = length(func.fit.p$coefficients)                          
    if (family != "cox") {  nonzero =  nonzero - 1 }
    xb     = predict( func.fit.p, as.data.frame( xs) )          
    step_perf2 = c(perf_gen(y__, xb, family), nonzero )
    xbetas0[,2] = xb
    
    step.devian.naive[c(1,2)] = c( step_perf1[1], step_perf2[1] )
    step.agree.naive[c(1,2)]  = c( step_perf1[3], step_perf2[3] )
    step.nzero[c(1,2)]        = c( step_perf1[6], step_perf2[6] )
    step.p                    = cv.stepreg.fit.all$best.p 
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

    predaic = predict(func.fit.aic, as.data.frame(xs))   
    aic_perf1 = perf_gen( y__ , predaic   , family )
    if (family != "cox") { nonzero = func.fit.aic$rank -1
    } else { nonzero = length(func.fit.aic$coefficients) }
    xbetas0[,3] = predaic
    
    step.devian.naive[c(3)]     = c( aic_perf1[1] )
    step.cal.devian.naive[c(3)] = c( aic_perf1[2] )
    step.agree.naive[c(3)]      = c( aic_perf1[3] )
    step.nzero[3]  = nonzero 

    if (track >= 1) { time_last = diff_time(time_start, time_last) } 
  }
  
  if ((dostep == 1) | (doaic == 1)) { xbetas = cbind(xbetas, xbetas0) }
  
  ##### track time before outer nested loop ###############################################################
  if (track >= 2) { 
    cat(paste0("\n", " ########## full data set analysis completed #####################################" , "\n")) 
  } 

  #########################################################################################################
  ##### FOLDS #############################################################################################
  ##### FOLDS #############################################################################################
  
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
      df.trainxs = as.data.frame(trainxs) 
      df.testxs  = as.data.frame(testxs) 
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
        if (is.null(start)) { 
          trainy__ = Surv(trainy_, trainevent)
          testy__  = Surv(testy_ , testevent)
        } else { 
          trainy__ = Surv(trainstart, trainy_, trainevent) 
          testy__  = Surv(teststart , testy_ , testevent ) 
        }
      } else {
        trainy__ = trainy_
        testy__  = testy_ 
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
      
      ##### get seeds for CV foldids and pass to CV functions #####
      
      seedr_ = seed$seedr[i_+1] 
      if (is.na(seedr_)) {  
        cat(paste( " seedr_ should already be set for CV ", i_, "\n")) 
        seedr_ = round(runif(1)*1e9)       ## duplicate 
        if (seedr_ == 0) { seedr_ = 1 } 
        seed$seedr[i_+1] = seedr_          ## duplicate 
      }
      seedt_ = seed$seedt[i_+1] 
      if (is.na(seedt_)) { 
        cat(paste( " seedt_ should already be set for CV ", i_, "\n")) 
        seedt_ = round(runif(1)*1e9)       ## duplicate 
        if (seedt_ == 0) { seedt_ = 1 } 
        seed$seedt[i_+1] = seedt_          ## duplicate 
      }

      ##### get CV foldid #####
      set.seed(seedr_) 
      foldidcv = get.foldid(trainy_, trainevent, family, folds_n, stratified) 
      table(foldidcv)   

      ##### get CV foldid for XGB #####      
      if (doxgb_ == 1) { 
        if (folds_xgb_n == folds_n) {
          foldidcv_xgb_v = foldidcv
        } else {
          if ( (folds_n / folds_xgb_n) == 2 ) {
            foldidcv_xgb_v = (foldidcv + 1) %/% 2 
          } else if ( (folds_n / folds_xgb_n) == 3 ) {
            foldidcv_xgb_v = (foldidcv + 1) %/% 3 
          } else {
            set.seed(seedr_) 
            foldidcv_xgb_v = get.foldid(trainy_, trainevent, family, folds_xgb_n, stratified) 
          } 
        }
        foldidcv_xgb = list()
        for (i1 in c(1:folds_xgb_n)) { foldidcv_xgb[[i1]] = c(1:length(trainy_))[foldidcv_xgb_v == i1] } 
        table( foldidcv_xgb_v, foldidcv)
        doxgb$folds = foldidcv_xgb             ## pass to xgb_tuned(), should be overwritten by folds=foldidcv_xgb, but just in case...
      }
      
      ##### get CV foldid for ANN #####      
      if (doann_ == 1) { 
       if (folds_ann_n == folds_n) {
          foldidcv_ann = foldidcv
        } else {
          if ( (folds_n / folds_ann_n) == 2 ) {
            foldidcv_ann = (foldidcv + 1) %/% 2 
          } else if ( (folds_n / folds_ann_n) == 3 ) {
            foldidcv_ann = (foldidcv + 1) %/% 3 
          } else {
            set.seed(seedr_) 
            foldidcv_ann = get.foldid(trainy_, trainevent, family, folds_ann_n, stratified) 
          } 
        }
        table( foldidcv_ann, foldidcv)
      } 
      
      ##### get null deviance for individual folds #####
      if (family == "cox") { 
        if (is.null(start)) { fit0 = coxph( Surv(testy_, testevent) ~ 1, ties=ties )
        } else              { fit0 = coxph( Surv(teststart, testy_, testevent) ~ 1, ties=ties ) } 
        null.m2LogLik.cv[i_] =  -2*fit0$loglik[1] / fit0$nevent
        if (ties == "breslow") {
          sat.m2LogLik.cv[i_] = - 2 * cox.sat.dev(testy_, testevent)[2] / fit0$nevent
        } else {
          sat.m2LogLik.cv[i_] = - 2 * cox.sat.dev(testy_, testevent)[1] / fit0$nevent
        }
        n.cv[i_] = fit0$nevent 
      } else if (family == "binomial") { 
        fit0 = glm( testy_ ~ 1 , family=family)
        null.m2LogLik.cv[i_] = fit0$null.deviance / length( testy_ )
        sat.m2LogLik.cv[i_] = 0 
        n.cv[i_] = length( testy_ )
      } else if (family == "gaussian") { 
        null.m2LogLik.cv[i_] = var(testy_) * (length(testy_)-1) / length(testy_)
#        var(testy_) * (length(testy_)-1) / length(testy_)
#        mean(testy_^2) - mean(testy_)^2
        sat.m2LogLik.cv[i_] = 0 
        n.cv[i_] = length( testy_ )
      }
      
      if (family == "cox") { testxbnull = 0 
      } else if (family == "binomial") {
        p = mean(testy_) 
        if (p > (1-tol_)) { p = 1/tol_
        } else if (p < tol_) { p = tol_ }
        testxbnull = log(p/(1-p))
      } else {
        testxbnull = mean(testy_) 
      }
      xbetas.cv[(foldid==i_),1] = testxbnull
      
      ##### LASSO fit ##########################################################################################
      if (dolasso==1) { 
        if (track >= 1) { cat( " ## Starting LASSO model fits ##", "\n" ) } 
        
        if (relax) { 
          #  xs=trainxs; start=trainstart; y_=trainy_; event=trainevent; lambda=lambda; gamma=gamma; folds_n=folds_n; limit=limit; fine=fine; track=1; family=family ;         
          cv_glmnet_fit = cv.glmnetr( trainxs, trainstart, trainy_, trainevent, family=family, 
                                      lambda=lambda, gamma=gamma, folds_n=folds_n, foldid=foldidcv, limit=limit, fine=fine, 
                                      track=0, ties=ties, ... )  
        } else {
          if (family=="cox") {   
            if ( is.null(start) ) { 
              cv_glmnet_fit = cv.glmnet( trainxs, Surv(trainy_, trainevent)            , nfolds=folds_n, foldid=foldidcv, family="cox", lambda=lambda, relax=FALSE, ... )  
            } else                { 
              cv_glmnet_fit = cv.glmnet( trainxs, Surv(trainstart, trainy_, trainevent), nfolds=folds_n, foldid=foldidcv, family="cox", lambda=lambda, relax=FALSE, ... )  
            }
          } else {
            cv_glmnet_fit = cv.glmnet( trainxs, trainy_, family=family, lambda=lambda, nfolds=folds_n, foldid=foldidcv, relax=FALSE, ... )  
          }
        }

        if (folds_n >= 3) {
          if ((family == "cox") & (is.null(trainstart))) {
            cv_ridge_fit = cv.glmnet( trainxs, Surv(trainy_, trainevent), family=family, alpha=0, nfolds=folds_n, foldid=foldidcv, ... )  
          } else if ((family == "cox") & (!is.null(trainstart))) {
            cv_ridge_fit = cv.glmnet( trainxs, Surv(trainstart, trainy_, trainevent), family=family, alpha=0, nfolds=folds_n, foldid=foldidcv, ... )  
          } else {
            cv_ridge_fit = cv.glmnet( trainxs, trainy_, family=family, alpha=0, nfolds=folds_n, foldid=foldidcv, ... )  
          }
        }
        
        lasso.nzero.cv[i_,1] = cv_glmnet_fit$nzero [ cv_glmnet_fit$index[2] ]
        lasso.nzero.cv[i_,2] = cv_glmnet_fit$nzero [ cv_glmnet_fit$index[1] ]     
        if (relax) {
          lasso.nzero.cv[i_,3] = cv_glmnet_fit$relaxed$nzero.1se 
          lasso.nzero.cv[i_,4] = cv_glmnet_fit$relaxed$nzero.min 
          lasso.nzero.cv[i_,5] = cv_glmnet_fit$nzero [ cv_glmnet_fit$relaxed$index.g0[2] ]
          lasso.nzero.cv[i_,6] = cv_glmnet_fit$nzero [ cv_glmnet_fit$relaxed$index.g0[1] ]
        }
        if (folds_n >= 3) { lasso.nzero.cv[i_,7] = cv_ridge_fit$nzero[ cv_ridge_fit$index ][1] 
        } else { lasso.nzero.cv[i_,7] = 0 }
        
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
        if (folds_n >= 3) { predridge = predict(cv_ridge_fit , testxs)
        } else { predridge = rep(0,dim(testxs)[1]) }
        
        perfm1 = perf_gen( testy__ , pred1se  , family )
        perfm2 = perf_gen( testy__ , predmin  , family )
        perfm3 = perf_gen( testy__ , pred1seR , family )
        perfm4 = perf_gen( testy__ , predminR , family )
        perfm5 = perf_gen( testy__ , pred1seR0, family )
        perfm6 = perf_gen( testy__ , predminR0, family )
        perfm7 = perf_gen( testy__ , predridge, family )
        lasso.devian.cv[i_,]  = c( perfm1[1] , perfm2[1] , perfm3[1] , perfm4[1] , perfm5[1] , perfm6[1] , perfm7[1] )
        lasso.cal.devian.cv[i_,]=c(perfm1[2] , perfm2[2] , perfm3[2] , perfm4[2] , perfm5[2] , perfm6[2] , perfm7[2] )    
        lasso.agree.cv[i_,]   = c( perfm1[3] , perfm2[3] , perfm3[3] , perfm4[3] , perfm5[3] , perfm6[3] , perfm7[3] )
        lasso.intcal.cv[i_,]  = c( perfm1[4] , perfm2[4] , perfm3[4] , perfm4[4] , perfm5[4] , perfm6[4] , perfm7[4] )
        lasso.lincal.cv[i_,]  = c( perfm1[5] , perfm2[5] , perfm3[5] , perfm4[5] , perfm5[5] , perfm6[5] , perfm7[5] )

#        xbetas0.cv.lasso[(foldid==i_),] = cbind(pred1se, predmin, pred1seR, predminR, pred1seR0, predminR0, predridge) 
        xbetas0.cv.lasso[(foldid==i_),1] = pred1se
        xbetas0.cv.lasso[(foldid==i_),2] = predmin
        xbetas0.cv.lasso[(foldid==i_),3] = pred1seR
        xbetas0.cv.lasso[(foldid==i_),4] = predminR
        xbetas0.cv.lasso[(foldid==i_),5] = pred1seR0
        xbetas0.cv.lasso[(foldid==i_),6] = predminR0
        xbetas0.cv.lasso[(foldid==i_),7] = predridge
        
        ## calibrate lasso ##
        if (sum(ensemble[c(2:4,6:8)]) >= 1) { 
          trainpredminR = predict(cv_glmnet_fit, trainxs, lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso as an offset #####
          testpredminR  = predminR 
          ofst0 = trainpredminR 
          if (family=="cox") {
            fit0 = coxph(Surv(trainy_, trainevent) ~ ofst0) 
          } else if (family %in% c("binomial", "gaussian")) {
            fit0 = glm(trainy_ ~ ofst0, family=family) 
          }
#          ofst0 = trainpredminR ; 
          trainofst = predict(fit0, as.data.frame(ofst0 )) 
          ofst0 = testpredminR ; 
          testofst  = predict(fit0, as.data.frame(ofst0 ))  
          trainxsO  = cbind(trainxs, ofst=trainofst) 
          testxsO   = cbind(testxs , ofst=testofst) 
          df.trainxsO = as.data.frame(trainxsO) 
          df.testxsO  = as.data.frame(testxsO ) 
          c(length(trainofst) , length(testofst) ) 
        }
        
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
        
        if (sum(ensemble[c(1,5)]) >= 1) {
          if (family=="cox") {
            train.xgb.dat <- xgb.DMatrix(data = trainxs, label = trainSurv.xgb)
            test.xgb.dat  <- xgb.DMatrix(data = testxs , label = testSurv.xgb )
          } else {
            train.xgb.dat <- xgb.DMatrix(data = trainxs, label = trainy_) 
            test.xgb.dat  <- xgb.DMatrix(data = testxs , label = testy_ ) 
          }
        }
        
        if (sum(ensemble[c(2,6)]) >= 1) { 
          if (family=="cox") {
            train.xgb.datF <- xgb.DMatrix(data = cbind(trainxs, ofst=trainofst), label = trainSurv.xgb)
            test.xgb.datF  <- xgb.DMatrix(data = cbind(testxs , ofst=testofst ), label = testSurv.xgb )
          } else {
            train.xgb.datF <- xgb.DMatrix(data = cbind(trainxs, ofst=trainofst), label = trainy_ ) 
            test.xgb.datF  <- xgb.DMatrix(data = cbind(testxs , ofst=testofst ), label = testy_  ) 
          }
        }
        
        if (sum(ensemble[c(3,4,7,8)]) >= 1) { 
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
          if (sum(ensemble[c(1,5)]) >= 1) {
            if (track >= 2) { cat("  XGB Simple, i_=", i_,"\n") }
            ## xgb_perform(xgb_model, xs_, y=y__, family, tol=tol_) 
            xgb.simple.train = xgb.simple(train.xgb.dat, objective = objective, seed=seedr_, folds=foldidcv_xgb, doxgb=doxgb, track=track )
            xgbperf1 = xgb_perform(xgb_model=xgb.simple.train, xs_=test.xgb.dat, y=testy__, family=family ) 
            xbetas0.cv.xgb[(foldid==i_),1] = xgb_xbhat(xgb.simple.train, test.xgb.dat, family)
          } else { xgbperf1 = c(1,1,0,0,0) }
          if (sum(ensemble[c(2,6)]) >= 1)  {
            if (track >= 2) { cat("  XGB Simple factor, i_=", i_,"\n") 
              if (track >= 4) { 
                print( summary(trainofst) ) 
                print( summary(testofst ) ) }
            } 
            xgb.simple.trainF = xgb.simple(train.xgb.datF, objective = objective, seed=seedr_, folds=foldidcv_xgb, doxgb=doxgb, track=track )
            xgbperf2 = xgb_perform(xgb.simple.trainF, test.xgb.datF, y=testy__, family=family) 
            xbetas0.cv.xgb[(foldid==i_),2] = xgb_xbhat(xgb.simple.trainF, test.xgb.datF, family)
          } else { xgbperf2 = c(1,1,0,0,0) }
          if (sum(ensemble[c(3,4,7,8)]) >= 1) {
            if (track >= 2) { cat("  XGB simple offset, i_=", i_,"\n") }
            xgb.simple.trainO = xgb.simple(train.xgb.datO, objective = objective, seed=seedr_, folds=foldidcv_xgb, doxgb=doxgb, track=track )
            xgbperf3 = xgb_perform(xgb.simple.trainO, test.xgb.datO, y=testy__, family=family) 
            xbetas0.cv.xgb[(foldid==i_),3] = xgb_xbhat(xgb.simple.trainO, test.xgb.datO, family)
          } else { xgbperf3 = c(1,1,0,0,0) }
          ## TUNED model fit train #######################################################
          if (sum(ensemble[c(1,5)]) >= 1) {
            if (track >= 2) { cat("\n  XGB tuned, i_=", i_) }
            xgb.tuned.train = xgb.tuned(train.xgb.dat, objective = objective, seed=seedr_, folds=foldidcv_xgb, doxgb=doxgb, track=track )
            xgbperf4 = xgb_perform(xgb.tuned.train, test.xgb.dat, y=testy__, family=family) 
            xbetas0.cv.xgb[(foldid==i_),4] = xgb_xbhat(xgb.tuned.train, test.xgb.dat, family)
          } else { xgbperf4 = c(1,1,0,0,0) }
          if (sum(ensemble[c(2,6)]) >= 1) {
            if (track >= 2) { cat("\n  XGB tuned factor, i_=", i_,"\n") }
            xgb.tuned.trainF = xgb.tuned(train.xgb.datF, objective = objective, seed=seedr_, folds=foldidcv_xgb, doxgb=doxgb, track=track )
            xgbperf5 = xgb_perform(xgb.tuned.trainF, test.xgb.datF, y=testy__, family=family) 
            xbetas0.cv.xgb[(foldid==i_),5] = xgb_xbhat(xgb.tuned.trainF, test.xgb.datF, family)
          } else { xgbperf5 = c(1,1,0,0,0) }
          if (sum(ensemble[c(3,4,7,8)]) >= 1) {
            if (track >= 2) { cat("\n  XGB tuned offset, i_=", i_,"\n") }
            xgb.tuned.trainO = xgb.tuned(train.xgb.datO, objective = objective, seed=seedr_, folds=foldidcv_xgb, doxgb=doxgb, track=track )
            xgbperf6 = xgb_perform(xgb.tuned.trainO, test.xgb.datO, y=testy__, family=family) 
            xbetas0.cv.xgb[(foldid==i_),6] = xgb_xbhat(xgb.tuned.trainO, test.xgb.datO, family)
          } else { xgbperf6 = c(1,1,0,0,0) }
          ################################################################################
          xgb.devian.cv[i_,]     = c( xgbperf1[1] , xgbperf2[1] , xgbperf3[1] , xgbperf4[1] , xgbperf5[1] , xgbperf6[1] )
          xgb.cal.devian.cv[i_,] = c( xgbperf1[2] , xgbperf2[2] , xgbperf3[2] , xgbperf4[2] , xgbperf5[2] , xgbperf6[2] )
          xgb.agree.cv[i_,]      = c( xgbperf1[3] , xgbperf2[3] , xgbperf3[3] , xgbperf4[3] , xgbperf5[3] , xgbperf6[3] )
          xgb.intcal.cv[i_,]     = c( xgbperf1[4] , xgbperf2[4] , xgbperf3[4] , xgbperf4[4] , xgbperf5[4] , xgbperf6[4] )
          xgb.lincal.cv[i_,]     = c( xgbperf1[5] , xgbperf2[5] , xgbperf3[5] , xgbperf4[5] , xgbperf5[5] , xgbperf6[5] )
          ################################################################################
        }
        
      xgb.nzero.cv[i_,]      = rep( xs_ncol , 6)
      doxgb$folds = foldid_xgb             ## foldidcv_xgb passed do xgb_.tuned(), should be returned to whole sample value 
      if (track >= 1) { time_last = diff_time(time_start, time_last) } 
        #      if (track >= 1) { cat(paste0(" (Continuing Nested Cross Validation outer fold  ", i_, "  of  " , folds_n , " )" , "\n")) }
      }
      
      ##### RandomForest fit ##########################################################################################
      if (dorf_==1) { 
        if (track >= 1) { cat(" ## Starting Random Forest model fits ##" , "\n") }
        
        if (sum(ensemble[c(1,5)]) >= 1) {
          rf_tuned_train = rf_tune(xs=trainxs, y_=trainy_, event=trainevent, family=family, mtryc=mtryc, ntreec = ntreec, seed=seedr_)
          rf_perf1 = rf_perform(rf_model=rf_tuned_train$rf_tuned, dframe=as.data.frame(testxs), ofst=NULL, y__=testy__, family=family, tol=tol_)
          xbetas0.cv.rf[(foldid==i_),1] = rf_xbhat(rf_model=rf_tuned_train$rf_tuned, dframe=as.data.frame(testxs), ofst=NULL, family=family, tol=tol_)
          if (track >= 2) { time_lasti = diff_time(time_start, time_last) } 
        } else { rf_perf1 = c(1,1,0,0,0,0) }
        
        if (sum(ensemble[c(2,6)]) >= 1) { 
          rf_tuned_trainF = rf_tune(xs=trainxsO, y_=trainy_, event=trainevent, family=family, mtryc=mtryc, ntreec = ntreec, seed=seedr_)
          rf_perf2 = rf_perform(rf_model=rf_tuned_trainF$rf_tuned, dframe=df.testxsO, ofst=NULL, y__=testy__, family=family, tol=tol_)
          xbetas0.cv.rf[(foldid==i_),2] = rf_xbhat(rf_model=rf_tuned_trainF$rf_tuned, dframe=df.testxsO, ofst=NULL, family=family, tol=tol_)
          if (track >= 2) { time_lasti = diff_time(time_start, time_lasti) } 
        } else { rf_perf2 = c(1,1,0,0,0,0) }
        
        if ( (sum(ensemble[c(3,4,7,8)]) >= 1) & (family %in% c("gaussian")) ) {
          trainyo = trainy_ - trainofst 
          c(var(trainyo), var(trainy__))
          rf_tuned_trainO = rf_tune(xs=trainxs, y_=trainyo, event=NULL, family=family, mtryc=mtryc, ntreec = ntreec, seed=seedr_)
          rf_perf3 = rf_perform(rf_model=rf_tuned_trainO$rf_tuned, dframe=df.testxs, ofst=testofst, y__=testy__, family=family, tol=tol_)
          xbetas0.cv.rf[(foldid==i_),3] = rf_xbhat(rf_model=rf_tuned_trainO$rf_tuned, dframe=df.testxs, ofst=testofst, family=family, tol=tol_)
          if (track >= 2) { time_lasti = diff_time(time_start, time_lasti) } 
        } else { rf_perf3 = c(1,1,0,0,0,0) }
        
        rf.devian.cv[i_,]     = c( rf_perf1[1] , rf_perf2[1] , rf_perf3[1] )
        rf.cal.devian.cv[i_,] = c( rf_perf1[2] , rf_perf2[2] , rf_perf3[2] )
        rf.agree.cv[i_,]      = c( rf_perf1[3] , rf_perf2[3] , rf_perf3[3] )
        rf.intcal.cv[i_,]     = c( rf_perf1[4] , rf_perf2[4] , rf_perf3[4] )
        rf.lincal.cv[i_,]     = c( rf_perf1[5] , rf_perf2[5] , rf_perf3[5] )
        rf.mtry.cv[i_,]       = c( rf_perf1[6] , rf_perf2[6] , rf_perf3[6] )
        
        if (track >= 1) { time_last = diff_time(time_start, time_last) } 
      }
      
      ##### RPART fit ##########################################################################################
      if (dorpart==1) { 
        if (track >= 1) { cat(" ## Starting RPART model fits ##" , "\n") }
        
        #  nms = c("cp=0.00", "cp=0.01", "cp=0.02", "F cp=0.00", "F cp=0.01", "F cp=0.02", "O cp=0.00", "O cp=0.01", "O cp=0.02" ) 

        if (family == "cox") {
          rpmethod = "exp" 
        } else if (family == "binomial") {
          rpmethod = "class" 
        } else if (family == "gaussian") {
          rpmethod = "anova"
        }
        
        if (sum(ensemble[c(1,5)]) >= 1) {
          #        colnames(df.train)[1] = c("trainy_") 
          set.seed(seedr_) 
          rpart.train.00 <- rpart( trainy__ ~ ., data=df.trainxs, method=rpmethod, cp = 0 )
          rpart.train.01 <- prune(rpart.train.00, cp = 0.01)
          rpart.train.02 <- prune(rpart.train.00, cp = 0.02)
          rp_perf1 = rpart_perform(rpart.train.00 , df.testxs, NULL, testy__, family )
          rp_perf2 = rpart_perform(rpart.train.01 , df.testxs, NULL, testy__, family )
          rp_perf3 = rpart_perform(rpart.train.02 , df.testxs, NULL, testy__, family )
          xbetas0.cv.rpart[(foldid==i_),1] = rpart_xbhat(rpart.train.00 , df.testxs, NULL, family )
          xbetas0.cv.rpart[(foldid==i_),2] = rpart_xbhat(rpart.train.01 , df.testxs, NULL, family )
          xbetas0.cv.rpart[(foldid==i_),3] = rpart_xbhat(rpart.train.02 , df.testxs, NULL, family )
        } else { rp_perf1 = rep(0,6) ; rp_perf2 = rep(0,6) ; rp_perf3 = rep(0,6) }
        
        if (sum(ensemble[c(2,6)]) >= 1)  { 
          set.seed(seedr_) 
          rpart.trainF.00 <- rpart( trainy__ ~ ., data=df.trainxsO, method=rpmethod, cp = 0 )
          rpart.trainF.01 <- prune(rpart.trainF.00, cp = 0.01)
          rpart.trainF.02 <- prune(rpart.trainF.00, cp = 0.02)
          rp_perf4 = rpart_perform(rpart.trainF.00 , df.testxsO, NULL, testy__, family )
          rp_perf5 = rpart_perform(rpart.trainF.01 , df.testxsO, NULL, testy__, family )
          rp_perf6 = rpart_perform(rpart.trainF.02 , df.testxsO, NULL, testy__, family )
          xbetas0.cv.rpart[(foldid==i_),4] = rpart_xbhat(rpart.trainF.00 , df.testxsO, NULL, family )
          xbetas0.cv.rpart[(foldid==i_),5] = rpart_xbhat(rpart.trainF.01 , df.testxsO, NULL, family )
          xbetas0.cv.rpart[(foldid==i_),6] = rpart_xbhat(rpart.trainF.02 , df.testxsO, NULL, family )
        } else { rp_perf4 = rep(0,6) ; rp_perf5 = rep(0,6) ; rp_perf6 = rep(0,6) }
        
        if ((sum(ensemble[c(3,4,7,8)]) >= 1) & (!(family %in% c("binomial"))))  {
          #        xsnames = names(df.trainxsO)
          xsnames = names(df.trainxs)
          form1 = formula( paste( "trainy__ ~ ", paste(xsnames, collapse = " + " ), " + offset(ofst)" ) ) 
          set.seed(seedr_) 
          rpart.trainO.00 <- rpart( form1 , data=df.trainxsO, method=rpmethod, cp = 0 )
          rpart.trainO.01 <- prune(rpart.trainO.00, cp = 0.01)
          rpart.trainO.02 <- prune(rpart.trainO.00, cp = 0.02)
          rp_perf7 = rpart_perform(rpart.trainO.00 , df.testxsO, testofst, testy__, family )
          rp_perf8 = rpart_perform(rpart.trainO.01 , df.testxsO, testofst, testy__, family )
          rp_perf9 = rpart_perform(rpart.trainO.02 , df.testxsO, testofst, testy__, family )
          xbetas0.cv.rpart[(foldid==i_),7] = rpart_xbhat(rpart.trainO.00 , df.testxsO, testofst, family )
          xbetas0.cv.rpart[(foldid==i_),8] = rpart_xbhat(rpart.trainO.01 , df.testxsO, testofst, family )
          xbetas0.cv.rpart[(foldid==i_),9] = rpart_xbhat(rpart.trainO.02 , df.testxsO, testofst, family )
        } else { rp_perf7 = rep(0,6) ; rp_perf8 = rep(0,6) ; rp_perf9 = rep(0,6) }

        rpart.devian.cv[i_,]     = c( rp_perf1[1] , rp_perf2[1] , rp_perf3[1] , rp_perf4[1] , rp_perf5[1] , rp_perf6[1] , rp_perf7[1] , rp_perf8[1] , rp_perf9[1] )
        rpart.cal.devian.cv[i_,] = c( rp_perf1[2] , rp_perf2[2] , rp_perf3[2] , rp_perf4[2] , rp_perf5[2] , rp_perf6[2] , rp_perf7[2] , rp_perf8[2] , rp_perf9[2] )
        rpart.agree.cv[i_,]      = c( rp_perf1[3] , rp_perf2[3] , rp_perf3[3] , rp_perf4[3] , rp_perf5[3] , rp_perf6[3] , rp_perf7[3] , rp_perf8[3] , rp_perf9[3] )
        rpart.intcal.cv[i_,]     = c( rp_perf1[4] , rp_perf2[4] , rp_perf3[4] , rp_perf4[4] , rp_perf5[4] , rp_perf6[4] , rp_perf7[4] , rp_perf8[4] , rp_perf9[4] )
        rpart.lincal.cv[i_,]     = c( rp_perf1[5] , rp_perf2[5] , rp_perf3[5] , rp_perf4[5] , rp_perf5[5] , rp_perf6[5] , rp_perf7[5] , rp_perf8[5] , rp_perf9[5] )
        rpart.nzero.cv[i_,]      = c( rp_perf1[6] , rp_perf2[6] , rp_perf3[6] , rp_perf4[6] , rp_perf5[6] , rp_perf6[6] , rp_perf7[6] , rp_perf8[6] , rp_perf9[6] )
        
        if (track >= 1) { time_last = diff_time(time_start, time_last) } 
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
        testxs_z0    = testxs_z [, (trainxs_sds_ > 0) ]                         ## add lasso prediction as first column 

        ## calibrate lasso ##
        if (sum(ensemble[c(2:8)]) > 0.5) {        
          lassopredtrain0 = predict(cv_glmnet_fit, trainxs, lam="lambda.min" , gam="gamma.min" )
          lassopredtest0  = predict(cv_glmnet_fit, testxs , lam="lambda.min" , gam="gamma.min" )
          lassobeta = predict(cv_glmnet_fit, lam="lambda.min" , gam="gamma.min" )
          lassopred0 = lassopredtrain0 
          if (family=="cox") {
            fit0 = coxph(Surv(trainy_, trainevent) ~ lassopred0) 
          } else if (family %in% c("binomial", "gaussian")) {
            fit0 = glm(trainy_ ~ lassopred0, family=family) 
          }
          lassopred0 = lassopredtrain0 
          lassopredtrain = predict(fit0, as.data.frame(lassopred0)) ;  length(lassopredtrain)
          lassopred0 = lassopredtest0 
          lassopredtest  = predict(fit0, as.data.frame(lassopred0)) ;  length(lassopredtest)
          trainxs_z0L = cbind(lasso=lassopredtrain,trainxs_z0)                  ## add lasso prediction as first column 
          testxs_z0L  = cbind(lasso=lassopredtest,testxs_z0)                    ## add lasso prediction as first column 
          dim(trainxs_z0L)
          dim(testxs_z0L)
          
          if (family=="cox") {
            trainxs_z1 = as.matrix( trainxs_z[,(lassobeta[[1]] != 0)[1:length(lassobeta[[1]])]] )   ## pick up non zero features, remove intercept from this list 
            testxs_z1  = as.matrix( testxs_z [,(lassobeta[[1]] != 0)[1:length(lassobeta[[1]])]] )   ## pick up non zero features, remove intercept from this list                         ## add lasso prediction as first column 
          } else if (family %in% c("binomial", "gaussian")) {
            trainxs_z1 = as.matrix( trainxs_z[,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])]] )   ## pick up non zero features, remove intercept from this list 
            testxs_z1  = as.matrix( testxs_z [,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])]] )   ## pick up non zero features, remove intercept from this list 
          }          
          trainxs_z1L = cbind(lasso=lassopredtrain,trainxs_z1)                        ## add lasso prediction as first column 
          testxs_z1L  = cbind(lasso=lassopredtest,testxs_z1) 
          
          # table(round(diag(cov(trainxs_z1L)),digits=8))
        }
        
        if (folds_n >= 3) {
          if (getwd  == 1) { wd  = cv_ridge_fit$lambda.min } 
          if (getwd2 == 1) { wd2 = cv_ridge_fit$lambda.min }
        }
        
        if (getl1  == 1) { l1  = cv_glmnet_fit$relaxed$lambda.min.g0 } 
        if (getl12 == 1) { l12 = cv_glmnet_fit$relaxed$lambda.min.g0 }
        
        if (family %in% c("cox", "binomial", "gaussian")) {
          if (ensemble[1] == 1) { 
            if (eppr >= -2) { cat(paste0("\n  ** fitting uninformed ANN **\n")) }
            ann_fit_1 = ann_tab_cv_best(myxs=trainxs_z0, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=folds_ann_n, family=family, 
                                        epochs=epochs, eppr=eppr, lenz1=lenz1, lenz2=lenz2, mylr=mylr, actv=actv, 
                                        drpot=drpot, wd=wd, l1=l1, scale=scale, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                        seed = c(seedr_, seedt_) ) 
            ann_perf1 = ann_perform(ann_fit_1, testxs_z0, testy_, family, teststart, testevent)
            xbetas0.cv.ann[(foldid==i_),1] = ann_predict(ann_fit_1, testxs_z0, family)
          } else { ann_perf1 = c(1,1,0,0,0) }
          
          if (ensemble[2] == 1) { 
            if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso feature **\n")) }
            lenz1_ = lenz1 + 1 ; lenz2_ = lenz2 + 1 ;
            ann_fit_2 = ann_tab_cv_best(myxs=trainxs_z0L, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=folds_ann_n, family=family, 
                                        epochs=epochs, eppr=eppr, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr, actv=actv, 
                                        drpot=drpot, wd=wd, l1=l1, scale=scale, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                        seed = c(seedr_, seedt_) ) 
            ann_perf2 = ann_perform(ann_fit_2, testxs_z0L, testy_, family, teststart, testevent)
            xbetas0.cv.ann[(foldid==i_),2] = ann_predict(ann_fit_2, testxs_z0L, family)
          } else { ann_perf2 = c(1,1,0,0,0) }
          
          if (ensemble[3] == 1) { 
            if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso weights **\n")) } 
            lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
            ann_fit_3 = ann_tab_cv_best(myxs=trainxs_z0L, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=folds_ann_n, family=family, 
                                        epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, 
                                        drpot=drpot, wd=wd, l1=l1, lasso=1, lscale=lscale, scale=scale, resetlw=0, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                        seed = c(seedr_, seedt_) ) 
            ann_perf3 = ann_perform(ann_fit_3, testxs_z0L, testy_, family, teststart, testevent)
            xbetas0.cv.ann[(foldid==i_),3] = ann_predict(ann_fit_3, testxs_z0L, family)
          } else { ann_perf3 = c(1,1,0,0,0) }
          
          if (ensemble[4] == 1) { 
            if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso weights updated each epoch **\n")) } 
            lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
            ann_fit_4 = ann_tab_cv_best(myxs=trainxs_z0L, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=folds_ann_n, family=family, 
                                        epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, 
                                        drpot=drpot, wd=wd2, l1=l12, lasso=1, lscale=lscale, scale=scale, resetlw=1, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                        seed = c(seedr_, seedt_) ) 
            ann_perf4 = ann_perform(ann_fit_4, testxs_z0L, testy_, family, teststart, testevent)
            xbetas0.cv.ann[(foldid==i_),4] = ann_predict(ann_fit_4, testxs_z0L, family)
          } else { ann_perf4 = c(1,1,0,0,0) }
          
          if (ensemble[5] == 1) { 
            if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso terms **\n")) }
            ann_fit_5 = ann_tab_cv_best(myxs=trainxs_z1, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=folds_ann_n, family=family, 
                                        epochs=epochs, eppr=eppr, lenz1=lenz1, lenz2=lenz2, mylr=mylr, actv=actv, 
                                        drpot=drpot, wd=wd, l1=l1, scale=scale, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                        seed = c(seedr_, seedt_) ) 
            ann_perf5 = ann_perform(ann_fit_5, testxs_z1, testy_, family, teststart, testevent)
            xbetas0.cv.ann[(foldid==i_),5] = ann_predict(ann_fit_5, testxs_z1, family)
          } else { ann_perf5 = c(1,1,0,0,0) }
          
          if (ensemble[6] == 1) { 
            if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso terms, lasso feature **\n")) }
            lenz1_ = lenz1 + 1 ; lenz2_ = lenz2 + 1 ;
            ann_fit_6 = ann_tab_cv_best(myxs=trainxs_z1L, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=folds_ann_n, family=family, 
                                        epochs=epochs, eppr=eppr, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr, actv=actv, 
                                        drpot=drpot, wd=wd, l1=l1, scale=scale, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                        seed = c(seedr_, seedt_) ) 
            ann_perf6 = ann_perform(ann_fit_6, testxs_z1L, testy_, family, teststart, testevent)
            xbetas0.cv.ann[(foldid==i_),6] = ann_predict(ann_fit_6, testxs_z1L, family)
          } else { ann_perf6 = c(1,1,0,0,0) }
          
          if (ensemble[7] == 1) { 
            if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso terms, lasso weights **\n")) } 
            lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
            ann_fit_7 = ann_tab_cv_best(myxs=trainxs_z1L, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=folds_ann_n, family=family, 
                                        epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, 
                                        drpot=drpot, wd=wd, l1=l1, lasso=1, lscale=lscale, scale=scale, resetlw=0, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                        seed = c(seedr_, seedt_) ) 
            ann_perf7 = ann_perform(ann_fit_7, testxs_z1L, testy_, family, teststart, testevent)
            xbetas0.cv.ann[(foldid==i_),7] = ann_predict(ann_fit_7, testxs_z1L, family)
          } else { ann_perf7 = c(1,1,0,0,0) }
          
          if (ensemble[8] == 1) { 
            if (eppr >= -2) { cat(paste0("  ** fitting ANN lasso terms, lasso weights updated each epoch **\n")) } 
            lenz1_ = lenz1 + 2 ; lenz2_ = lenz2 + 2 ; 
            ann_fit_8 = ann_tab_cv_best(myxs=trainxs_z1L, mystart=trainstart, myy=trainy_, myevent=trainevent, fold_n=folds_ann_n, family=family, 
                                        epochs=epochs2, eppr=eppr2, lenz1=lenz1_, lenz2=lenz2_, mylr=mylr2, actv=actv, 
                                        drpot=drpot, wd=wd2, l1=l12, lasso=1, lscale=lscale, scale=scale, resetlw=1, minloss=minloss, gotoend=gotoend, bestof=bestof, 
                                        seed = c(seedr_, seedt_) ) 
            ann_perf8 = ann_perform(ann_fit_8, testxs_z1L, testy_, family, teststart, testevent)
            xbetas0.cv.ann[(foldid==i_),8] = ann_predict(ann_fit_8, testxs_z1L, family)
          } else { ann_perf8 = c(1,1,0,0,0) }
          
          ann.devian.cv[i_,]  = c( ann_perf1[1], ann_perf2[1], ann_perf3[1], ann_perf4[1], ann_perf5[1], ann_perf6[1], ann_perf7[1], ann_perf8[1] )
          ann.cal.devian.cv[i_,]=c(ann_perf1[2], ann_perf2[2], ann_perf3[2], ann_perf4[2], ann_perf5[2], ann_perf6[2], ann_perf7[2], ann_perf8[2] )    
          ann.agree.cv [i_,]  = c( ann_perf1[3], ann_perf2[3], ann_perf3[3], ann_perf4[3], ann_perf5[3], ann_perf6[3], ann_perf7[3], ann_perf8[3] )
          ann.intcal.cv[i_,]  = c( ann_perf1[4], ann_perf2[4], ann_perf3[4], ann_perf4[4], ann_perf5[4], ann_perf6[4], ann_perf7[4], ann_perf8[4] )
          ann.lincal.cv[i_,]  = c( ann_perf1[5], ann_perf2[5], ann_perf3[5], ann_perf4[5], ann_perf5[5], ann_perf6[5], ann_perf7[5], ann_perf8[5] )
          ann.nzero.cv[i_,]   = c( rep(dim(xs)[2],4), rep(lasso.nzero.cv[i_,4], 4) )
        } 
        if (track == 1) { eppr = eppr_t ; eppr2 = eppr2_t }
        if (track >= 1) { time_last = diff_time(time_start, time_last) } 
      } 
      
      ###### STEPWISE fits #####################################################################################
      
      if (dostep == 1) { 
        if (track >= 1) { cat(paste0(" ## Starting STEPWISE model fits ## ", "\n")) }
        
        if (family == "cox") {
          if (is.null(start)) { 
            testy__ = Surv(testy_, testevent)
          } else {
            testy__ = Surv(teststart , testy_, testevent)
          } 
        } else {
          testy__ = testy_ 
        }

        ## xs_cv=trainxs ; start_cv=trainstart ; y_cv=trainy_ ; event_cv=trainevent ; steps_n=steps_n ; folds_n_cv=folds_n ; method=method ; family=family ; 
        
        cv.stepreg.fit = cv.stepreg(trainxs, trainstart, trainy_, trainevent, steps_n, folds_n, method=method, family=family, foldid=foldidcv, track=0) 
        stepreg.fit.all.best = cv.stepreg.fit$stepreg.fit.all.best 
        
        #### stepwise tuned by nzero ###########################################
        func.fit.df = cv.stepreg.fit$func.fit.df                              
        nonzero      = length(func.fit.df$coefficients)                        
        if (family != "cox") {  nonzero =  nonzero - 1 }
        testxb  = predict( func.fit.df, as.data.frame( testxs) )              
        step_perf1 = c(perf_gen(testy__, testxb, family) )
        xbetas0.cv.step[(foldid==i_),1] = testxb
        step.nzero.cv [i_,1] = nonzero 
        
        #### stepwise tuned by p (critical value) ##############################
        func.fit.p = cv.stepreg.fit$func.fit.p                                
        nonzero     = length(func.fit.p$coefficients)                          
        if (family != "cox") {  nonzero =  nonzero - 1 }
        testxb     = predict( func.fit.p, as.data.frame( testxs) )          
        step_perf2 = c(perf_gen(testy__, testxb, family) )
        xbetas0.cv.step[(foldid==i_),2] = testxb
        step.nzero.cv [i_,2] = nonzero 
        step.p.cv  [i_,1] = cv.stepreg.fit$best.p 
        
        ########################################################################
        
        step.devian.cv[i_,c(1,2)] = c(step_perf1[1] , step_perf2[1])
        step.cal.devian.cv[i_,c(1,2)] = c(step_perf1[2] , step_perf2[2])
        step.agree.cv[i_,c(1,2)]  = c(step_perf1[3] , step_perf2[3])
        step.intcal.cv[i_,c(1,2)] = c(step_perf1[4] , step_perf2[4])
        step.lincal.cv[i_,c(1,2)] = c(step_perf1[5] , step_perf2[5])
        
        if (track >= 1) { time_last = diff_time(time_start, time_last) } 
      }
      
      ###### AIC fits ##########################################################################################
      if (doaic==1) {
        if (track >= 1) { cat(" ## Starting AIC model fits ##" , "\n") }
        
        if (dostep==1) {
          stepreg.fit.best = as.matrix( cv.stepreg.fit$stepreg.fit.all.best ) 
        } else {
          stepreg.fit = stepreg(trainxs, trainstart, trainy_, trainevent, steps_n=steps_n)  
          class(stepreg.fit) = "data.frame" 
          stepreg.fit.best = stepreg.fit[(stepreg.fit$best==1) , ] 
        }
        aic = 2*stepreg.fit.best[,1] - 2*stepreg.fit.best[,5]                     
        nzero  = which.min( aic )                                                    ## may be different from first p > 0.15 XX
        
        if (family=="cox") {
          beta = stepreg.fit.best [ nzero , (xs_ncol+8):(2*xs_ncol+7) ]              ## get beta, the regression estimate
        } else {
          beta = stepreg.fit.best [ nzero , (xs_ncol+8):(2*xs_ncol+8) ]              ## get beta, the regression estimate
        }
        
        beta = as.numeric(beta)
        if (family=="cox") { testxb  = testxs  %*% beta
        } else             { testxb  = testxs1 %*% beta }                      ## get predicteds 
        
        if (family == "cox") {
          if (is.null(start)) { 
            testy__ = Surv(testy_, testevent)
          } else {
            testy__ = Surv(teststart , testy_, testevent)
          } 
        } else {
          testy__ = testy_ 
        }
        
        step_perf3 = c(perf_gen(testy__, testxb, family), nzero )
        xbetas0.cv.step[(foldid==i_),3] = testxb

        step.devian.cv[i_,c(3)]     = c( step_perf3[1] )
        step.cal.devian.cv[i_,c(3)] = c( step_perf3[2] )
        step.agree.cv[i_,c(3)]      = c( step_perf3[3] )
        step.intcal.cv[i_,c(3)]     = c( step_perf3[4] )
        step.lincal.cv[i_,c(3)]     = c( step_perf3[5] )
        step.nzero.cv [i_,3] = nzero 
        
        if (track >= 1) { time_last = diff_time(time_start, time_last) } 
      }
      ## end within fold AIC FIT #################################################
    }
  }
  
  ##### END FOLDS ##########################################################################################
  ##### END FOLDS ##########################################################################################
  ##########################################################################################################
  
  if (track >= 1) { cat(paste0("\n", " ################################################################################################" , "\n")) } 
  
  if ( family=="cox" ) { 
    nevents = sum( event )
    if (is.null(start)) { fit0 = coxph( Surv(y_, event) ~ 1, ties=ties )
    } else              { fit0 = coxph( Surv(start, y_, event) ~ 1, ties=ties ) } 
    if (ties == "breslow") {
      sat.m2LogLik = - 2 * cox.sat.dev(y_, event)[2] / fit0$nevent 
    } else {
      sat.m2LogLik = - 2 * cox.sat.dev(y_, event)[1] / fit0$nevent 
    }
    null.m2LogLik = - 2*(fit0$loglik[1])/ fit0$nevent 
    null.deviance = null.m2LogLik - sat.m2LogLik 
  } else if (family=="binomial") { 
    nevents=sum(y_) 
    fit0 = glm( y_ ~ 1 , family=family)
    null.m2LogLik = fit0$null.deviance / length( y_ )
    sat.m2LogLik = 0 
    null.deviance = null.m2LogLik
  } else { 
    nevents=NA     
    null.m2LogLik = (length(y_)-1)*var(y_)/length(y_) 
    sat.m2LogLik  = 0 
    null.deviance = null.m2LogLik
  } 

  time_stop=Sys.time()
  
  if (family == "cox") { 
    dep_names = c(start_name, y_name, event_name) 
    names(dep_names) = c("start","y_","event") 
    if ( (dep_names[1] == "NULL") & (dep_names[2] == "NULL") & (dep_names[3] == "NULL") ) { dep_names = NULL }
  } else { 
    dep_names = c(y_name) 
    names(dep_names) = c("y_") 
    if (dep_names == "NULL") { dep_names = NULL }
  }
  
  if (dolasso  == 1) { xbetas.cv = cbind( xbetas.cv, xbetas0.cv.lasso ) } 
  if (doxgb_   == 1) { xbetas.cv = cbind( xbetas.cv, xbetas0.cv.xgb   ) } 
  if (dorf_    == 1) { xbetas.cv = cbind( xbetas.cv, xbetas0.cv.rf    ) } 
  if (doann_   == 1) { xbetas.cv = cbind( xbetas.cv, xbetas0.cv.ann   ) } 
  if (dorpart  == 1) { xbetas.cv = cbind( xbetas.cv, xbetas0.cv.rpart ) } 
  if ((dostep == 1) | (doaic==1)) { xbetas.cv = cbind( xbetas.cv, xbetas0.cv.step ) } 
  
  rver = R.Version()$version.string

  Call <- match.call()
#  indx <- match(c("xs","start","y_","event","family"), names(Call), nomatch=0)
  
  nestedcv = list( Call = Call, 
                   sample = c (family=family, n=dim(xs)[1], nevent=nevents, xs.columns=dim(xs)[2], xs.df=rankMatrix(xs)[1], 
                               null.dev.n=null.deviance, null.m2LogLik.n=null.m2LogLik, sat.m2LogLik.n=sat.m2LogLik ), 
                   dep_names = dep_names, 
                   xvars = colnames(xs), 
                   tuning = c( folds_n=folds_n, stratified=stratified, limit=limit, fine=fine, ties=ties, method=method, steps_n=steps_n ),
                   fits   = c( dolasso=dolasso, doxgb=doxgb_, dorf=dorf_, dorpart=dorpart, doann=doann_, dostep=dostep, doaic=doaic ), 
                   seed   = seed, 
                   ensemble = ensemble, 
                   doxgb  = doxgb, 
                   dorf   = dorf, 
                   doann  = doann, 
                   do_ncv = do_ncv, 
                   foldid = foldid, 
                   y_     = y__ , 
                   xbetas = xbetas , 
                   xbetas.cv = xbetas.cv , 
                   times  = c(time_start=time_start, time_stop=time_stop),  
                   hr_mn_sc = diff_time1(time_stop, time_start),
                   version = c(R=rver, glmnetr=pver) )
#  print( keepxbetas )
  if ( keepxbetas == 0) { nestedcv$xbetas = NULL ; nestedcv$xbetas.cv = NULL ; }
  if ( do_ncv == 0) { nestedcv$xbetas.cv = NULL }
  
  if (family %in% c("cox")) { 
    names(nestedcv$sample)[c(6,7,8)] = c("null.dev/nevent","null.m2LogLik/nevent","sat.m2LogLik/nevent") 
  } else { 
    names(nestedcv$sample)[c(6,7,8)] = c("null.dev/n","null.m2LogLik/n","sat.m2LogLik/n") 
  }

  nestedcv$ null.m2LogLik.cv =  null.m2LogLik.cv
  nestedcv$sat.m2LogLik.cv   = sat.m2LogLik.cv 
  nestedcv$n.cv           = n.cv 
  
  #============================================================================#
  
  if (dolasso == 1) { 
    if (do_ncv == 1) {
      nestedcv$lasso.devian.cv    = lasso.devian.cv
      nestedcv$lasso.cal.devian.cv = lasso.cal.devian.cv
      nestedcv$lasso.intcal.cv    = lasso.intcal.cv
      nestedcv$lasso.lincal.cv    = lasso.lincal.cv
      nestedcv$lasso.agree.cv     = lasso.agree.cv
      nestedcv$lasso.nzero.cv     = lasso.nzero.cv 
    }
    nestedcv$lasso.devian.naive = lasso.devian.naive 
    nestedcv$lasso.intcal.naive = lasso.intcal.naive 
    nestedcv$lasso.lincal.naive = lasso.lincal.naive 
    nestedcv$lasso.agree.naive  = lasso.agree.naive 
    nestedcv$lasso.nzero        = lasso.nzero
    nestedcv$cv_glmnet_fit      = cv_glmnet_fit_f 
    if (folds_n >= 3) { nestedcv$cv_ridge_fit = cv_ridge_fit_f } 
  } 
  
  #-----------------------------------------------------------------------------
  
  if (doxgb_ == 1) { 
    if (do_ncv == 1) {
      nestedcv$xgb.devian.cv    = xgb.devian.cv
      nestedcv$xgb.intcal.cv    = xgb.intcal.cv 
      nestedcv$xgb.lincal.cv    = xgb.lincal.cv
      nestedcv$xgb.agree.cv     = xgb.agree.cv
      nestedcv$xgb.nzero.cv     = xgb.nzero.cv 
    }
    nestedcv$xgb.devian.naive = xgb.devian.naive
    nestedcv$xgb.intcal.naive = xgb.intcal.naive 
    nestedcv$xgb.lincal.naive = xgb.lincal.naive 
    nestedcv$xgb.agree.naive  = xgb.agree.naive 
    nestedcv$xgb.nzero        = xgb.nzero 
    if (xgbkeep == 1) {
      if (sum(ensemble[c(1,5)])     >= 1) { 
        nestedcv$xgb.simple.fit  = xgb.simple.fit
        nestedcv$xgb.tuned.fit   = xgb.tuned.fit
##        xgb.tuned.fit$param.final
##        xgb.tuned.fit$doxgb$nrounds.final
##        xgb.tuned.fit2 = xgb.train(params = xgb.tuned.fit$param.final, data = full.xgb.dat, nrounds = xgb.tuned.fit$doxgb$nrounds.final)        
      } 
      if (sum(ensemble[c(2,6)])     >= 1) {
        nestedcv$xgb.simple.fitF = xgb.simple.fitF
        nestedcv$xgb.tuned.fitF  = xgb.tuned.fitF 
      } 
      if (sum(ensemble[c(3,4,7,8)]) >= 1) {
        nestedcv$xgb.simple.fitO = xgb.simple.fitO
        nestedcv$xgb.tuned.fitO  = xgb.tuned.fitO 
      } 
    } else { 
      if (sum(ensemble[c(1,5)])     >= 1) { 
        nestedcv$doxgb_simple           = doxgb_simple 
        nestedcv$xgb.simple.param.final = xgb.simple.fit$param.final
        nestedcv$doxgb_tuned            = xgb.tuned.fit$doxgb
        nestedcv$xgb.tuned.param.final  = xgb.tuned.fit$param.final
#        nestedcv$xgb.nrounds.final = xgb.tuned.fit$doxgb$nrounds.final         
      }
      if (sum(ensemble[c(2,6)])     >= 1) {
        nestedcv$doxgb_simpleF           = doxgb_simpleF
        nestedcv$xgb.simple.param.finalF = xgb.simple.fitF$param.final
        nestedcv$doxgb_tunedF            = xgb.tuned.fitF$doxgb
        nestedcv$xgb.tuned.param.finalF  = xgb.tuned.fitF$param.final
      } 
      if (sum(ensemble[c(3,4,7,8)]) >= 1) {
        nestedcv$doxgb_simpleO           = doxgb_simpleO
        nestedcv$xgb.simple.param.finalO = xgb.simple.fitO$param.final
        nestedcv$doxgb_tunedO            = xgb.tuned.fitO$doxgb
        nestedcv$xgb.tuned.param.finalO  = xgb.tuned.fitO$param.final
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  
  if (dorf_ == 1) {
    if (do_ncv == 1) {
      nestedcv$rf.devian.cv   = rf.devian.cv
      nestedcv$rf.intcal.cv   = rf.intcal.cv
      nestedcv$rf.lincal.cv   = rf.lincal.cv
      nestedcv$rf.agree.cv    = rf.agree.cv 
      nestedcv$rf.mtry.cv     = rf.mtry.cv 
    }
    nestedcv$rf.devian.naive = rf.devian.naive
    nestedcv$rf.intcal.naive = rf.intcal.naive
    nestedcv$rf.lincal.naive = rf.lincal.naive
    nestedcv$rf.agree.naive  = rf.agree.naive
    nestedcv$rf.mtry         = rf.mtry.naive  
    if (rfkeep == 1) {
      if  (sum(ensemble[c(1,5)]) >= 1) { 
        if( keepdata == 0) { 
          rf_tuned$rf_tuned$xvar = NULL
          rf_tuned$rf_tuned$yvar = NULL
        }
        nestedcv$rf_tuned_fit  = rf_tuned 
      } 
      if  (sum(ensemble[c(2,6)]) >= 1) { 
        if ( keepdata == 0) { 
          rf_tunedF$rf_tuned$xvar = NULL
          rf_tunedF$rf_tuned$yvar = NULL
        }
        nestedcv$rf_tuned_fitF = rf_tunedF 
      }
      if ((sum(ensemble[c(3,4,7,8)]) >= 1) & (family == "gaussian")) {
        if( keepdata == 0) { 
          rf_tunedO$rf_tuned$xvar = NULL
          rf_tunedO$rf_tuned$yvar = NULL
        }
        nestedcv$rf_tuned_fitO = rf_tunedO
      } 
    } else {
      if  (sum(ensemble[c(1,5)]) >= 1) { nestedcv$dorf_base = dorf_base }
      if  (sum(ensemble[c(2,6)]) >= 1) { nestedcv$dorf_feat = dorf_feat }
      if ((sum(ensemble[c(3,4,7,8)]) >= 1) & (family == "gaussian")) { nestedcv$dorf_offs = dorf_offs }
    }
  }
  
  #-----------------------------------------------------------------------------
  
  if (doann_ == 1) {
    if (do_ncv == 1) {
      nestedcv$ann.devian.cv = ann.devian.cv
      nestedcv$ann.intcal.cv = ann.intcal.cv
      nestedcv$ann.lincal.cv = ann.lincal.cv
      nestedcv$ann.agree.cv  = ann.agree.cv
      nestedcv$ann.nzero.cv  = ann.nzero.cv
    }
    nestedcv$ann.devian.naive = ann.devian.naive
    nestedcv$ann.intcal.naive = ann.intcal.naive
    nestedcv$ann.lincal.naive = ann.lincal.naive
    nestedcv$ann.agree.naive  = ann.agree.naive 
    nestedcv$ann.nzero        = ann.nzero
    nestedcv$ann_zb      = ann_zb 
    nestedcv$ann_fit_1   = ann_fit_1_f 
    nestedcv$ann_fit_2   = ann_fit_2_f 
    nestedcv$ann_fit_3   = ann_fit_3_f 
    nestedcv$ann_fit_4   = ann_fit_4_f
    nestedcv$ann_fit_5   = ann_fit_5_f
    nestedcv$ann_fit_6   = ann_fit_6_f 
    nestedcv$ann_fit_7   = ann_fit_7_f 
    nestedcv$ann_fit_8   = ann_fit_8_f 
  } 
  
  #-----------------------------------------------------------------------------
  
  if (dorpart == 1) {
    if (do_ncv == 1) {
      nestedcv$rpart.devian.cv = rpart.devian.cv
      nestedcv$rpart.intcal.cv = rpart.intcal.cv
      nestedcv$rpart.lincal.cv = rpart.lincal.cv
      nestedcv$rpart.agree.cv  = rpart.agree.cv 
      nestedcv$rpart.nzero.cv  = rpart.nzero.cv 
    }
    nestedcv$rpart.devian.naive = rpart.devian.naive
    nestedcv$rpart.intcal.naive = rpart.intcal.naive
    nestedcv$rpart.lincal.naive = rpart.lincal.naive
    nestedcv$rpart.agree.naive  = rpart.agree.naive
    nestedcv$rpart.nzero        = rpart.nzero
    if  (sum(ensemble[c(1,5)])      >= 1) { nestedcv$rpart.fit.00  = rpart.fit.00  }
    if  (sum(ensemble[c(2,6)])      >= 1) { nestedcv$rpart.fitF.00 = rpart.fitF.00 }
    if ((sum(ensemble[c(3,4,7,8)])  >= 1) & (family %in% c("cox"))) { 
      nestedcv$rpart.fitO.00 = rpart.fitO.00 }
  } 
  
  #-----------------------------------------------------------------------------
  
  if ((dostep == 1) | (doaic==1)) { 
    if (do_ncv == 1) {
      nestedcv$step.devian.cv = step.devian.cv 
      nestedcv$step.cal.devian.cv = step.cal.devian.cv 
      nestedcv$step.intcal.cv = step.intcal.cv 
      nestedcv$step.lincal.cv = step.lincal.cv 
      nestedcv$step.agree.cv  = step.agree.cv
      nestedcv$step.nzero.cv  = step.nzero.cv 
      nestedcv$step.p.cv      = step.p.cv 
    }
  nestedcv$step.devian.naive = step.devian.naive 
  nestedcv$step.agree.naive  = step.agree.naive 
  nestedcv$step.nzero        = step.nzero 
  if (dostep == 1) { nestedcv$step.p = step.p }
  }
  
  if (dostep  == 1) { nestedcv$cv.stepreg.fit = cv.stepreg.fit.all }
  
  if (doaic == 1) { 
    if( keepdata == 0) { func.fit.aic$data = NULL }
    nestedcv$func.fit.aic = func.fit.aic 
  } 
  
  #============================================================================#
  
  class(nestedcv) <- c("nested.glmnetr")
  
  names (nestedcv)
  
  return( nestedcv )
  
  if (track >= 1) { cat("\n Program End \n") }
  if (track >= 2) { time_last = diff_time(time_start, time_last) } 
  
}

###############################################################################################################################################
###############################################################################################################################################
