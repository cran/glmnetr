################################################################################
##### predict_yymmdd.R #########################################################
################################################################################

#' get XBeta from an XGB.train output object
#'
#' @param xgb_model an output from xgb.train() 
#' @param xs_ and XGBoost data matrix 
#' @param family one of c("cox", "binomial", "gaussian") 
#' @param tol a small number, a lower bound to avoid division by 0
#'
#' @return XBeta hat vector
#' 
#' @noRd

xgb_xbhat = function(xgb_model, xs_, family, tol=1e-5) {
  pred = predict(xgb_model, xs_) 
  if (family == "cox") {
    pred[(pred < tol)] = tol
    pred[(pred >  (1/tol))] =  (1/tol)
    xbhat = log( pred )
  } else if (family == "binomial") {
    pred[(pred < tol)] = tol
    pred[(pred > 1 - tol)] = 1 - tol
    xbhat = log( pred / (1-pred) )
  } else if (family == "gaussian") {
    xbhat = pred 
  }
  return( xbhat ) 
}

################################################################################
#' get XBeta from an rfsrc output object
#'
#' @param rf_model an output from rfsrc() 
#' @param dframe a data fraem of new data
#' @param ofst an offset on the XB scale
#' @param family one of c("cox", "binomial", "gaussian") 
#' @param tol a small number, a lower bound to avoid division by 0
#'
#' @return XBeta hat vector
#' 
#' @noRd

rf_xbhat = function(rf_model, dframe, ofst=NULL, family, tol=1e-05, track=0) {
  if (track >= 2) { cat("     In rf_xbhat \n") }
  if (family == "cox") {
    hrhat = predict(object=rf_model, newdata=dframe)$predicted
    hrhat[(hrhat<tol)] = tol
    hrhat[(hrhat>(1/tol))] = 1/tol 
    xbhat = log(hrhat) 
  } else if (family == "binomial") {
    phat_  = predict(object=rf_model, newdata=dframe)$predicted
    if (length(dim(phat_)) == 2) { phat_ = phat_[,2] } 
    phat_t = phat_ 
    phat_t[(phat_ < tol)] = tol
    phat_t[(phat_ > (1-tol))] = 1 - tol
    xbhat  = log(phat_t/(1-phat_t))
  } else if (family == "gaussian") {
    if (track >= 2) { cat("     In rf_xbhat  gaussian \n") }
    xbhat = predict(object=rf_model, newdata=dframe)$predicted
    if (track >= 2) {
      print(summary( xbhat ))
      cat("     In rf_xbhat  gaussian \n" )
      cat("  class(xbhat) = ", class(xbhat), "\n")
      cat("  length(xbhat) = ", length(xbhat), "\n")
      cat("  class(ofst) = " , class(ofst), "\n")
      cat("  length(ofst) = ", length(ofst), "\n")
    }
  }
  
  if (!is.null(ofst)) { xbhat = xbhat + ofst } 
  return (xbhat)
}

################################################################################
#' get XBeta from an rfsrc output object
#'
#' @param orf_model an output from rfsrc() 
#' @param dframe a data fraem of new data
#' @param ofst an offset on the XB scale
#' @param family one of c("cox", "binomial", "gaussian") 
#' @param tol a small number, a lower bound to avoid division by 0
#'
#' @return XBeta hat vector
#' 
#' @noRd

orf_xbhat = function(orf_model, dframe, ofst=NULL, family, tol=1e-05) {
  if (family == "cox") {
    predorf = predict(orf_model, dframe, pred_type='chf')
    hrhat = predorf / mean(predorf)
    hrhat[(hrhat<tol)] = tol
    hrhat[(hrhat>(1/tol))] = 1/tol 
    xbhat = log(hrhat) 
  } else if (family == "binomial") {
    phat_  = predict(object=orf_model, dframe , pred_type='prob') 
    if (length(dim(phat_)) == 2) { phat_ = phat_[,2] } 
#    summary(phat_)
    phat_t = phat_ 
    phat_t[(phat_ < tol)] = tol
    phat_t[(phat_ > (1-tol))] = 1 - tol
    xbhat  = log(phat_t/(1-phat_t))
#    xbhat  = phat_ 
  } else if (family == "gaussian") {
    xbhat = predict(object=orf_model, dframe) 
  }
  if (!is.null(ofst)) { xbhat = xbhat + ofst } 
  return (xbhat)
}

################################################################################

#' get XBeta from an rpart output object
#'
#' @param rpart_model an rpart() output object
#' @param dframe A data frame as input
#' @param ofst An offset term
#' @param family family, one of "cox", "binomial" or "gaussian" 
#' @param tol a cap at small or inverse of large numbers to avoid numeric overflow
#'
#' @return XBeta hat vector 
#' 
#' @noRd

rpart_xbhat = function(rpart_model, dframe, ofst=NULL, family, tol=1e-5) {
  if (family == "cox") {
    hrhat  = predict(rpart_model, dframe, type="vector")  
    hrhat[(hrhat<tol)] = tol
    hrhat[(hrhat>(1/tol))] = 1/tol 
    xbhat = log(hrhat) 
  } else if (family == "binomial") {
    phat_  = predict(rpart_model, dframe, type="prob")[,2]  
    phat_t = phat_ 
    phat_t[(phat_ < tol)] = tol
    phat_t[(phat_ > (1-tol))] = 1 - tol
    xbhat  = log(phat_t/(1-phat_t))
  } else if (family == "gaussian") {
    xbhat = predict(rpart_model, dframe, type="vector")                         ## predict.rpart disregards values of offset 
  }
  if (!is.null(ofst)) { xbhat = xbhat + ofst }
  return( xbhat )
}

################################################################################
#' get XBeta from an ann_tab_cv() or ann_tab_cv_best() output object
#'
#' @param object an ann_tab_cv() or ann_tab_cv_best() output object
#' @param newdat A new feature data matrix as input
#' @param family family, one of "cox", "binomial" or "gaussian" 
#' @param tol a cap at small or inverse of large numbers to avoid numeric overflow
#'
#' @return XBeta hat vector
#' 
#' @noRd

ann_xbhat = function(object, newdat, family="binomial", tol=1e-5) {
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

################################################################################

#' Title
#'
#' @param pred     predicteds for calibration
#' @param trainy__ training y
#' @param pred.tr  predicteds from the training data  
#' @param family   model family, i.e. of (c("cox","binomial","gaussian"))
#'
#' @return S Betas 
#' 
#' @noRd
#' 
cal_train_xbhat = function( pred, trainy__, pred.tr, family ) { 
  perf.tr  = perf_gen( trainy__ , pred.tr  , family )
  pred.cal = perf.tr[4] + pred * perf.tr[5]
  return( pred.cal )
}

################################################################################


