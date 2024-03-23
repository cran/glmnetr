################################################################################
##### predict_yymmdd.R #########################################################
################################################################################
#' Get predicteds for a XGB fit from a nested.glmnetr() output object
#'
#' @param object an output object from a nested.glmnetr() call 
#' @param xs an new predictor matrix of same format as used in the nested.glmnetr() call 
#' @param type NULL for the standard model, "feat" for the model including the 
#' realaxed lasso predicteds as a featuere,  "offs" for the model including the 
#' realaxed lasso predicteds as an offset.
#' @param tuned 1 (default) for tuned xgb models, 0 for the simple untuned models
#' 
#' @return predicteds from a XGB model fit in a nested.glmnetr() output object   
#' 
#' @seealso 
#'    \code{\link{predict_ann_tab}} , \code{\link{predict_nested_rf}} , \code{\link{predict.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
predict_nested_xgb = function(object, xs, type=NULL, tuned=1) {

  if (is.null(type)) { type = "base" }
  
  if ( type %in% c("feat","offs") ) {
    predminR = predict(object, xs, comment=0)
    ofst = object$lasso.intcal.naive[4] + object$lasso.lincal.naive[4] * predminR
    family = object$sample[1]
    if (family=="cox") { ofst  = ofst - mean(ofst) }
  } else {
    type = "base" 
  }
  
  if (is.null(type)) { type = "base" }
  if ( type == "base") {
    xgb.dat <- xgb.DMatrix(data = xs)
    if (tuned >= 1) {  preds = predict( object$xgb.tuned.fit , xgb.dat) 
    } else          {  preds = predict( object$xgb.simple.fit, xgb.dat)  }
  } else if (type == "feat") {
    xgb.datF <- xgb.DMatrix(data = cbind(xs,ofst))
    if (tuned >= 1) {  preds = predict( object$xgb.tuned.fitF , xgb.datF) 
    } else          {  preds = predict( object$xgb.simple.fitF, xgb.datF)  }
  } else if (type == "offs") {
    xgb.datO <- xgb.DMatrix(data = xs, base_margin=ofst)
    if (tuned >= 1) {  preds = predict( object$xgb.tuned.fitO , xgb.datO) 
    } else          {  preds = predict( object$xgb.simple.fitO, xgb.datO)  }
  }
  return(preds)
}

###############################################################################################################
###############################################################################################################
#' Get predicteds for a rf fit from nested.glmnetr() output object
#'
#' @param object an output object from a nested.glmnetr() call 
#' @param xs an new predictor matrix of same format as used in the nested.glmnetr() call 
#' @param type NULL for the standard model, "feat" for the model including the 
#' realaxed lasso predicteds as a featuere, "offs" for the model including the 
#' realaxed lasso predicteds as an offset.
#' 
#' @return predicteds from a Random Forest model fit in a nested.glmnetr() output object   
#' 
#' @seealso 
#'    \code{\link{predict_ann_tab}}, \code{\link{predict_nested_xgb}} , \code{\link{predict.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
predict_nested_rf = function(object, xs, type=NULL) {

  if (is.null(type)) { type = "base" }
  if ( type %in% c("feat","offs") ) {
    predminR = predict(object, xs, comment=0)
    ofst = object$lasso.intcal.naive[4] + object$lasso.lincal.naive[4] * predminR
  } else {
    type = "base" 
  }

  ensemble = object$ensemble 
  dorf_ = object$fits[3]
  preds = NULL 
  if (dorf_ == 1) {
    if ( (type == "base") &  (sum(ensemble[c(1,5)]) >= 1) ) {
      rf.df <- data.frame(xs)
      preds = predict( object$rf_tuned_fit$rf_tuned, rf.df )$predicted 
      
    } else if ( (type == "feat") & (sum(ensemble[c(2,6)]) >= 1) ) {
      rf.dfF <- data.frame(cbind(xs,ofst))
      preds = predict( object$rf_tuned_fitF$rf_tuned, rf.dfF )$predicted  
      
    } else if ( (type == "offs") & (sum(ensemble[c(3,4,7,8)]) >= 1) & (object$sample[1] == "gaussian") ) {
      rf.dat <- data.frame(xs)
      preds = predict( object$rf_tuned_fitO$rf_tuned, rf.dat)$predicted  
      preds = preds + ofst 
    } else if ( (type == "offs") & (sum(ensemble[c(3,4,7,8)]) >= 1) & (object$sample[1] != "gaussian") ) {
      warning(" Random Forest program does not support offset for non-gaussian errors")
    }
  } else {
    warning(" The object might not have the model stored for Random Forest")
  }
  return(preds)
}

###############################################################################################################
###############################################################################################################

#' get XBeta from an XGB.train output object
#'
#' @param xgb_model an output from xgb.train() 
#' @param xs_ and XGBoost data matrix 
#' @param family one of c("cox", "binomial", "gaussian") 
#' @param tol a small number, a lower bound to avoid division by 0
#'
#' @return the XBeta term from an XGB.train() output object 
#'
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

###############################################################################################################
#' get XBeta from an rfsrc output object
#'
#' @param rf_model an output from rfsrc() 
#' @param dframe a data fraem of new data
#' @param ofst an offset on the XB scale
#' @param family one of c("cox", "binomial", "gaussian") 
#' @param tol a small number, a lower bound to avoid division by 0
#'
#' @return the XBeta term from an rfsrc() output object
#'
rf_xbhat = function(rf_model, dframe, ofst=NULL, family, tol=1e-05) {
  if (family == "cox") {
    hrhat = predict(object=rf_model, newdata=dframe)$predicted
    hrhat[(hrhat<tol)] = tol
    hrhat[(hrhat>(1/tol))] = 1/tol 
    xbhat = log(hrhat) 
  } else if (family == "binomial") {
    phat_  = predict(object=rf_model, newdata=dframe)$predicted
    phat_t = phat_ 
    phat_t[(phat_ < tol)] = tol
    phat_t[(phat_ > (1-tol))] = 1 - tol
    xbhat  = log(phat_t/(1-phat_t))
  } else if (family == "gaussian") {
    xbhat = predict(object=rf_model, newdata=dframe)$predicted
  }
  if (!is.null(ofst)) { xbhat = xbhat + ofst } 
  return (xbhat)
}

###############################################################################################################

#' get XBeta from an rpart output object
#'
#' @param rpart_model an rpart() output object
#' @param dframe A data frame as input
#' @param ofst An offset term
#' @param family family, one of "cox", "binomial" or "gaussian" 
#' @param tol a cap at small or inverse of large numbers to avoid numeric overflow
#'
#' @return XBeta hat for processing to rpart_perf()
#'
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

###############################################################################################################



