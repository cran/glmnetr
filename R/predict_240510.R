################################################################################
##### predict_yymmdd.R #########################################################
################################################################################

#' Give predicteds based upon the cv.glmnet output object contained in the 
#' nested.glmnetr output object.
#'
#' @description This is essentially a redirect to the summary.cv.glmnetr
#' function for nested.glmnetr output objects, based uopn the cv.glmnetr
#' output object contained in the nested.glmnetr output object.  
#' 
#' @param object  A nested.glmnetr output object.
#' @param xs_new The predictor matrix.  If NULL, then betas are provided.  
#' @param lam The lambda value for choice of beta.  If NULL, then 
#' lambda.min is used from the cross validation informed relaxed model.  We
#' use the term lam instead of lambda as lambda usually denotes a vector 
#' in the package.    
#' @param gam The gamma value for choice of beta.  If NULL, then 
#' gamma.min is used from the cross validation informed relaxed model.  We
#' use the term gam instead of gamma as gamma usually denotes a vector 
#' in the package.    
#' @param comment Default of TRUE to write to console information on lam and gam selected for output.
#' FALSE will suppress this write to console.  
#' @param ... Additional arguments passed to the predict function.  
#'
#' @return Either the xs_new*Beta estimates based upon the predictor matrix, 
#' or model coefficients. 
#' 
#' @seealso 
#'    \code{\link{predict.cv.glmnetr}} , \code{\link{predict_ann_tab}} , \code{\link{nested.glmnetr}}
#    , \code{\link{predict_nested_xgb}} , \code{\link{predict_nested_rf}} 
#' 
#' @export 
#'
#' @examples
#' \donttest{
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$yt
#' event=sim.data$event
#' # for this example we use a small number for folds_n to shorten run time 
#' fit3 = nested.glmnetr(xs, NULL, y_, event, family="cox", folds_n=3) 
#' betas = predict(fit3)
#' betas$beta
#' }
#' 
predict.nested.glmnetr = function( object, xs_new=NULL, lam=NULL, gam=NULL, comment=TRUE, ...) {
  if (class(object) %in% c("nested.glmnetr","nested.stepreg")) { object = object$cv_glmnet_fit }
  retobj = predict.cv.glmnetr( object=object, xs_new=xs_new, lam=lam, gam=gam, comment=comment) 
  return(retobj)
}

################################################################################
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
#' @return  XBeta hat vector
#' 
# @seealso 
#    \code{\link{predict_ann_tab}} , \code{\link{predict_nested_rf}} , \code{\link{predict_nested_orf}} , 
#    \code{\link{predict.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
# 
# @export
#
#' @noRd

predict_nested_xgb = function(object, xs, type=NULL, tuned=1) {
  
  family = object$sample[1]

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
#' @return  XBeta hat vector
#' 
# @seealso 
#    \code{\link{predict_ann_tab}}, \code{\link{predict_nested_xgb}} , \code{\link{predict_nested_orf}} , 
#    \code{\link{predict.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
# 
# @export
#
#' @noRd

predict_nested_rf = function(object, xs, type=NULL) {
  
  family = object$sample[1]

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
      preds = rf_xbhat(object$rf_tuned_fit$rf_tuned, rf.df, ofst=NULL, family, tol=1e-05) 
    } else if ( (type == "feat") & (sum(ensemble[c(2,6)]) >= 1) ) {
      rf.dfF <- data.frame(cbind(xs,ofst))
      preds = rf_xbhat(object$rf_tuned_fitF$rf_tuned, rf.dfF, ofst=NULL, family, tol=1e-05) 
    } else if ( (type == "offs") & (sum(ensemble[c(3,4,7,8)]) >= 1) & (object$sample[1] == "gaussian") ) {
      rf.dat <- data.frame(xs)
      preds = rf_xbhat(object$rf_tuned_fitO$rf_tuned, rf.dat, ofst=ofst, family, tol=1e-05) 
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
#' Get predicteds for an orf fit from nested.glmnetr() output object
#'
#' @param object an output object from a nested.glmnetr() call 
#' @param xs an new predictor matrix of same format as used in the nested.glmnetr() call 
#' @param type NULL for the standard model, "feat" for the model including the 
#' realaxed lasso predicteds as a featuere, "offs" for the model including the 
#' realaxed lasso predicteds as an offset.
#' 
#' @return  XBeta hat vector
#' 
# @seealso 
#    \code{\link{predict_ann_tab}}, \code{\link{predict_nested_xgb}} , \code{\link{predict_nested_rf}} , 
#    \code{\link{predict.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
# 
# @export
#
#' @noRd

predict_nested_orf = function(object, xs, type=NULL) {
  
  family = object$sample[1]
  
  if (is.null(type)) { type = "base" }
  if ( type %in% c("feat","offs") ) {
    predminR = predict(object, xs, comment=0)
    ofst = object$lasso.intcal.naive[4] + object$lasso.lincal.naive[4] * predminR
  } else {
    type = "base" 
  }
  
  ensemble = object$ensemble 
  doorf_ = object$fits[8]
  preds = NULL 
  if (doorf_ == 1) {
    if ( (type == "base") &  (sum(ensemble[c(1,5)]) >= 1) ) {
      orf.df <- data.frame(xs)
      preds = orf_xbhat(object$orf_tuned_fit$orf_tuned, orf.df, ofst=NULL, family, tol=1e-05) 
    } else if ( (type == "feat") & (sum(ensemble[c(2,6)]) >= 1) ) {
      orf.dfF <- data.frame(cbind(xs,ofst))
      preds = orf_xbhat(object$orf_tuned_fitF$orf_tuned, orf.dfF, ofst=NULL, family, tol=1e-05) 
    } else if ( (type == "offs") & (sum(ensemble[c(3,4,7,8)]) >= 1) & (object$sample[1] == "gaussian") ) {
      orf.dat <- data.frame(xs)
      preds = orf_xbhat(object$orf_tuned_fitO$orf_tuned, orf.dat, ofst=ofst, family, tol=1e-05) 
    } else if ( (type == "offs") & (sum(ensemble[c(3,4,7,8)]) >= 1) & (object$sample[1] != "gaussian") ) {
      warning(" Oblique Random Forest program does not support offset for non-gaussian errors")
    }
  } else {
    warning(" The object might not have the model stored for Oblique Random Forest")
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
#' @return XBeta hat vector
#' 
# @export
#
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

###############################################################################################################
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
# @export
#
#' @noRd

rf_xbhat = function(rf_model, dframe, ofst=NULL, family, tol=1e-05) {
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
    xbhat = predict(object=rf_model, newdata=dframe)$predicted
  }
  if (!is.null(ofst)) { xbhat = xbhat + ofst } 
  return (xbhat)
}

###############################################################################################################
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
# @export
#
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

###############################################################################################################

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
# @export
#
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

###############################################################################################################
#' get XBeta from an ann_tab_cv() or ann_tab_cv_best() output object
#'
#' @param object an ann_tab_cv() or ann_tab_cv_best() output object
#' @param newdat A new feature data matrix as input
#' @param family family, one of "cox", "binomial" or "gaussian" 
#' @param tol a cap at small or inverse of large numbers to avoid numeric overflow
#'
#' @return XBeta hat vector
#' 
# @export
#
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

###############################################################################################################



