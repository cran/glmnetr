################################################################################
##### predict.nested.glmnetr_yymmdd.R #########################################################
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
#' @param alpha A specific value of alpha for plotting.  Used only when type is 
#' set to "elastic". Specifies which alpha is to be used for deviance plots.
#' Default is "alpha.min", else must be an element of the alpha vector used in 
#' running the elastic net model. This can be reviewed using summary(fit) 
#' where fit is a nested.glmnetr() output object. Note, alpha is 1 for the 
#' lasso model and alpha is 0 for the ridge model.
#' @param gamma The gamma value for choice of beta.  If NULL, then 
#' gamma.min is used from the cross validation informed relaxed model.  We
#' use the term gam instead of gamma as gamma usually denotes a vector 
#' in the package.    
#' @param lambda The lambda value for choice of beta.  If NULL, then 
#' lambda.min is used from the cross validation informed relaxed model.  We
#' use the term lam instead of lambda as lambda usually denotes a vector 
#' in the package.    
#' @param type type of model on which to base predictds. One of "lasso", "ridge" 
#' and "elastic" if elastic net model is fit. 
#' @param comment Default of TRUE to write to console information on lam and gam selected for output.
#' FALSE will suppress this write to console.  
#' @param ... Additional arguments passed to the predict function.  
#'
#' @return Either the xs_new*Beta estimates based upon the predictor matrix, 
#' or model coefficients. 
#' 
#' @seealso 
#'    \code{\link{predict.cv.glmnetr}} , \code{\link{predict_ann_tab}} , \code{\link{nested.glmnetr}}
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
predict.nested.glmnetr = function( object, xs_new=NULL, alpha=NULL, gamma=NULL, lambda=NULL, type="lasso", comment=TRUE, ...) {
#  cat( "  in predict.nested.glmnetr   class(object) = ", class(object), "\n")
  retobj = predict.cv.glmnetr( object=object, xs_new=xs_new, alpha=alpha, gamma=gamma, lambda=lambda, type=type, comment=comment ) 
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
# @seealso 
#    \code{\link{predict_ann_tab}}, \code{\link{predict_nested_xgb}} , \code{\link{predict_nested_rf}} , 
#    \code{\link{predict.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
# 
# @export
#
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
