###############################################################################################################

#' Calculate performance metrics for binomial family
#'
#' @param yy outcome (0 or 1) variable
#' @param xbhat predicted values
#'
#' @return a vector of performance metrics
#' 
#' @noRd

perf_bin = function(yy, xbhat) {
  fit0 = glm( yy ~ xbhat , family="binomial") 
  p_ = 1/(1+exp(-xbhat)) 
  retvec = c( dev1 = -2*( t(log(p_))%*%yy + t(log(1-p_))%*%(1-yy) ) / length(yy) , 
              dev2 = fit0$deviance / length(yy) ,
              agree  = concordance(fit0)[[1]] ,          
              intcal = fit0$coefficients[1] , 
              lincal = fit0$coefficients[2] ) 
  return( retvec )
}

################################################################################

#' Calculate performance metrics for cox family
#'
#' @param SURV a Surv object as outcome variable
#' @param xbhat predicted values
#'
#' @return a vector of performance metrics
#' 
#' @noRd

perf_cox = function(SURV, xbhat) {
  fit0 = coxph( SURV ~ xbhat, init=c(1)) 
  retvec = c(dev1 = -2*fit0$loglik[1] / fit0$nevent ,
             dev2 = -2*fit0$loglik[2] / fit0$nevent , 
             agree = fit0$concordance[6] ,
             intcal = 0 , 
             lincal = fit0$coefficients ) 
  return( retvec )
}

################################################################################

#' Calculate performance metrics for guassian family
#'
#' @param yy outcome variable
#' @param xbhat predicted values
#'
#' @return a vector of performance metrics
#' 
#' @noRd

perf_gau = function(yy, xbhat) {
  fit0 = glm( yy ~ xbhat , family="gaussian")
  returnvec = c(  dev1 = sum((yy - xbhat)^2) / length(yy) ,
                  dev2 = fit0$deviance / length(yy) ,
                  agree = ifelse( var(xbhat)>0 , cor(x=yy, y=xbhat) , 0 ) ,
                  intcal = fit0$coefficients[1] , 
                  lincal = fit0$coefficients[2] )
  return( returnvec )
}

################################################################################

#' Calculate performance metrics for binomial, cox or gaussian family
#'
#' @param yy outcome variable 
#' @param xbhat predicted values
#' @param family one of "binomial", "cox" or "gaussian" 
#'
#' @return a vector of performance metrics
#' 
#' @noRd

perf_gen = function(yy, xbhat, family) {
  if        (family == "cox")      { returnvec = perf_cox(yy, xbhat) 
  } else if (family == "binomial") { returnvec = perf_bin(yy, xbhat) 
  } else if (family == "gaussian") { returnvec = perf_gau(yy, xbhat) } 
  return(returnvec)
}

###############################################################################################################
#' Get model performance metrics for XGB model
#'
#' @param xgb_model a model from xgb_simple() or xgb.tuned() 
#' @param dframe a data frame for new data
#' @param ofst an offset term 
#' @param y__ the outcome variable
#' @param family one of "cox", "binomial" or "gaussian" 
#' @param tol a small tolerance value 
#'
#' @return a numeric (vector) for model performance
#' 
#' @noRd

xgb_perform = function(xgb_model, xs_, y, family, tol=1e-5) {
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

###############################################################################################################
#' Get model performance metrics for Random Forest SRC model
#'
#' @param rf_model a model from rf_tune() 
#' @param dframe a data frame for new data
#' @param ofst an offset term which only works for "gaussian" 
#' @param y__ the outcome variable
#' @param family one of "cox", "binomial" or "gaussian" 
#' @param tol a small tolerance value 
#'
#' @return a numeric (vector) for model performance
#' 
#' @noRd

rf_perform = function(rf_model, dframe, ofst=NULL, y__, family, tol=1e-5) {
  xbhat = rf_xbhat(rf_model, dframe, ofst, family, tol) 
  summary(xbhat)
  if (family == "cox") {
    returnvec = perf_cox(y__, xbhat)
  } else if (family == "binomial") {
    returnvec = perf_bin(y__, xbhat) 
  } else if (family == "gaussian") {
    returnvec = perf_gau(y__, xbhat ) 
  }
  returnvec = c(returnvec, rf_model$mtry)
}

###############################################################################################################
#' Get model performance metrics for Oblique Random Forest model
#'
#' @param orf_model a model from orf_tune() 
#' @param dframe a data frame for new data
#' @param ofst an offset term which only works for "gaussian" 
#' @param y__ the outcome variable
#' @param family one of "cox", "binomial" or "gaussian" 
#' @param tol a small tolerance value 
#'
#' @return a numeric (vector) for model performance
#' 
#' @noRd

orf_perform = function(orf_model, dframe, ofst=NULL, y__, family, tol=1e-5) {
  xbhat = orf_xbhat(orf_model, dframe, ofst, family, tol) 
  if (family == "cox") {
    returnvec = perf_cox(y__, xbhat)
  } else if (family == "binomial") {
#    print(summary(xbhat))
    if (is.factor( y__ ) ) { y__ = as.numeric(y__) - 1 }  
#    print( table(y__) ) 
    returnvec = perf_bin(y__, xbhat) 
  } else if (family == "gaussian") {
    returnvec = perf_gau(y__, xbhat ) 
  }
  returnvec = c(returnvec, orf_model$mtry)
}

###############################################################################################################

#' Get number of leaves in an RPART model 
#'
#' @param obj a RPART output object
#'
#' @return number of leaves
#' 
#' @noRd

rpnz = function(obj) { 
  frame = obj$frame ;  leaves = frame$var == "<leaf>" ;  used <- unique(frame$var[!leaves]) ; 
  return( length(used) )
}

###############################################################################################################


#' Get model performance metrics for RPART
#'
#' @param rpart_model a model from rpart() 
#' @param dframe a data frame for new data
#' @param ofst an offset term which only works for "binomial" 
#' @param y__ the outcome variable
#' @param family one of "cox", "binomial" or "gaussian" 
#' @param tol a small tolerance value 
#'
#' @return a numeric (vector) for model performance
#' 
#' @noRd

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

###############################################################################################################
#' Get model performance metrics for ANN
#'
#' @param object a model from ann_tab_cv() or ann_tab_cv_best()
#' @param newdat a data matrix for new feature data
#' @param newy the outcome variable
#' @param family one of "cox", "binomial" or "gaussian" 
#' @param start start time for a Cox model framework, not yet implemented
#' @param event event variable for survival data
#' @param tol a small tolerance value 
#'
#' @return a numeric (vector) for model performance
#' 
#' @noRd

ann_perform = function(object, newdat, newy, family="binomial", start=NULL, event=NULL, tol=1e-5) {
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

###############################################################################################################



