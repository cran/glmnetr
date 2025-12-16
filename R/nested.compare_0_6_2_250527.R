###############################################################################################################
###############################################################################################################
#' Calculate performance differences with CI and p
#' 
#' @description 
#' Perform a paired t-test as called from nested.compare().  
#'
#' @param a One term
#' @param b A second term 
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4
#' @param txt 1 (default) to include inline text for estimated, 95 percent CI and p
#' @param pow Power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  Pow is ignored for the family of "cox" and "binomial".
#' @param bootstrap 1 if a bootstrap fit was done as opposed to a (neted) CV. 
#' 
#' @return An estimate, 95% CI and p for agreement comparison 
#' 
#' @importFrom stats t.test qt pt var 
#'
#' @noRd
#' 
nested.compare0_0_6_2 = function(df, group, desc, a, b, pow=1, bootstrap=0, digits=4, txt=0 ) { 
  if ( pow != 2) { pow = 1 } 
  if (pow == 1) {
    dff = (a-b)[!is.na(a-b)]
    meandiff = mean(dff) 
    n_ = length(dff)
    if (bootstrap == 0) {
      se = sd(dff)/sqrt(n_) 
    } else if (bootstrap >= 1) {
      se = sd(dff) 
    }
    lo = meandiff - qt(0.975,df=(n_-1))*se
    up = meandiff + qt(0.975,df=(n_-1))*se
    t_ = meandiff/se
    p_ = 2*pt(-abs(t_), df=(n_-1)) 
  } else if ( pow == 2) {
    n_ = length(a)
    deltalo1 = rep(0,n_)
    for ( i_ in c(1:n_)) {
      deltalo1[i_] = mean(a[-i_])^2 - mean(b[-i_])^2 
    }
    deltasd = sqrt( (n_+1) * var(deltalo1) ) 
    corr1 = mean(a) 
    corr2 = mean(b) 
    mean =  mean(a)^2 - mean(b)^2 
    qt_ = qt(0.975,(n_-1))
    lo = mean - qt_ * deltasd 
    up = mean + qt_ * deltasd 
    t_ = mean / deltasd  
    p_ = 2*pt(-abs(t_), df=(n_-1)) 
  }

  dfr = as.data.frame(list(group=group, desc=desc, meandiff=round(meandiff, digits=digits), lower95=round(lo, digits=digits), upper95=round(up, digits=digits), p=round(p_, digits=digits)))
  if (is.null(df)) { df = dfr 
  } else { df = rbind( df, dfr) }
  return( df )
}

###############################################################################################################
###############################################################################################################

#' Compare cross validation fit performances from a nested.glmnetr output.
#'
#' @description 
#' Compare cross-validation model fits in terms of average performances from the 
#' nested cross validation fits. In general the standard deviations for the
#' performance measures evaluated on the leave-out samples may be biased. While 
#' the standard deviations of the paired within fold differences of 
#' performances intuitively might be less biased this has not been shown. See 
#' the package vignettes for more discussion.    
#'
#' @param object A nested.glmnetr output object.
#' @param type determines what type of nested cross validation performance measures are 
#' compared.  Possible values are "devrat" to compare the deviance ratios, i.e. the 
#' fractional reduction in deviance relative to the null model deviance, 
#' "agree" to compare agreement, "lincal" to compare the linear calibration 
#' slope coefficients, "intcal" to compare the linear calibration intercept 
#' coefficients, from the nested cross validation. 
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4 
#' @param pow the power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  pow is ignored for the family of "cox" and "binomial".  
#' @param table 1 to print table to console, 0 to output the tabled information to a data frame
#'  
#' @return A printout to the R console. 
#' 
#' @noRd
#' 
nested.compare_0_6_2 = function( object, type="devrat", digits=4, pow=1, table=1 ) {
  family =  object$sample[1]
  tuning  = object$tuning
  bootstrap = tuning[8] 
  fits    = object$fits 
  dolasso = fits[1]
  doxgb   = fits[2]
  dorf    = fits[3]
  dorpart = fits[4]
  doann   = fits[5]
  dostep  = fits[6]
  doaic   = fits[7]
  doorf   = fits[8]
  dofull  = fits[9]
#  if (is.na(doorf)) { doorf = 0 }
  
  ensemble  = object$ensemble 
  
  if (!(type %in% c('devrat','devian','intcal','lincal','agree'))) { type = 'devrat'}
  
  m2.ll.null = object$null.m2LogLik.rep
  m2.ll.sat  = object$sat.m2LogLik.rep
  n__ = object$n.rep
  
  if (dolasso == 1) {
    if (type == "agree") {
      lasso.perf.rep  = object$lasso.agree.rep
    } else if (type == "devrat") {
      object1  = object$lasso.devian.rep
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      lasso.perf.rep = object1
    } else if (type == "devian") { lasso.perf.rep = object$lasso.devian.rep
    } else if (type == "lincal") { lasso.perf.rep = object$lasso.lincal.rep
    } else if (type == "intcal") { lasso.perf.rep = object$lasso.intcal.rep }
  }
  
  if (doxgb == 1) {
    if (type == "agree") {
      xgb.perf.rep  = object$xgb.agree.rep
    } else if (type == "devrat") {
      object1  = object$xgb.devian.rep
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      xgb.perf.rep = object1
    } else if (type == "devian") { xgb.perf.rep  = object$xgb.devian.rep
    } else if (type == "lincal") { xgb.perf.rep  = object$xgb.lincal.rep
    } else if (type == "intcal") { xgb.perf.rep  = object$xgb.intcal.rep }
  }
  
  if (dorf == 1) {
    if (type == "agree") {
      rf.perf.rep  = object$rf.agree.rep
    } else if (type == "devrat") {
      object1  = object$rf.devian.rep
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      rf.perf.rep = object1
    } else if (type == "devian") { rf.perf.rep = object$rf.devian.rep
    } else if (type == "lincal") { rf.perf.rep = object$rf.lincal.rep
    } else if (type == "intcal") { rf.perf.rep = object$rf.intcal.rep }
  }
  
  if (doorf == 1) {
    if (type == "agree") {
      orf.perf.rep  = object$orf.agree.rep
    } else if (type == "devrat") {
      object1  = object$orf.devian.rep
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      orf.perf.rep = object1
    } else if (type == "devian") { orf.perf.rep  = object$orf.devian.rep
    } else if (type == "lincal") { orf.perf.rep  = object$orf.lincal.rep
    } else if (type == "intcal") { orf.perf.rep  = object$orf.intcal.rep }
  }
  
  if (doann == 1) {
    if (type == "agree") {
      ann.perf.rep  = object$ann.agree.rep
    } else if (type == "devrat") {
      object1  = object$ann.devian.rep
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      ann.perf.rep = object1
    } else if (type == "devian") { ann.perf.rep  = object$ann.devian.rep
    } else if (type == "lincal") { ann.perf.rep  = object$ann.lincal.rep
    } else if (type == "intcal") { ann.perf.rep  = object$ann.intcal.rep }
  }
  
  if (dorpart == 1) {
    if (type == "agree") {
      rpart.perf.rep  = object$rpart.agree.rep
    } else if (type == "devrat") {
      object1  = object$rpart.devian.rep
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      rpart.perf.rep = object1
    } else if (type == "devian") { rpart.perf.rep  = object$rpart.devian.rep 
    } else if (type == "lincal") { rpart.perf.rep  = object$rpart.lincal.rep 
    } else if (type == "intcal") { rpart.perf.rep  = object$rpart.intcal.rep }
  }
  
  if ((dostep == 1) | (dofull == 1))  {
    if (type == "agree") {
      step.perf.rep  = object$step.agree.rep
    } else if (type == "devrat") {
      object1  = object$step.devian.rep
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      step.perf.rep = object1
    } else if (type == "devian") { step.perf.rep  = object$step.devian.rep
    } else if (type == "lincal") { step.perf.rep  = object$step.lincal.rep
    } else if (type == "intcal") { step.perf.rep  = object$step.intcal.rep }
  }   

  #  pow = 1 
  
  if (type == "agree") {
    if (family == "gaussian") { 
      if (pow == 2) { 
        pm = "** R-square **"
      } else { 
        pow = 1 
        pm = "** Correlation **" 
      }
    } else if (family %in% c("cox","binomial")) { 
      pow = 1 
      pm = "** Concordance **" 
    }
  } else if (type=="devrat") {
    pm = "** Deviance Ratio **"
  } else if (type=="devian") {
    pm = "** Deviance **"
  } else if (type=="lincal") {
    pm = "** linear calibration slope **"
  } else if (type=="intcal") {
    pm = "** linear calibration intercept **"
  }
  
  #  if ( pow == 2 ) { gagree = "R-square" } else { gagree = "Correlation" }
  #  if (family %in% c("cox","binomial")) { gagree = "Concordance" ; pow = 1 }
  
  if (sum(ensemble[c(2:8)]) > 0 ) {
    cat ("  Ensemble option used when fitting models : ")
    # ensemble\n" ) 
    cat(paste0("(", ensemble[1],",", ensemble[2],",", ensemble[3],",", ensemble[4],", ",
               ensemble[5],",", ensemble[6],",", ensemble[7],",", ensemble[8],")\n\n")) 
  } 
  
  if ( ensemble[c(1)] == 0 ) {
    cat ("\n Simple models with information from losso not run.  Output is abbreviated. \n\n" ) 
    doxgb = 0 ; dorf = 0 ; doorf = 0 ; doann = 0 ; dorpart = 0 ; dostep = 0 ; doaic = 0 ; dofull = 0 ; 
  }
  
  cat ("  Model performance comparisons in terms of", pm, "\n\n" ) 
  
  df = NULL 
  
  if (dolasso == 1) {
    group = "01:lasso"
    df = nested.compare0_0_6_2(df, group, "lasso relaxed - lasso"           , lasso.perf.rep[,2] , lasso.perf.rep[,1], pow, bootstrap, digits) 
    df = nested.compare0_0_6_2(df, group, "lasso relaxed - lasso relaxed G0", lasso.perf.rep[,2] , lasso.perf.rep[,3], pow, bootstrap, digits)   
    if (!is.null(object$cv_elastic_fit)) {
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - elastic net"   , lasso.perf.rep[,2] , lasso.perf.rep[,4], pow, bootstrap, digits)   
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - elastic net G0", lasso.perf.rep[,2] , lasso.perf.rep[,6], pow, bootstrap, digits)   
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - elastic net G1", lasso.perf.rep[,2] , lasso.perf.rep[,7], pow, bootstrap, digits)   
      df = nested.compare0_0_6_2(df, group, "ridge - elastic net G1"         , lasso.perf.rep[,5] , lasso.perf.rep[,7], pow, bootstrap, digits)   
    }
    df = nested.compare0_0_6_2(df, group, "lasso relaxed - ridge"    , lasso.perf.rep[,2] , lasso.perf.rep[,5], pow, bootstrap, digits)   
    df = nested.compare0_0_6_2(df, group, "lasso - lasso relaxed G0"  , lasso.perf.rep[,1] , lasso.perf.rep[,3], pow, bootstrap, digits)   
  }
  
  if (doxgb == 1) {
    group = "02:xgb"
    df = nested.compare0_0_6_2(df, group, "XGBoost (tuned) - XGBoost (simple)", xgb.perf.rep[,4] , xgb.perf.rep[,1], pow, bootstrap, digits)  
    if (sum(ensemble[c(2,6)])> 0) {
      df = nested.compare0_0_6_2(df, group, "XGBoost (tuned) lasso feature - no feature", xgb.perf.rep[,5] , xgb.perf.rep[,4], pow, bootstrap, digits)   
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      df = nested.compare0_0_6_2(df, group, "XGBoost (tuned) lasso offset - no offset"  , xgb.perf.rep[,6] , xgb.perf.rep[,4], pow, bootstrap, digits)   
    }
  }
  
  if (dorf == 1) {
    group = "03:RF"
    if (sum(ensemble[c(2,6)])> 0) {
      df = nested.compare0_0_6_2(df, group, "RF with lasso feature - no feature", rf.perf.rep[,2] , rf.perf.rep[,1], pow, bootstrap, digits)   
    }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) {
      df = nested.compare0_0_6_2(df, group, "RF with lasso offset - no offset", rf.perf.rep[,3] , rf.perf.rep[,1], pow, bootstrap, digits)   
    }
  }
  
  if (doorf == 1) {
    group = "04:ORF"
    if (sum(ensemble[c(2,6)])> 0) {
      df = nested.compare0_0_6_2(df, group, "ORF with lasso feature - no feature", orf.perf.rep[,2] , orf.perf.rep[,1], pow, bootstrap, digits)   
    }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) {
      df = nested.compare0_0_6_2(df, group, "ORF with lasso offset - no offset"  , orf.perf.rep[,3] , orf.perf.rep[,1], pow, bootstrap, digits) 
    }
  }
  
  if (doann == 1) {
    group = "05:ANN"
    if (sum(ensemble[6])> 0) {
      df = nested.compare0_0_6_2(df, group, "ANN with with lasso feature - no feature", ann.perf.rep[,6] , ann.perf.rep[,1], pow, bootstrap, digits)  
    } else if (sum(ensemble[2])> 0) {
      df = nested.compare0_0_6_2(df, group, "ANN with with lasso feature - no feature", ann.perf.rep[,2] , ann.perf.rep[,1], pow, bootstrap, digits) 
    } 
    if (sum(ensemble[8])> 0) { 
      df = nested.compare0_0_6_2(df, group, "ANN with with lasso offset - no offset", ann.perf.rep[,8] , ann.perf.rep[,1], pow, bootstrap, digits)  
    } else if (sum(ensemble[7])> 0) { 
      df = nested.compare0_0_6_2(df, group, "ANN with with lasso offset - no offset", ann.perf.rep[,7] , ann.perf.rep[,1], pow, bootstrap, digits)  
    } else if (sum(ensemble[4])> 0) { 
      df = nested.compare0_0_6_2(df, group, "ANN with with lasso offset - no offset", ann.perf.rep[,4] , ann.perf.rep[,1], pow, bootstrap, digits) 
    } else if (sum(ensemble[3])> 0) { 
      df = nested.compare0_0_6_2(df, group, "ANN with with lasso offset - no offset", ann.perf.rep[,3] , ann.perf.rep[,1], pow, bootstrap, digits) 
    } 
  }
  
  if (dostep == 1) {
    group = "06:step"
    df = nested.compare0_0_6_2(df, group, "step (df) - step (p)", step.perf.rep[,1]      , step.perf.rep[,2], pow, bootstrap, digits) ;  cat("\n")
  }

  if ((dolasso == 1) & (dostep == 1)) {
    group = "07:lasso-step"
    df = nested.compare0_0_6_2(df, group, "lasso relaxed G0 - step (p)", lasso.perf.rep[,3] , step.perf.rep[,2], pow, bootstrap, digits)   
  }

  if ((dolasso == 1) & (dofull == 1)) {
    group = "08:lasso-full"
    df = nested.compare0_0_6_2(df, group, "lasso relaxed - full model", lasso.perf.rep[,2] , step.perf.rep[,4], pow, bootstrap, digits) 
  }
  
    
  if ((dolasso == 1) & (doxgb == 1)) {
    group = "08:lasso-full"
    df = nested.compare0_0_6_2(df, group, "lasso relaxed - XGB (tuned)", lasso.perf.rep[,2] , xgb.perf.rep[,4], pow, bootstrap, digits) 
    if (sum(ensemble[c(2,6)])> 0) {
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - XGB with lasso feature", lasso.perf.rep[,2] , xgb.perf.rep[,5], pow, bootstrap, digits) 
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      df =nested.compare0_0_6_2(df, group, "lasso relaxed - XGB with lasso offset", lasso.perf.rep[,2] , xgb.perf.rep[,6], pow, bootstrap, digits)   
    }
  }
  
  if ((dolasso == 1) & (dorf == 1)) {
    group = "09:lasso-RF"
    df = nested.compare0_0_6_2(df, group, "lasso relaxed - Random Forest", lasso.perf.rep[,2] , rf.perf.rep[,1], pow, bootstrap, digits) 
    if (sum(ensemble[c(2,6)])> 0) {
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - RF with lasso feature", lasso.perf.rep[,2] , rf.perf.rep[,2], pow, bootstrap, digits) 
    }
    if ( (sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian") ) {
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - RF with lasso offset", lasso.perf.rep[,2] , rf.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if ((dolasso == 1) & (doorf == 1)) {
    group = "10:lasso-ORF"
    df = nested.compare0_0_6_2(df, group, "lasso relaxed - Oblique Random Forest", lasso.perf.rep[,2] , orf.perf.rep[,1], pow, bootstrap, digits) 
    if (sum(ensemble[c(2,6)])> 0) {
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - ORF with lasso feature", lasso.perf.rep[,2] , orf.perf.rep[,2], pow, bootstrap, digits) 
    }
    if ( (sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian") ) {
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - ORF with lasso offset", lasso.perf.rep[,2] , orf.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if ((dolasso == 1) & (doann == 1)) {
    group = "11:lasso-ANN"
    df = nested.compare0_0_6_2(df, group, "lasso relaxed - ANN", lasso.perf.rep[,2] , ann.perf.rep[,1], pow, bootstrap, digits) 
    if (ensemble[6]) { 
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - ANN l lasso feature", lasso.perf.rep[,2] , ann.perf.rep[,6], pow, bootstrap, digits)   
    } else if (ensemble[2]) { 
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - ANN lasso feature"  , lasso.perf.rep[,2] ,  ann.perf.rep[,2], pow, bootstrap, digits)   
    }  
    if (ensemble[8]) { 
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - ANN l lasso offset (upated)", lasso.perf.rep[,2] , ann.perf.rep[,8], pow, bootstrap, digits)   
    } else if (ensemble[4]) { 
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - ANN lasso offset"           , lasso.perf.rep[,2] , ann.perf.rep[,4], pow, bootstrap, digits)  
    } else if (ensemble[7]) { 
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - ANN l lasso offset (upated)", lasso.perf.rep[,2] , ann.perf.rep[,7], pow, bootstrap, digits)   
    }  else if (ensemble[3]) { 
      df = nested.compare0_0_6_2(df, group, "lasso relaxed - ANN lasso offset"           , lasso.perf.rep[,2] , ann.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if ((doxgb == 1) & (dorf == 1)) {
    group = "12:XGB-RF"
    df = nested.compare0_0_6_2(df, group, "XGBoost (tuned) - RF", xgb.perf.rep[,4] ,  rf.perf.rep[,1], pow, bootstrap, digits)   
    if (sum(ensemble[c(2,6)]) > 0) {
      df = nested.compare0_0_6_2(df, group, "XGBoost lasso feature-RF with lasso feature", xgb.perf.rep[,5] ,  rf.perf.rep[,2], pow, bootstrap, digits)   
    } 
    if ( (sum(ensemble[c(3,4,7,8)]) > 0)  & (family == "gaussian") ) {
      df = nested.compare0_0_6_2(df, group, "XGBoost lasso offset-  RF with lasso offset", xgb.perf.rep[,6] ,  rf.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if ((doxgb == 1) & (doorf == 1)) {
    group = "13:XGB-ORF"
    df = nested.compare0_0_6_2(df, group, "XGBoost (tuned) - ORF", xgb.perf.rep[,4] ,  orf.perf.rep[,1], pow, bootstrap, digits)   
    if (sum(ensemble[c(2,6)]) > 0) {
      df = nested.compare0_0_6_2(df, group, "XGBoost lasso feature-ORF with lasso feature", xgb.perf.rep[,5] ,  orf.perf.rep[,2], pow, bootstrap, digits)   
    } 
    if ( (sum(ensemble[c(3,4,7,8)]) > 0)  & (family == "gaussian") ) {
      df = nested.compare0_0_6_2(df, group, "XGBoost lasso offset-  ORF with lasso offset", xgb.perf.rep[,6] ,  orf.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if ((doxgb == 1) & (doann == 1)) {
    group = "14:XGB-ANN"
    df = nested.compare0_0_6_2(df, group, "XGBoost (tuned) - ANN", xgb.perf.rep[,4] ,  ann.perf.rep[,1], pow, bootstrap, digits)   
    if (ensemble[6]) { 
      df = nested.compare0_0_6_2(df, group, "XGBoost lasso feature - ANN, l lasso feature", xgb.perf.rep[,5] ,  ann.perf.rep[,6], pow, bootstrap, digits)   
    } else if (ensemble[2]) { 
      df = nested.compare0_0_6_2(df, group, "XGBoost lasso feature - ANN lasso feature"   , xgb.perf.rep[,5] ,  ann.perf.rep[,2], pow, bootstrap, digits)   
    } 
    if (family == "gaussian") {
      if (ensemble[8]) { 
        df = nested.compare0_0_6_2(df, group, "XGBoost lasso offset-ANN l lasso offset(upated)", xgb.perf.rep[,6] ,  ann.perf.rep[,8], pow, bootstrap, digits)   
      } else if (ensemble[4]) { 
        df = nested.compare0_0_6_2(df, group, "XGBoost lasso offset - ANN, lasso offset"       , xgb.perf.rep[,6] ,  ann.perf.rep[,4], pow, bootstrap, digits)   
      } else if (ensemble[7]) { 
        df = nested.compare0_0_6_2(df, group, "XGBoost lasso offset - ANN l lasso offset"      , xgb.perf.rep[,6] ,  ann.perf.rep[,7], pow, bootstrap, digits)   
      } else if (ensemble[3]) { 
        df = nested.compare0_0_6_2(df, group, "XGBoost offset - ANN lasso offset"              , xgb.perf.rep[,6] ,  ann.perf.rep[,3], pow, bootstrap, digits)   
      }  
    }
  }
  
  if ((dorf == 1) & (doorf == 1)) {
    group = "15:RF-ORF"
    df = nested.compare0_0_6_2(df, group, "RF - ORF", rf.perf.rep[,1] ,  orf.perf.rep[,1], pow, bootstrap, digits) 
    if (sum(ensemble[c(2,6)]) > 0) {
      df = nested.compare0_0_6_2(df, group, "RF lasso feature - ORF l lasso feature", rf.perf.rep[,2] , orf.perf.rep[,2], pow, bootstrap, digits) 
    }  
    if (sum(ensemble[c(3,4,7,8)]) > 0) {
      df = nested.compare0_0_6_2(df, group, "RF lasso feature - ORF lasso feature"  , rf.perf.rep[,3] , orf.perf.rep[,3], pow, bootstrap, digits)  
    }
  }
  
  if ((dorf == 1) & (doann == 1)) {
    group = "16:RF-ANN"
    df = nested.compare0_0_6_2(df, group, "RF - ANN", rf.perf.rep[,1] ,  ann.perf.rep[,1], pow, bootstrap, digits) 
    if (ensemble[6]) {
      df = nested.compare0_0_6_2(df, group, "RF lasso feature - ANN  l lasso feature", rf.perf.rep[,2] ,  ann.perf.rep[,6], pow, bootstrap, digits) 
    } else if (ensemble[2]) {
      df = nested.compare0_0_6_2(df, group, "RF lasso feature - ANN lasso feature", rf.perf.rep[,2] ,  ann.perf.rep[,2], pow, bootstrap, digits)  
    }
    if (ensemble[8]) {
      df = nested.compare0_0_6_2(df, group, "RF lasso offset - ANN l lasso offset (upated)", rf.perf.rep[,3] ,  ann.perf.rep[,8], pow, bootstrap, digits)   
    } else if (ensemble[4]) {
      df = nested.compare0_0_6_2(df, group, "RF lasso offset - ANN lasso offset", rf.perf.rep[,3] ,  ann.perf.rep[,4], pow, bootstrap, digits)  
    } else if (ensemble[7]) {
      df = nested.compare0_0_6_2(df, group, "RF lasso offset - ANN, l lasso offset (upated)", rf.perf.rep[,3] ,  ann.perf.rep[,7], pow, bootstrap, digits) 
    } else if (ensemble[3]) {
      df = nested.compare0_0_6_2(df, group, "RF lasso offset - ANN, lasso offset", rf.perf.rep[,3] ,  ann.perf.rep[,3], pow, bootstrap, digits)  
    }
  }
  
  if ((doorf == 1) & (doann == 1)) {
    group = "17:ORF-ANN"
    df = nested.compare0_0_6_2(df, group, "ORF - ANN", orf.perf.rep[,1] ,  ann.perf.rep[,1], pow, bootstrap, digits) 
    if (ensemble[6]) {
      df = nested.compare0_0_6_2(df, group, "ORF lasso feature - ANN  l lasso feature", orf.perf.rep[,2] ,  ann.perf.rep[,6], pow, bootstrap, digits) 
    } else if (ensemble[2]) {
      df = nested.compare0_0_6_2(df, group, "ORF lasso feature - ANN lasso feature", orf.perf.rep[,2] ,  ann.perf.rep[,2], pow, bootstrap, digits)  
    }
    if (ensemble[8]) {
      df = nested.compare0_0_6_2(df, group, "ORF lasso offset - ANN l lasso offset (upated)", orf.perf.rep[,3] ,  ann.perf.rep[,8], pow, bootstrap, digits)   
    } else if (ensemble[4]) {
      df = nested.compare0_0_6_2(df, group, "ORF lasso offset - ANN lasso offset", orf.perf.rep[,3] ,  ann.perf.rep[,4], pow, bootstrap, digits)  
    } else if (ensemble[7]) {
      df = nested.compare0_0_6_2(df, group, "ORF lasso offset - ANN, l lasso offset (upated)", orf.perf.rep[,3] ,  ann.perf.rep[,7], pow, bootstrap, digits) 
    } else if (ensemble[3]) {
      df = nested.compare0_0_6_2(df, group, "ORF lasso offset - ANN, lasso offset", orf.perf.rep[,3] ,  ann.perf.rep[,3], pow, bootstrap, digits)  
    }
  }
  
#  print(df)
  
  if (table %in% c(1,2)) { 
    comma = ","
    df2 = cbind(df[,c(1:3)], " (", df[,4, drop=FALSE], comma, df[,5, drop=FALSE], ") ", df[,6, drop=FALSE])
    names(df2)[c(4,6,8)] = ""
    names(df2)[c(2,3,5,7,9)] = c("Comparison", type, "lower 95%", "upper 95%", "  p")
    
    ugroup = unique(df2$group)
    for (i_ in c(1:length(ugroup))) {
      df3 = df2[df2$group==ugroup[i_],-1]
      names(df3)[c(5,7)] = "" 
      print(df3, right=FALSE, row.names=FALSE)
      cat("\n")
    }
  } 
  if (table %in% c(0,2))  {
    row.names(df) = NULL 
    names(df) = c("Group", "Comparison", paste("mean",type), "lower 95%", "upper 95%", "  p")
    return(df)
  }
  
}

###############################################################################################################
###############################################################################################################
