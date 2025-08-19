################################################################################
################################################################################

#' Calculate performance measure CI's and p's
#' 
#' @description 
#' Calculate "nominal" performances as called from nested.cis(). In general 
#' the standard deviations for the performance measures evaluated 
#' on the leave-out samples may be biased. The confidence intervals and 
#' p values presented by this function should not be interpreted exactly 
#' but may still be of value. See the package vignettes for more discussion.  
#'
#' @param a One term
#' @param mu null value
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4
#' @param txt 1 (default) to include inline text for estimated, 95 percent CI and p
#' @param pow Power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be 1 to 
#' show correlations.  Pow is ignored for the family of "cox" and "binomial".
#' @param alldevrat 1 for optional derivation of deviance ratio. 
#' @param bootstrap 1 to indicate a bootstrap sample instead of cross validation.
#' 
#' @return An estimate, 95% CI and p
#' 
#' @importFrom stats t.test  
#' 
#' @noRd

nested.cis0_0_6_2 = function(df, group, desc, a, mu, digits=4, pow=1, alldevrat=NULL, bootstrap=0) { 
  if ( pow != 2) { pow = 1 } 

    azs = (a)[!is.na(a)]
    mean = mean(azs) 
    n_ = length(azs)
    if (bootstrap == 0) {
      se = sd(azs)/sqrt(n_) 
    } else if (bootstrap >= 1) {
      se = sd(azs) 
    }
    lo = mean - qt(0.975,df=(n_-1))*se
    up = mean + qt(0.975,df=(n_-1))*se
    t_ = (mean-mu)/se 
    p_ = 2*pt(-abs(t_), df=(n_-1)) 

  if ( pow == 2) {
    mean = sign(mean) * mean^2 
    lo = sign(lo) * lo^2 
    up = sign(up) * up^2 
  }
  
  if (!is.null(alldevrat)) {
    adj = alldevrat - mean                          ## bootstrap bias adjustment 
    mean = mean + adj
    lo   = lo + adj
    up   = up + adj
  }
  
#    cat ( paste0( round(mean, digits=digits), " (", round(lo, digits=digits), ", ", round(up, digits=digits), ")   ", round(p_, digits=digits) ) )

#    dfr = as.data.frame(list(group=group))  
    dfr = as.data.frame(list(group=group, desc=desc, mean=round(mean, digits=digits), lower95=round(lo, digits=digits), upper95=round(up, digits=digits), p=round(p_, digits=digits)))
    if (is.null(df)) { df = dfr 
    } else { df = rbind( df, dfr) }
    return( df )
}

################################################################################
################################################################################

#' Calculate performance measure "nominal" CI's and p's
#'
#' @description 
#' Calculate overall estimates and "nominal" confidence intervals for 
#' performance measures based upon stored cross validation performance 
#' measures in a nested.glmnetr() output object. The simple standard errors 
#' derived here from cross-validation are questionable and the confidence 
#' intervals coverage probabilities and p-values may be inaccruate.
#' See the Vignette references. 
#'
#' @param object A nested.glmnetr output object.
#' @param type determines what type of nested cross validation performance measures are 
#' compared.  Possible values are "devrat" to compare the deviance ratios, i.e. the 
#' fractional reduction in deviance relative to the null model deviance, 
#' "agree" to compare agreement, "lincal" to compare the linear calibration 
#' slope coefficients, "intcal" to compare the linear calibration intercept 
#' coefficients, from the nested cross validation. 
#' @param pow the power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  pow is ignored for the family of "cox" and "binomial".  When
#' pow = 2, calculations are made using correlations and the final estimates and 
#' confidence intervals are raised to the power of 2.  A negative sign before an 
#' R-square estimate or confidence limit indicates the estimate or confidence 
#' limit was negative before being raised to the power of 2.       
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4 
#' @param returnd 1 to return the deviance ratios in a list, 0 to not return.  The
#' deviances are stored in the nested.glmnetr() output object but not the deviance
#' ratios.  This function provides a simple mechanism to obtain the cross validated
#' deviance ratios. 
#' @param table 1 to print table to console, 0 to output the tabled information 
#' to a data frame
#'  
#' @return A printout to the R console 
#' 
#' @seealso
#'   \code{\link{nested.compare}} , \code{\link{summary.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
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
#' nested.cis(fit3)
#' }
#' 

## object = nested_glmnetr_fit_kp_ind3 ; digits=4 ; type="devrat" ; pow=1 ; returnd = 0 ; 
nested.cis = function(object, type="devrat", pow=1, digits=4, returnd=0, table=1) {
  
  cat(paste(" The standard errors (SEs) derived from cross-validation are questionable",
            "\n and the actual CI coverage probabilities and p-values may be inaccurate.", 
            "\n See references in the vignettes.\n\n"))
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
  if (is.na(doorf)) { doorf = 0 }
  ensemble = object$ensemble 
  
  null.m2LogLik.rep = object$null.m2LogLik.rep
  sat.m2LogLik.rep  = object$sat.m2LogLik.rep
  n.rep = object$n.rep
  
  IndDevRat = function(devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep ) {
    ( null.m2LogLik.rep - devian.rep ) / (null.m2LogLik.rep - sat.m2LogLik.rep)
  }
  
  AllDevRat = function( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap ) {
    if (bootstrap == 0) {
      AllDevRat = colSums( n.rep*( null.m2LogLik.rep - devian.rep ))/sum((null.m2LogLik.rep - sat.m2LogLik.rep)*n.rep )
    } else {
      AllDevRat = colMeans( ( null.m2LogLik.rep - devian.rep ) / (null.m2LogLik.rep - sat.m2LogLik.rep) ) 
    }
    return( AllDevRat )
  }
  
  se_perf = function(perf.rep, bootstrap) {
   if (bootstrap == 0) { 
    sqrt( apply(perf.rep, MARGIN =2, var)/length(perf.rep[,1]) ) 
   } else {
     sqrt( apply(perf.rep, MARGIN =2, var) )
   }
  }
  
  if (!(type %in% c("devrat", "devian", "agree", "intcal", "lincal"))) {
    cat( ' type must be one of c("devrat", "devian", "agree", "intcal", "lincal"). type is set to "devrat".\n\n')
    type = "devrat" 
  }

  if (type == "agree") { 
    if (family %in% c("cox", "binomial")) { nullval = 0.5 
    } else { nullval = 0 } 
  } else if (type == "devrat") {
    nullval = 0 
    devrat = list()
  } else if (type == "devian") {
    nullval = 0 
    devrat = list()
  } else if (type == "lincal") {
    nullval = 1 
  } else if (type == "intcal") {
    nullval = 0
  } else { 
    type == "devrat"
    nullval = 0 
    devrat = list()
  }
  
  if (dolasso == 1) {
    if (type == "agree") {
      lasso.perf.rep  = object$lasso.agree.rep
    } else if (type == "devrat") {
      devian.rep = object$lasso.devian.rep  
      lasso.perf.rep  = IndDevRat(devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep)
      lasso.alldevrat = AllDevRat( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap )
      lasso.avedevrat = colMeans(lasso.perf.rep)
      lasso.sedevrat  = se_perf(lasso.perf.rep, bootstrap)
      if (returnd == 1) { devrat$lasso = lasso.perf.rep ; devrat$lasso.all=lasso.alldevrat ; }
    } else if (type == "devian") {
      lasso.perf.rep  = object$lasso.devian.rep
    } else if (type == "lincal") {
      lasso.perf.rep  = object$lasso.lincal.rep
    } else if (type == "intcal") {
      lasso.perf.rep  = object$lasso.intcal.rep
    }
  }

  if (doxgb == 1) {
    if (type == "agree") {
      xgb.perf.rep  = object$xgb.agree.rep
    } else if (type == "devrat") {
      devian.rep    = object$xgb.devian.rep 
      xgb.perf.rep  = IndDevRat(devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep)
      xgb.alldevrat = AllDevRat( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap )
      xgb.avedevrat = colMeans(xgb.perf.rep)
      xgb.sedevrat  = se_perf(xgb.perf.rep, bootstrap)
      if (returnd == 1) { 
        devrat$xgb = xgb.perf.rep 
        devrat$xgb.all = xgb.alldevrat 
        if (sum(ensemble[c(2,6)]) == 0 ) { devrat$xgb[,c(2)] = 0 ; devrat$xgb.all[c(2)] = 0 ; }
        if (sum(ensemble[c(3,4,7,8)]) == 0 ) { devrat$xgb[,c(3)] = 0 ; devrat$xgb.all[c(3)] = 0 ; }
      }
    } else if (type == "devian") {
      xgb.perf.rep  = object$xgb.devian.rep
    } else if (type == "lincal") {
      xgb.perf.rep  = object$xgb.lincal.rep
    } else if (type == "intcal") {
      xgb.perf.rep  = object$xgb.intcal.rep
    }
  }
  
  if (dorf == 1) {
    if (type == "agree") {
      rf.perf.rep  = object$rf.agree.rep
    } else if (type == "devrat") {
      devian.rep    = object$rf.devian.rep 
      rf.perf.rep  = IndDevRat(devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep)
      rf.alldevrat = AllDevRat( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap )
      rf.avedevrat = colMeans(rf.perf.rep)
      rf.sedevrat  = se_perf(rf.perf.rep, bootstrap)
      if (returnd == 1) { 
        devrat$rf = rf.perf.rep 
        devrat$rf.all=rf.alldevrat 
        if (sum(ensemble[c(2,6)]) == 0 ) { devrat$rf[,c(2)] = 0 ; devrat$rf.all[c(2)] = 0 ; }
        if ( (sum(ensemble[c(3,4,7,8)]) == 0 ) | (family != "gaussian") ) { devrat$rf[,c(3)] = 0 ; devrat$rf.all[c(3)] = 0 ; }
      }
    } else if (type == "devian") {
      rf.perf.rep  = object$rf.devian.rep
    }  else if (type == "lincal") {
      rf.perf.rep  = object$rf.lincal.rep
    } else if (type == "intcal") {
      rf.perf.rep  = object$rf.intcal.rep
    }
  }
  
  if (doorf == 1) {
    if (type == "agree") {
      orf.perf.rep  = object$orf.agree.rep
    } else if (type == "devrat") {
      devian.rep    = object$orf.devian.rep 
      orf.perf.rep  = IndDevRat(devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep)
      orf.alldevrat = AllDevRat( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap )
      orf.avedevrat = colMeans(orf.perf.rep)
      orf.sedevrat  = se_perf(orf.perf.rep, bootstrap)
      if (returnd == 1) { 
        devrat$orf = orf.perf.rep 
        devrat$orf.all=orf.alldevrat 
        if (sum(ensemble[c(2,6)]) == 0 ) { devrat$orf[,c(2)] = 0 ; devrat$orf.all[c(2)] = 0 ; }
        if ( (sum(ensemble[c(3,4,7,8)]) == 0 ) | (family != "gaussian") ) { devrat$orf[,c(3)] = 0 ; devrat$orf.all[c(3)] = 0 ; }
      }
    } else if (type == "devian") {
      orf.perf.rep  = object$orf.devian.rep
    } else if (type == "lincal") {
      orf.perf.rep  = object$orf.lincal.rep
    } else if (type == "intcal") {
      orf.perf.rep  = object$orf.intcal.rep
    }
  }
  
  if (doann == 1) {
    if (type == "agree") {
      ann.perf.rep  = object$ann.agree.rep
    } else if (type == "devrat") {
      devian.rep    = object$ann.devian.rep 
      ann.perf.rep  = IndDevRat(devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep)
      ann.alldevrat = AllDevRat( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap )
      ann.avedevrat = colMeans(ann.perf.rep)
      ann.sedevrat  = se_perf(ann.perf.rep, bootstrap)
      if (returnd == 1) { 
        devrat$ann = ann.perf.rep 
        devrat$ann.all=ann.alldevrat 
        for ( i_ in c(2:8)) { if (sum(ensemble[c(i_)]) == 0 ) { devrat$ann[,c(i_)] = 0 ; devrat$ann.all[c(i_)] = 0 ; } }
      }
    } else if (type == "devian") {
      ann.perf.rep  = object$ann.devian.rep
    }  else if (type == "lincal") {
      ann.perf.rep  = object$ann.lincal.rep
    } else if (type == "intcal") {
      ann.perf.rep  = object$ann.intcal.rep
    }
  }
  
  if (dorpart == 1) {
    if (type == "agree") {
      rpart.perf.rep  = object$rpart.agree.rep
    } else if (type == "devrat") {
      devian.rep    = object$rpart.devian.rep 
      rpart.perf.rep  = IndDevRat(devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep)
      rpart.alldevrat = AllDevRat( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap )
      rpart.avedevrat = colMeans(rpart.perf.rep)
      rpart.sedevrat  = se_perf(rpart.perf.rep, bootstrap)
      if (returnd == 1) { 
        devrat$rpart = rpart.perf.rep ; 
        devrat$rpart.all=rpart.alldevrat
        if (sum(ensemble[c(2,6)])     == 0 ) { devrat$rpart[,c(2)] = 0 ; devrat$rpart.all[c(2)] = 0 ; }
        if (sum(ensemble[c(3,4,7,8)]) == 0 ) { devrat$rpart[,c(3)] = 0 ; devrat$rpart.all[c(3)] = 0 ; }
        }
    } else if (type == "devian") {
      rpart.perf.rep  = object$rpart.devian.rep
    } else if (type == "lincal") {
      rpart.perf.rep  = object$rpart.lincal.rep
    } else if (type == "intcal") {
      rpart.perf.rep  = object$rpart.intcal.rep
    }
  }
  
  if ((dostep == 1) | (doaic==1) | (dofull == 1)) {
    if (type == "agree") {
      step.perf.rep  = object$step.agree.rep
    } else if (type == "devrat") {
      devian.rep    = object$step.devian.rep 
      step.perf.rep  = IndDevRat(devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep)
      step.alldevrat = AllDevRat( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap )
      step.avedevrat = colMeans(step.perf.rep)
      step.sedevrat  = se_perf(step.perf.rep, bootstrap)
      if (returnd == 1) {
        devrat$step = step.perf.rep ; 
        devrat$step.all=step.alldevrat ; 
        if ((dostep != 1) ) { devrat$step[,c(1,2)] = 0 ; devrat$step.all[c(1,2)] = 0 ; }
        if ((doaic  != 1) ) { devrat$step[,c( 3 )] = 0 ; devrat$step.all[c( 3 )] = 0 ; }
        if ((dofull != 1) ) { devrat$step[,c( 4 )] = 0 ; devrat$step.all[c( 4 )] = 0 ; }
        if ((dostep != 1) ) { devrat$step[,c(1,2)] = 0 ; devrat$step.all[c(1,2)] = 0 ; }
        if ((doaic  != 1) ) { devrat$step[,c( 3 )] = 0 ; devrat$step.all[c( 3 )] = 0 ; }
        if ((dofull != 1) ) { devrat$step[,c( 4 )] = 0 ; devrat$step.all[c( 4 )] = 0 ; }
      }
    } else if (type == "devian") {
      step.perf.rep  = object$step.devian.rep
    }  else if (type == "lincal") {
      step.perf.rep  = object$step.lincal.rep
    } else if (type == "intcal") {
      step.perf.rep  = object$step.intcal.rep
    }
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

  if (table %in% c(1,2)) {  
    if (sum(ensemble[c(2:8)]) > 0 ) {
      cat ("  Ensemble option used when fitting models : ") 
      cat(paste0("(", ensemble[1],",", ensemble[2],",", ensemble[3],",", ensemble[4],", ",
                 ensemble[5],",", ensemble[6],",", ensemble[7],",", ensemble[8],")\n\n")) 
    }
    
    if ( ensemble[c(1)] == 0 ) {
      cat ("\n Simple models with informaiton from losso not run.  Output is abbreviated. \n\n" ) 
      doxgb = 0 ; dorf = 0 ; doann = 0 ; dorpart = 0 ; dostep = 0 ; doaic = 0 ; dofull = 0 ; 
    }
    
    cat ("  Model performance evaluation in terms of", pm, "\n\n" )   
    cat ("  Null hypothesis value is", nullval, "\n\n") 
  } 
  
  df = NULL 

  if (dolasso == 1) {
    group = "01:lasso"
    if (type == "devrat") { alldevrat = lasso.alldevrat 
    } else { alldevrat = NULL }
    if (substr(object$version[2],17, 21) %in% c("0.6-2")) {
#          nested.cis0_0_6_2(df, group, desc, a, mu, digits=4, pow=1, alldevrat=NULL, bootstrap=0) 
      df = nested.cis0_0_6_2(df, group, "lasso"           , lasso.perf.rep[,1], nullval, pow=pow, digits=digits, alldevrat=alldevrat[1], bootstrap=bootstrap) 
      df = nested.cis0_0_6_2(df, group, "lasso relaxed"   , lasso.perf.rep[,2], nullval, pow=pow, digits=digits, alldevrat=alldevrat[2], bootstrap=bootstrap)   
      df = nested.cis0_0_6_2(df, group, "lasso relaxed G0", lasso.perf.rep[,3], nullval, pow=pow, digits=digits, alldevrat=alldevrat[3], bootstrap=bootstrap)
      if (!is.null(object$cv_elastic_fit)) {
        df = nested.cis0_0_6_2(df, group, "elastic net", lasso.perf.rep[,4], nullval, pow=pow, digits=digits, alldevrat=alldevrat[4], bootstrap=bootstrap) }
      df = nested.cis0_0_6_2(df, group, "ridge"        , lasso.perf.rep[,5], nullval, pow=pow, digits=digits, alldevrat=alldevrat[5], bootstrap=bootstrap)
      if (!is.null(object$cv_elastic_fit)) {
        df = nested.cis0_0_6_2(df, group, "elastic net G0", lasso.perf.rep[,6], nullval, pow=pow, digits=digits, alldevrat=alldevrat[6], bootstrap=bootstrap)
        df = nested.cis0_0_6_2(df, group, "elastic net G1", lasso.perf.rep[,7], nullval, pow=pow, digits=digits, alldevrat=alldevrat[7], bootstrap=bootstrap)
      }
    } else if (substr(object$version[2],17, 21) %in% c("0.6-1")) {
      df = nested.cis0_0_6_2(df, group, "lasso.min"       , lasso.perf.rep[,1], nullval, pow=pow, digits=digits, alldevrat=alldevrat[1], bootstrap=bootstrap) 
      df = nested.cis0_0_6_2(df, group, "lasso relaxed"   , lasso.perf.rep[,2], nullval, pow=pow, digits=digits, alldevrat=alldevrat[2], bootstrap=bootstrap)   
      df = nested.cis0_0_6_2(df, group, "lasso relaxed G0", lasso.perf.rep[,3], nullval, pow=pow, digits=digits, alldevrat=alldevrat[3], bootstrap=bootstrap)
      if (!is.null(object$cv_elastic_fit)) {
        df =nested.cis0_0_6_2(df, group, "elastic net", lasso.perf.rep[,4], nullval, pow=pow, digits=digits, alldevrat=alldevrat[4], bootstrap=bootstrap)
      }
      df = nested.cis0_0_6_2(df, group, "ridge"      , lasso.perf.rep[,5], nullval, pow=pow, digits=digits, alldevrat=alldevrat[5], bootstrap=bootstrap)
    } else {
      df = nested.cis0_0_6_2(df, group, "lasso"           , lasso.perf.rep[,2], nullval, pow=pow, digits=digits, alldevrat=alldevrat[2], bootstrap=bootstrap) 
      df = nested.cis0_0_6_2(df, group, "lasso relaxed"   , lasso.perf.rep[,4], nullval, pow=pow, digits=digits, alldevrat=alldevrat[4], bootstrap=bootstrap)   
      df = nested.cis0_0_6_2(df, group, "lasso relaxed G0", lasso.perf.rep[,6], nullval, pow=pow, digits=digits, alldevrat=alldevrat[6], bootstrap=bootstrap)   
      df = nested.cis0_0_6_2(df, group, "ridge"           , lasso.perf.rep[,7], nullval, pow=pow, digits=digits, alldevrat=alldevrat[7], bootstrap=bootstrap)   ## alldevrat[2] or alldevrat[7] ??
    }
  }
  
  if (doxgb == 1) {
    group = "02:xgb"
    alldevrat = NULL
    if (type == "devrat") { alldevrat = xgb.alldevrat }
    df =                                     nested.cis0_0_6_2(df, group, "XGBoost (simple)"              , xgb.perf.rep[,1], nullval, pow=pow, digits=digits, alldevrat=alldevrat[1], bootstrap=bootstrap) 
    if (sum(ensemble[c(2,6)])    > 0) { df = nested.cis0_0_6_2(df, group, "XGBoost (simple) lasso feature", xgb.perf.rep[,2], nullval, pow=pow, digits=digits, alldevrat=alldevrat[2], bootstrap=bootstrap) }
    if (sum(ensemble[c(3,4,7,8)])> 0) { df = nested.cis0_0_6_2(df, group, "XGBoost (simple) lasso offset" , xgb.perf.rep[,3], nullval, pow=pow, digits=digits, alldevrat=alldevrat[3], bootstrap=bootstrap) }
    df =                                     nested.cis0_0_6_2(df, group, "XGBoost (tuned)"               , xgb.perf.rep[,4], nullval, pow=pow, digits=digits, alldevrat=alldevrat[4], bootstrap=bootstrap) 
    if (sum(ensemble[c(2,6)])    > 0) { df = nested.cis0_0_6_2(df, group, "XGBoost (tuned) lasso feature" , xgb.perf.rep[,5], nullval, pow=pow, digits=digits, alldevrat=alldevrat[5], bootstrap=bootstrap) }
    if (sum(ensemble[c(3,4,7,8)])> 0) { df = nested.cis0_0_6_2(df, group, "XGBoost (tuned) lasso offset"  , xgb.perf.rep[,6], nullval, pow=pow, digits=digits, alldevrat=alldevrat[6], bootstrap=bootstrap) }
  }
  
  if (dorf == 1) {
    group = "03:RF"
    alldevrat = NULL
    if (type == "devrat") { alldevrat = rf.alldevrat }
    df =                                                                 nested.cis0_0_6_2(df, group, "RF"                   , rf.perf.rep[,1] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[1], bootstrap=bootstrap) ;  
    if (sum(ensemble[c(2,6)] )> 0) {                                df = nested.cis0_0_6_2(df, group, "RF with lasso feature", rf.perf.rep[,2] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[2], bootstrap=bootstrap) }
    if ((sum(ensemble[c(3,4,7,8)] )> 0) & (family == "gaussian")) { df = nested.cis0_0_6_2(df, group, "RF with lasso offset" , rf.perf.rep[,3] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[3], bootstrap=bootstrap) }
  }

  if (doorf == 1) {
    group = "04:ORF"
    alldevrat = NULL
    if (type == "devrat") { alldevrat = orf.alldevrat }
    df =                                                                 nested.cis0_0_6_2(df, group, "ORF"                   , orf.perf.rep[,1] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[1], bootstrap=bootstrap) ;  
    if (sum(ensemble[c(2,6)] )> 0) {                                df = nested.cis0_0_6_2(df, group, "ORF with lasso feature", orf.perf.rep[,2] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[2], bootstrap=bootstrap) ;  }
    if ((sum(ensemble[c(3,4,7,8)] )> 0) & (family == "gaussian")) { df = nested.cis0_0_6_2(df, group, "ORF with lasso offset" , orf.perf.rep[,3] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[3], bootstrap=bootstrap) ; }
  }
    
  if (doann == 1) {
    group = "05:ANN"
    alldevrat = NULL
    if (type == "devrat") { alldevrat = ann.alldevrat }
    df =                             nested.cis0_0_6_2(df, group, "ANN, no lasso info"               , ann.perf.rep[,1] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[1], bootstrap=bootstrap)  
    if (sum(ensemble[2] )> 0) { df = nested.cis0_0_6_2(df, group, "ANN, with lasso feature"          , ann.perf.rep[,2] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[2], bootstrap=bootstrap)  } 
    if (sum(ensemble[3] )> 0) { df = nested.cis0_0_6_2(df, group, "ANN, with lasso as offset"        , ann.perf.rep[,3] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[3], bootstrap=bootstrap)  } 
    if (sum(ensemble[4] )> 0) { df = nested.cis0_0_6_2(df, group, "ANN, updated lasso as offset"     , ann.perf.rep[,4] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[4], bootstrap=bootstrap)  } 
    if (sum(ensemble[5] )> 0) { df = nested.cis0_0_6_2(df, group, "ANN, lasso terms"                 , ann.perf.rep[,5] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[5], bootstrap=bootstrap)  } 
    if (sum(ensemble[6] )> 0) { df = nested.cis0_0_6_2(df, group, "ANN, lasso terms, lasso feature"  , ann.perf.rep[,6] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[6], bootstrap=bootstrap)  } 
    if (sum(ensemble[7] )> 0) { df = nested.cis0_0_6_2(df, group, "ANN, lasso terms, lasso as offset", ann.perf.rep[,7] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[7], bootstrap=bootstrap)  } 
    if (sum(ensemble[8] )> 0) { df = nested.cis0_0_6_2(df, group, "ANN, lasso terms, updated lasso"  , ann.perf.rep[,8] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[8], bootstrap=bootstrap)  } 
  }
  
  if (dostep == 1) {
    group = "05:step"
    alldevrat = NULL
    if (type == "devrat") { alldevrat = step.alldevrat }
    df = nested.cis0_0_6_2(df, group, "step (df)", step.perf.rep[,1] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[1], bootstrap=bootstrap) ;
    df = nested.cis0_0_6_2(df, group, "step (p)" , step.perf.rep[,2] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[2], bootstrap=bootstrap) ;
  }
  
  if (doaic == 1) {
    group = "05:step"
    df = nested.cis0_0_6_2(df, group, "AIC", step.perf.rep[,3] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[3], bootstrap=bootstrap) ;
  }
  
  if (dofull == 1) {
    group = "06:full"
    alldevrat = NULL
    if (type == "devrat") { alldevrat = step.alldevrat }
    df = nested.cis0_0_6_2(df, group, "full model, no penalty", step.perf.rep[,4] , nullval, pow=pow, digits=digits, alldevrat=alldevrat[4], bootstrap=bootstrap) ;
  }
#  if ((dostep == 1) | (doaic == 1)) { cat("\n") }
  
  if ( (returnd == 1) & (table %in% c(0,2)) ) {
    cat(" Only one of returnd and table can indicae return of an output object\n",
        " returnd is set to 0\n\n")
    returnd = 0 
  }
  
  if (table %in% c(1,2)) { 
    comma = ","
    df2 = cbind(df[,c(1:3)], " (", df[,4, drop=FALSE], comma, df[,5, drop=FALSE], ") ", df[,6, drop=FALSE])
    names(df2)[c(4,6,8)] = ""
    names(df2)[c(2,3,5,7,9)] = c("Model", type, "lower 95%", "upper 95%", "  p")
    
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
    names(df) = c("Model", "desc", paste("mean",type), "lower 95%", "upper 95%", "  p")
    return(df)
  }
  
  if ((type == "devrat") & (returnd == 1)) { return(devrat) }
}

####################################################################################################################################
####################################################################################################################################

################################################################################
################################################################################

#' A redirect to nested.cis() 
#'
#' @description 
#' See nested.cis(), glmnetr.cis() is depricated
#'
#' @param object A nested.glmnetr output object.
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4 
#' @param pow the power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  pow is ignored for the family of "cox" and "binomial".  When
#' pow = 2, calculations are made using correlations and the final estimates and 
#' confidence intervals are raised to the power of 2.  A negative sign before an 
#' R-square estimate or confidence limit indicates the estimate or confidence 
#' limit was negative before being raised to the power of 2.       
#' @param type determines what type of nested cross validation performance measures are 
#' compared.  Possible values are "devrat" to compare the deviance ratios, i.e. the 
#' fractional reduction in deviance relative to the null model deviance, 
#' "agree" to compare agreement, "lincal" to compare the linear calibration 
#' slope coefficients, "intcal" to compare the linear calibration intercept 
#' coefficients, from the nested cross validation. 
#' @param returnd 1 to return the deviance ratios in a list, 0 to not return.  The
#' deviances are stored in the nested.glmnetr() output object but not the deviance
#' ratios.  This function provides a simple mechanism to obtain the cross validated
#' deviance ratios. 
#'  
#' @return A printout to the R console 
#' 
#' @export
#'

glmnetr.cis = function(object, type="devrat", pow=1, digits=4, returnd=0) {
  nested.cis(object, type, pow, digits, returnd)  
} 
  