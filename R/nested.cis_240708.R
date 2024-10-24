################################################################################
################################################################################

#' Calculate performance measure CI's and p's
#' 
#' @description 
#' Calculate performances as called from nested.cis().  
#'
#' @param a One term
#' @param mu null value
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4
#' @param txt 1 (default) to include inline text for estimated, 95 percent CI and p
#' @param pow Power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be 1 to 
#' show correlations.  Pow is ignored for the family of "cox" and "binomial".
#' 
#' @return An estimate, 95% CI and p
#' 
#' @importFrom stats t.test  
#' 
#' @noRd

nested.cis0 = function(a, mu, digits=4, txt=0, pow=1, alldevrat=NULL, bootstrap=0) { 
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
  
  if (txt==1) {
    cat ( paste0(  " estimate (95% CI): ", round(mean, digits=digits), " (", round(lo, digits=digits), ", ", 
                   round(up, digits=digits), ") , p=", round(p_, digits=digits) ) )
  } else {
    cat ( paste0( round(mean, digits=digits), " (", round(lo, digits=digits), ", ", 
                  round(up, digits=digits), ")   ", round(p_, digits=digits) ) )
  }
  #  if ( pow == 2) {cat("   --", corr1, " - ", corr2)}
}

################################################################################
################################################################################

#' Calculate performance measure CI's and p's
#'
#' @description 
#' Calculate overall estimates and confidence intervals for performance measures
#' based upon stored cross validation performance measures in a nested.glmnetr()
#' output object.
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

## object = xx ; digits=4 ; type="devrat" ; pow=1 ; returnd = 0 ; 
nested.cis = function(object, type="devrat", pow=1, digits=4, returnd=0) {
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
    } else if (type == "devian") {
      rpart.perf.rep  = object$rpart.devian.rep
    } else if (type == "lincal") {
      rpart.perf.rep  = object$rpart.lincal.rep
    } else if (type == "intcal") {
      rpart.perf.rep  = object$rpart.intcal.rep
    }
  }
  
  if ((dostep == 1) | (doaic==1)) {
    if (type == "agree") {
      step.perf.rep  = object$step.agree.rep
    } else if (type == "devrat") {
      devian.rep    = object$step.devian.rep 
      step.perf.rep  = IndDevRat(devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep)
      step.alldevrat = AllDevRat( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap )
      step.avedevrat = colMeans(step.perf.rep)
      step.sedevrat  = se_perf(step.perf.rep, bootstrap)
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
  
  if (sum(ensemble[c(2:8)]) > 0 ) {
    cat ("  Ensemble option used when fitting models : ")
    # ensemble\n" ) 
    cat(paste0("(", ensemble[1],",", ensemble[2],",", ensemble[3],",", ensemble[4],", ",
                 ensemble[5],",", ensemble[6],",", ensemble[7],",", ensemble[8],")\n\n")) 
  } 
  
  if ( ensemble[c(1)] == 0 ) {
    cat ("\n Simple models with informaiton from losso not run.  Output is abbreviated. \n\n" ) 
    doxgb = 0 ; dorf = 0 ; doann = 0 ; dorpart = 0 ; dostep = 0 ; doaic = 0 ; 
  }
  
  cat ("  Model performance evaluation in terms of", pm, "\n" )   
  cat ("  Null hypothesis value is", nullval, "\n") 

  if (dolasso == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = lasso.alldevrat }
    cat ("\n lasso.min                      ") ;  nested.cis0(lasso.perf.rep[,2], nullval, pow=pow, alldevrat=alldevrat[2], bootstrap=bootstrap) 
    cat ("\n lasso.minR                     ") ;  nested.cis0(lasso.perf.rep[,4], nullval, pow=pow, alldevrat=alldevrat[4], bootstrap=bootstrap)   
    cat ("\n lasso.minR0                    ") ;  nested.cis0(lasso.perf.rep[,6], nullval, pow=pow, alldevrat=alldevrat[6], bootstrap=bootstrap)   
    cat ("\n ridge                          ") ;  nested.cis0(lasso.perf.rep[,7], nullval, pow=pow, alldevrat=alldevrat[2], bootstrap=bootstrap)   
    cat("\n")
  }
  
  if (doxgb == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = xgb.alldevrat }
    cat ("\n XGBoost (simple)               ") ;  nested.cis0(xgb.perf.rep[,1] , nullval, pow=pow, alldevrat=alldevrat[1], bootstrap=bootstrap) ; 
    if (sum(ensemble[c(2,6)])    > 0) { cat ("\n XGBoost (simple) lasso feature ") ;  nested.cis0(xgb.perf.rep[,2] , nullval,pow=pow, alldevrat=alldevrat[2], bootstrap=bootstrap) ; }
    if (sum(ensemble[c(3,4,7,8)])> 0) { cat ("\n XGBoost (simple) lasso offset  ") ;  nested.cis0(xgb.perf.rep[,3] , nullval,pow=pow, alldevrat=alldevrat[3], bootstrap=bootstrap) ; }
    cat ("\n XGBoost (tuned)                ") ;  nested.cis0(xgb.perf.rep[,4] , nullval, pow=pow, alldevrat=alldevrat[4], bootstrap=bootstrap) ; 
    if (sum(ensemble[c(2,6)])    > 0) { cat ("\n XGBoost (tuned) lasso feature  ") ;  nested.cis0(xgb.perf.rep[,5] , nullval,pow=pow, alldevrat=alldevrat[5], bootstrap=bootstrap) ; }
    if (sum(ensemble[c(3,4,7,8)])> 0) { cat ("\n XGBoost (tuned) lasso offset   ") ;  nested.cis0(xgb.perf.rep[,6] , nullval,pow=pow, alldevrat=alldevrat[6], bootstrap=bootstrap) ; }
    cat("\n")
  }
  
  if (dorf == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = rf.alldevrat }
    cat ("\n RF                             ") ;  nested.cis0(rf.perf.rep[,1] , nullval, pow=pow, alldevrat=alldevrat[1], bootstrap=bootstrap) ;  
    if (sum(ensemble[c(2,6)] )> 0) { cat ("\n RF with lasso feature           ") ;  nested.cis0(rf.perf.rep[,2] , nullval, pow=pow, alldevrat=alldevrat[2], bootstrap=bootstrap) ;  }
    if ((sum(ensemble[c(3,4,7,8)] )> 0) & (family == "gaussian")) { cat ("\n RF with lasso offset            ") ;  nested.cis0(rf.perf.rep[,3] , nullval, pow=pow, alldevrat=alldevrat[3], bootstrap=bootstrap) ; }
    cat("\n")
  }

  if (doorf == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = orf.alldevrat }
    cat ("\n ORF                            ") ;  nested.cis0(orf.perf.rep[,1] , nullval, pow=pow, alldevrat=alldevrat[1], bootstrap=bootstrap) ;  
    if (sum(ensemble[c(2,6)] )> 0) { cat ("\n ORF with lasso feature          ") ;  nested.cis0(orf.perf.rep[,2] , nullval, pow=pow, alldevrat=alldevrat[2], bootstrap=bootstrap) ;  }
    if ((sum(ensemble[c(3,4,7,8)] )> 0) & (family == "gaussian")) { cat ("\n ORF with lasso offset           ") ;  nested.cis0(orf.perf.rep[,3] , nullval, pow=pow, alldevrat=alldevrat[3], bootstrap=bootstrap) ; }
    cat("\n")
  }
    
  if (doann == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = ann.alldevrat }
    cat ("\n ANN, no lasso info             ") ;  nested.cis0(ann.perf.rep[,1] , nullval, pow=pow, alldevrat=alldevrat[1], bootstrap=bootstrap)  
    if (sum(ensemble[2] )> 0) { cat ("\n ANN, with lasso feature          ") ;  nested.cis0(ann.perf.rep[,2] , nullval, pow=pow, alldevrat=alldevrat[2], bootstrap=bootstrap)  } 
    if (sum(ensemble[3] )> 0) { cat ("\n ANN, with lasso as offset        ") ;  nested.cis0(ann.perf.rep[,3] , nullval, pow=pow, alldevrat=alldevrat[3], bootstrap=bootstrap)  } 
    if (sum(ensemble[4] )> 0) { cat ("\n ANN, updated lasso as offset     ") ;  nested.cis0(ann.perf.rep[,4] , nullval, pow=pow, alldevrat=alldevrat[4], bootstrap=bootstrap)  } 
    if (sum(ensemble[5] )> 0) { cat ("\n ANN, lasso terms                 ") ;  nested.cis0(ann.perf.rep[,5] , nullval, pow=pow, alldevrat=alldevrat[5], bootstrap=bootstrap)  } 
    if (sum(ensemble[6] )> 0) { cat ("\n ANN, lasso terms, lasso feature  ") ;  nested.cis0(ann.perf.rep[,6] , nullval, pow=pow, alldevrat=alldevrat[6], bootstrap=bootstrap)  } 
    if (sum(ensemble[7] )> 0) { cat ("\n ANN, lasso terms, lasso as offset") ;  nested.cis0(ann.perf.rep[,7] , nullval, pow=pow, alldevrat=alldevrat[7], bootstrap=bootstrap)  } 
    if (sum(ensemble[8] )> 0) { cat ("\n ANN, lasso terms, updated lasso  ") ;  nested.cis0(ann.perf.rep[,8] , nullval, pow=pow, alldevrat=alldevrat[8], bootstrap=bootstrap)  } 
    cat("\n")  
  }
  
  if (dostep == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = step.alldevrat }
    cat ("\n step (df)                      ") ;  nested.cis0(step.perf.rep[,1] , nullval, pow=pow, alldevrat=alldevrat[1], bootstrap=bootstrap) ;
    cat ("\n step (p)                       ") ;  nested.cis0(step.perf.rep[,2] , nullval, pow=pow, alldevrat=alldevrat[2], bootstrap=bootstrap) ;
  }
  if (doaic == 1) {
    cat ("\n AIC                            ") ;  nested.cis0(step.perf.rep[,3] , nullval, pow=pow, alldevrat=alldevrat[3], bootstrap=bootstrap) ;
  }
  if ((dostep == 1) | (doaic == 1)) { cat("\n") }
  
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
  