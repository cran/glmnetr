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

nested.cis0 = function(a, mu, digits=4, txt=0, pow=1, alldevrat=NULL) { 
  if ( pow != 2) { pow = 1 } 
  tdiff = t.test(a, mu=mu)
  mean = tdiff$estimate
  lo = tdiff$conf.int[1]
  up = tdiff$conf.int[2]
  p_ = tdiff$p.value
  if ( pow == 2) {
    mean = sign(mean) * mean^2 
    lo = sign(lo) * lo^2 
    up = sign(up) * up^2 
  }
  
  if (!is.null(alldevrat)) {
    dif = alldevrat - mean
    mean = mean + dif
    lo   = lo + dif
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
  
  m2.ll.null = object$null.m2LogLik.cv
  m2.ll.sat  = object$sat.m2LogLik.cv
  n__ = object$n.cv
  
  if (type == "agree") { 
    if (family %in% c("cox", "binomial")) { nullval = 0.5 
    } else { nullval = 0 } 
  } else if (type == "devrat") {
    nullval = 0 
    devrat = list()
  } else if (type == "lincal") {
    nullval = 1 
  } else if (type == "intcal") {
    nullval = 0
  }
  
  if (dolasso == 1) {
    if (type == "agree") {
      lasso.perf.cv  = object$lasso.agree.cv
    } else if (type == "devrat") {
      object1  = object$lasso.devian.cv ; object2 = 0 ; object3 = 0 ; object4 = 1 ;
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) 
        object1[,j_] = devratl[[j_]][[1]]
        object2[j_] = devratl[[j_]][[2]]
        object3[j_] = devratl[[j_]][[3]]
        object4[j_] = devratl[[j_]][[4]]
      }
      lasso.perf.cv   = object1
      lasso.alldevrat = object2
      lasso.avedevrat = object3
      lasso.sedevrat  = object4
    } else if (type == "lincal") {
      lasso.perf.cv  = object$lasso.lincal.cv
    } else if (type == "intcal") {
      lasso.perf.cv  = object$lasso.intcal.cv
    }
  }

  if (doxgb == 1) {
    if (type == "agree") {
      xgb.perf.cv  = object$xgb.agree.cv
    } else if (type == "devrat") {
      object1  = object$xgb.devian.cv ; object2 = 0 ; object3 = 0 ; object4 = 1 ;
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) 
        object1[,j_] = devratl[[j_]][[1]]
        object2[j_]  = devratl[[j_]][[2]]
      }
      xgb.perf.cv = object1
      xgb.alldevrat = object2
    } else if (type == "lincal") {
      xgb.perf.cv  = object$xgb.lincal.cv
    } else if (type == "intcal") {
      xgb.perf.cv  = object$xgb.intcal.cv
    }
  }
  
  if (dorf == 1) {
    if (type == "agree") {
      rf.perf.cv  = object$rf.agree.cv
    } else if (type == "devrat") {
      object1  = object$rf.devian.cv ;  object2 = 0 ; object3 = 0 ; object4 = 1 ;
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) 
        object1[,j_] = devratl[[j_]][[1]]
        object2[j_]  = devratl[[j_]][[2]]
      }
      rf.perf.cv = object1
      rf.alldevrat = object2
    }  else if (type == "lincal") {
      rf.perf.cv  = object$rf.lincal.cv
    } else if (type == "intcal") {
      rf.perf.cv  = object$rf.intcal.cv
    }
  }
  
  if (doorf == 1) {
    if (type == "agree") {
      orf.perf.cv  = object$orf.agree.cv
    } else if (type == "devrat") {
      object1  = object$orf.devian.cv ;  object2 = 0 ; object3 = 0 ; object4 = 1 ;
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) 
        object1[,j_] = devratl[[j_]][[1]]
        object2[j_]  = devratl[[j_]][[2]]
      }
      orf.perf.cv = object1
      orf.alldevrat = object2
    }  else if (type == "lincal") {
      orf.perf.cv  = object$orf.lincal.cv
    } else if (type == "intcal") {
      orf.perf.cv  = object$orf.intcal.cv
    }
  }
  
  if (doann == 1) {
    if (type == "agree") {
      ann.perf.cv  = object$ann.agree.cv
    } else if (type == "devrat") {
      object1  = object$ann.devian.cv ;  object2 = 0 ; object3 = 0 ; object4 = 1 ;
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) 
        object1[,j_] = devratl[[j_]][[1]]
        object2[j_]  = devratl[[j_]][[2]]
      }
      ann.perf.cv = object1
      ann.alldevrat = object2
    }  else if (type == "lincal") {
      ann.perf.cv  = object$ann.lincal.cv
    } else if (type == "intcal") {
      ann.perf.cv  = object$ann.intcal.cv
    }
  }
  
  if (dorpart == 1) {
    if (type == "agree") {
      rpart.perf.cv  = object$rpart.agree.cv
    } else if (type == "devrat") {
      object1  = object$rpart.devian.cv ;  object2 = 0 ; object3 = 0 ; object4 = 1 ;
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) 
        object1[,j_] = devratl[[j_]][[1]]
        object2[j_]  = devratl[[j_]][[2]]
      }
      rpart.perf.cv = object1
      rpart.alldevrat = object2
    } else if (type == "lincal") {
      rpart.perf.cv  = object$rpart.lincal.cv
    } else if (type == "intcal") {
      rpart.perf.cv  = object$rpart.intcal.cv
    }
  }
  
  if ((dostep == 1) | (doaic==1)) {
    if (type == "agree") {
      step.perf.cv  = object$step.agree.cv
    } else if (type == "devrat") {
      object1  = object$step.devian.cv ;  object2 = 0 ; object3 = 0 ; object4 = 1 ;
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) 
        object1[,j_] = devratl[[j_]][[1]]
        object2[j_]  = devratl[[j_]][[2]]
      }
      step.perf.cv = object1
      step.alldevrat = object2
    }  else if (type == "lincal") {
      step.perf.cv  = object$step.lincal.cv
    } else if (type == "intcal") {
      step.perf.cv  = object$step.intcal.cv
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
    cat ("\n Simple models with informaiton from lossa not run.  Output is abbreviated. \n\n" ) 
    doxgb = 0 ; dorf = 0 ; doann = 0 ; dorpart = 0 ; dostep = 0 ; doaic = 0 ; 
  }
  
  cat ("  Model performance evaluation in terms of", pm, "\n" )   
  cat ("  Null hypothesis value is", nullval, "\n") 

  if (dolasso == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = lasso.alldevrat }
    cat ("\n lasso.min                      ") ;  nested.cis0(lasso.perf.cv[,2], nullval, pow=pow, alldevrat=alldevrat[2]) 
    cat ("\n lasso.minR                     ") ;  nested.cis0(lasso.perf.cv[,4], nullval, pow=pow, alldevrat=alldevrat[4])   
    cat ("\n lasso.minR0                    ") ;  nested.cis0(lasso.perf.cv[,6], nullval, pow=pow, alldevrat=alldevrat[6])   
    cat ("\n ridge                          ") ;  nested.cis0(lasso.perf.cv[,7], nullval, pow=pow, alldevrat=alldevrat[2])   
    cat("\n")
  }
  
  if (doxgb == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = xgb.alldevrat }
    cat ("\n XGBoost (simple)               ") ;  nested.cis0(xgb.perf.cv[,1] , nullval, pow=pow, alldevrat=alldevrat[1]) ; 
    if (sum(ensemble[c(2,6)])    > 0) { cat ("\n XGBoost (simple) lasso feature ") ;  nested.cis0(xgb.perf.cv[,2] , nullval,pow=pow, alldevrat=alldevrat[2]) ; }
    if (sum(ensemble[c(3,4,7,8)])> 0) { cat ("\n XGBoost (simple) lasso offset  ") ;  nested.cis0(xgb.perf.cv[,3] , nullval,pow=pow, alldevrat=alldevrat[3]) ; }
    cat ("\n XGBoost (tuned)                ") ;  nested.cis0(xgb.perf.cv[,4] , nullval, pow=pow, alldevrat=alldevrat[4]) ; 
    if (sum(ensemble[c(2,6)])    > 0) { cat ("\n XGBoost (tuned) lasso feature  ") ;  nested.cis0(xgb.perf.cv[,5] , nullval,pow=pow, alldevrat=alldevrat[5]) ; }
    if (sum(ensemble[c(3,4,7,8)])> 0) { cat ("\n XGBoost (tuned) lasso offset   ") ;  nested.cis0(xgb.perf.cv[,6] , nullval,pow=pow, alldevrat=alldevrat[6]) ; }
    cat("\n")
  }
  
  if (dorf == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = rf.alldevrat }
    cat ("\n RF                             ") ;  nested.cis0(rf.perf.cv[,1] , nullval, pow=pow, alldevrat=alldevrat[1]) ;  
    if (sum(ensemble[c(2,6)])> 0) { cat ("\n RF with lasso feature           ") ;  nested.cis0(rf.perf.cv[,2] , nullval, pow=pow, alldevrat=alldevrat[2]) ;  }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) { cat ("\n RF with lasso offset            ") ;  nested.cis0(rf.perf.cv[,3] , nullval, pow=pow, alldevrat=alldevrat[3]) ; }
    cat("\n")
  }

  if (doorf == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = orf.alldevrat }
    cat ("\n ORF                            ") ;  nested.cis0(orf.perf.cv[,1] , nullval, pow=pow, alldevrat=alldevrat[1]) ;  
    if (sum(ensemble[c(2,6)])> 0) { cat ("\n ORF with lasso feature          ") ;  nested.cis0(orf.perf.cv[,2] , nullval, pow=pow, alldevrat=alldevrat[2]) ;  }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) { cat ("\n ORF with lasso offset           ") ;  nested.cis0(orf.perf.cv[,3] , nullval, pow=pow, alldevrat=alldevrat[3]) ; }
    cat("\n")
  }
    
  if (doann == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = ann.alldevrat }
    cat ("\n ANN, no lasso info             ") ;  nested.cis0(ann.perf.cv[,1] , nullval, pow=pow, alldevrat=alldevrat[1])  
    if (sum(ensemble[2])> 0) { cat ("\n ANN, with lasso feature          ") ;  nested.cis0(ann.perf.cv[,2] , nullval, pow=pow, alldevrat=alldevrat[2])  } 
    if (sum(ensemble[3])> 0) { cat ("\n ANN, with lasso as offset        ") ;  nested.cis0(ann.perf.cv[,3] , nullval, pow=pow, alldevrat=alldevrat[3])  } 
    if (sum(ensemble[4])> 0) { cat ("\n ANN, updated lasso as offset     ") ;  nested.cis0(ann.perf.cv[,4] , nullval, pow=pow, alldevrat=alldevrat[4])  } 
    if (sum(ensemble[5])> 0) { cat ("\n ANN, lasso terms                 ") ;  nested.cis0(ann.perf.cv[,5] , nullval, pow=pow, alldevrat=alldevrat[5])  } 
    if (sum(ensemble[6])> 0) { cat ("\n ANN, lasso terms, lasso feature  ") ;  nested.cis0(ann.perf.cv[,6] , nullval, pow=pow, alldevrat=alldevrat[6])  } 
    if (sum(ensemble[7])> 0) { cat ("\n ANN, lasso terms, lasso as offset") ;  nested.cis0(ann.perf.cv[,7] , nullval, pow=pow, alldevrat=alldevrat[7])  } 
    if (sum(ensemble[8])> 0) { cat ("\n ANN, lasso terms, updated lasso  ") ;  nested.cis0(ann.perf.cv[,8] , nullval, pow=pow, alldevrat=alldevrat[8])  } 
    cat("\n")  
  }
  
  if (dostep == 1) {
    alldevrat = NULL
    if (type == "devrat") { alldevrat = step.alldevrat }
    cat ("\n step (df)                      ") ;  nested.cis0(step.perf.cv[,1] , nullval, pow=pow, alldevrat=alldevrat[1]) ;
    cat ("\n step (p)                       ") ;  nested.cis0(step.perf.cv[,2] , nullval, pow=pow, alldevrat=alldevrat[2]) ;
  }
  if (doaic == 1) {
    cat ("\n AIC                            ") ;  nested.cis0(step.perf.cv[,3] , nullval, pow=pow, alldevrat=alldevrat[3]) ;
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
  