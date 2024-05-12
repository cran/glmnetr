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
#' 
#' @return An estimate, 95% CI and p for agreement comparison 
#' 
#' @importFrom stats t.test qt pt var 
#'
#' @noRd

nested.compare0 = function(a, b, digits=4, txt=0, pow=1) { 
  if ( pow != 2) { pow = 1 } 
  if (pow == 1) {
    tdiff = t.test(a-b)
    mean = tdiff$estimate
    lo = tdiff$conf.int[1]
    up = tdiff$conf.int[2]
    p_ = tdiff$p.value
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
    p_ = 2*min( pt(t_,n_-1), pt(-t_,n_-1) )
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

###############################################################################################################
###############################################################################################################

#' Compare cross validation fit performances from a nested.glmnetr output.
#'
#' @description 
#' Compare cross-validation model fits in terms of average performances from the 
#' nested cross validation fits.  
#'
#' @param object A nested.glmnetr output object.
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4 
#' @param pow the power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  pow is ignored for the family of "cox" and "binomial".  
#' @param type determines what type of nested cross validation performance measures are 
#' compared.  Possible values are "devrat" to compare the deviance ratios, i.e. the 
#' fractional reduction in deviance relative to the null model deviance, 
#' "agree" to compare agreement, "lincal" to compare the linear calibration 
#' slope coefficients, "intcal" to compare the linear calibration intercept 
#' coefficients, from the nested cross validation. 
#'  
#' @return A printout to the R console. 
#' 
#' @seealso
#'   \code{\link{nested.cis}} , \code{\link{summary.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
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
#' nested.compare(fit3)
#' }
#' 
nested.compare = function(object, digits=4, type="devrat", pow=1) {
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
  
  if (dolasso == 1) {
    if (type == "agree") {
      lasso.perf.cv  = object$lasso.agree.cv
    } else if (type == "devrat") {
      object1  = object$lasso.devian.cv
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      lasso.perf.cv = object1
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
      object1  = object$xgb.devian.cv
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      xgb.perf.cv = object1
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
      object1  = object$rf.devian.cv
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      rf.perf.cv = object1
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
      object1  = object$orf.devian.cv
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      orf.perf.cv = object1
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
      object1  = object$ann.devian.cv
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      ann.perf.cv = object1
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
      object1  = object$rpart.devian.cv
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      rpart.perf.cv = object1
    } else if (type == "lincal") {
      rpart.perf.cv  = object$rpart.lincal.cv
    } else if (type == "intcal") {
      rpart.perf.cv  = object$rpart.intcal.cv
    }
  }
  
  if (dostep == 1) {
    if (type == "agree") {
      step.perf.cv  = object$step.agree.cv
    } else if (type == "devrat") {
      object1  = object$step.devian.cv
      devratl = list() 
      for (j_ in c(1:dim(object1)[2])) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
      step.perf.cv = object1
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
    cat ("\n Simple models with information from losso not run.  Output is abbreviated. \n\n" ) 
    doxgb = 0 ; dorf = 0 ; doorf = 0 ; doann = 0 ; dorpart = 0 ; dostep = 0 ; doaic = 0 ; 
  }
  
  cat ("  Model performance comparison in terms of", pm, "\n\n" )   
  cat ("  Comparison                                estimate   (95% CI)         p\n") 
  
  if (dolasso == 1) {
    cat ("\n lasso.minR  - lasso.min                     ") ;  nested.compare0(lasso.perf.cv[,4] , lasso.perf.cv[,2],pow=pow) 
    cat ("\n lasso.minR  - lasso.minR0                   ") ;  nested.compare0(lasso.perf.cv[,4] , lasso.perf.cv[,6],pow=pow)   
    cat ("\n lasso.min   - lasso.minR0                   ") ;  nested.compare0(lasso.perf.cv[,2] , lasso.perf.cv[,6],pow=pow)   
    cat("\n")
  }
  
  #  print(xgb.perf.cv)
  
  if (doxgb == 1) {
    cat ("\n XGBoost (tuned) - XGBoost (simple)          ") ;  nested.compare0(xgb.perf.cv[,4] , xgb.perf.cv[,1],pow=pow) ; 
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n XGBoost (tuned) lasso feature - no feature  ") ;  nested.compare0(xgb.perf.cv[,5] , xgb.perf.cv[,4],pow=pow) ;  
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      cat ("\n XGBoost (tuned) lasso offset - no offset    ") ;  nested.compare0(xgb.perf.cv[,6] , xgb.perf.cv[,4],pow=pow) ;  
    }
    cat("\n")
  }
  
  if (dorf == 1) {
    lr = 0 
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n RF with lasso feature - no feature          ") ;  nested.compare0(rf.perf.cv[,2] , rf.perf.cv[,1],pow=pow) ;  
      lr = 1 
    }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) {
      cat ("\n RF with lasso offset - no offset            ") ;  nested.compare0(rf.perf.cv[,3] , rf.perf.cv[,1],pow=pow) ;  
      lr = 1 
    }
    if (lr == 1) { cat("\n") } 
  }
  
  if (doorf == 1) {
    lr = 0 
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n ORF with lasso feature - no feature          ") ;  nested.compare0(orf.perf.cv[,2] , orf.perf.cv[,1],pow=pow) ;  
      lr = 1 
    }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) {
      cat ("\n ORF with lasso offset - no offset            ") ;  nested.compare0(orf.perf.cv[,3] , orf.perf.cv[,1],pow=pow) ;  
      lr = 1 
    }
    if (lr == 1) { cat("\n") } 
  }
  
  if (doann == 1) {
    lr = 0 
    if (sum(ensemble[6])> 0) {
      cat ("\n ANN with with lasso feature - no feature    ") ;  nested.compare0(ann.perf.cv[,6] , ann.perf.cv[,1],pow=pow) ; lr = 1 
    } else if (sum(ensemble[2])> 0) {
      cat ("\n ANN with with lasso feature - no feature    ") ;  nested.compare0(ann.perf.cv[,2] , ann.perf.cv[,1],pow=pow) ; lr = 1 
      
    } 
    if (sum(ensemble[8])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  nested.compare0(ann.perf.cv[,8] , ann.perf.cv[,1],pow=pow) ; lr = 1 
    } else if (sum(ensemble[7])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  nested.compare0(ann.perf.cv[,7] , ann.perf.cv[,1],pow=pow) ; lr = 1  
    } else     if (sum(ensemble[4])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  nested.compare0(ann.perf.cv[,4] , ann.perf.cv[,1],pow=pow) ; lr = 1 
    } else     if (sum(ensemble[3])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  nested.compare0(ann.perf.cv[,3] , ann.perf.cv[,1],pow=pow) ; lr = 1 
    } 
    if (lr == 1) { cat("\n") } 
  }
  
  if (dostep == 1) {
    cat ("\n step (df) - step (p)                        ") ;  nested.compare0(step.perf.cv[,1]      , step.perf.cv[,2],pow=pow)    ;  cat("\n")
  }
  
  #  cat("\n")
  
  if ((dolasso == 1) & (doxgb == 1)) {
    cat ("\n lasso.minR - XGB (tuned)                    ") ;  nested.compare0(lasso.perf.cv[,4] , xgb.perf.cv[,4],pow=pow) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - XGB with lasso feature         ") ;  nested.compare0(lasso.perf.cv[,4] , xgb.perf.cv[,5],pow=pow) 
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      cat ("\n lasso.minR - XGB with lasso offset          ") ;  nested.compare0(lasso.perf.cv[,4] , xgb.perf.cv[,6],pow=pow)   
    }
  }
  
  if ((dolasso == 1) & (dorf == 1)) {
    cat ("\n lasso.minR - Random Forest                  ") ;  nested.compare0(lasso.perf.cv[,4] , rf.perf.cv[,1],pow=pow) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - RF with lasso feature          ") ;  nested.compare0(lasso.perf.cv[,4] , rf.perf.cv[,2],pow=pow) 
    }
    if ( (sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian") ) {
      cat ("\n lasso.minR - RF with lasso offset           ") ;  nested.compare0(lasso.perf.cv[,4] , rf.perf.cv[,3],pow=pow)   
    }
  }
  
  if ((dolasso == 1) & (doorf == 1)) {
    cat ("\n lasso.minR - Oblique Random Forest          ") ;  nested.compare0(lasso.perf.cv[,4] , orf.perf.cv[,1],pow=pow) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - ORF with lasso feature         ") ;  nested.compare0(lasso.perf.cv[,4] , orf.perf.cv[,2],pow=pow) 
    }
    if ( (sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian") ) {
      cat ("\n lasso.minR - ORF with lasso offset          ") ;  nested.compare0(lasso.perf.cv[,4] , orf.perf.cv[,3],pow=pow)   
    }
  }
  
  if ((dolasso == 1) & (doann == 1)) {
    cat ("\n lasso.minR - ANN                            ") ;  nested.compare0(lasso.perf.cv[,4] , ann.perf.cv[,1],pow=pow) 
    if (ensemble[6]) { 
      cat ("\n lasso.minR - ANN l lasso feature            ") ;  nested.compare0(lasso.perf.cv[,4] , ann.perf.cv[,6],pow=pow)   
    } else if (ensemble[2]) { 
      cat ("\n lasso.minR - ANN lasso feature              ") ;  nested.compare0(lasso.perf.cv[,4] ,  ann.perf.cv[,2],pow=pow)   
    }  
    if (ensemble[8]) { 
      cat ("\n lasso.minR - ANN l lasso offset (upated)    ") ;  nested.compare0(lasso.perf.cv[,4] , ann.perf.cv[,8],pow=pow)   
    } else if (ensemble[4]) { 
      cat ("\n lasso.minR - ANN lasso offset               ") ;  nested.compare0(lasso.perf.cv[,4] , ann.perf.cv[,4],pow=pow)  
    } else if (ensemble[7]) { 
      cat ("\n lasso.minR - ANN l lasso offset (upated)    ") ;  nested.compare0(lasso.perf.cv[,4] , ann.perf.cv[,7],pow=pow)   
    }  else if (ensemble[3]) { 
      cat ("\n lasso.minR - ANN lasso offset               ") ;  nested.compare0(lasso.perf.cv[,4] , ann.perf.cv[,3],pow=pow)   
    }
  }
  
  if (dolasso) { cat("\n") } 
  
  if ((doxgb == 1) & (dorf == 1)) {
    cat ("\n XGBoost (tuned) - RF                        ") ;  nested.compare0(xgb.perf.cv[,4] ,  rf.perf.cv[,1],pow=pow)   
    if (sum(ensemble[c(2,6)]) > 0) {
      cat ("\n XGBoost lasso feature-RF with lasso feature ") ;  nested.compare0(xgb.perf.cv[,5] ,  rf.perf.cv[,2],pow=pow)   
    } 
    if ( (sum(ensemble[c(3,4,7,8)]) > 0)  & (family == "gaussian") ) {
      cat ("\n XGBoost lasso offset-  RF with lasso offset ") ;  nested.compare0(xgb.perf.cv[,6] ,  rf.perf.cv[,3],pow=pow)   
    }
  }
  
  if ((doxgb == 1) & (doorf == 1)) {
    cat ("\n XGBoost (tuned) - ORF                       ") ;  nested.compare0(xgb.perf.cv[,4] ,  orf.perf.cv[,1],pow=pow)   
    if (sum(ensemble[c(2,6)]) > 0) {
      cat ("\n XGBoost lasso feature-ORF with lasso feature") ;  nested.compare0(xgb.perf.cv[,5] ,  orf.perf.cv[,2],pow=pow)   
    } 
    if ( (sum(ensemble[c(3,4,7,8)]) > 0)  & (family == "gaussian") ) {
      cat ("\n XGBoost lasso offset-  ORF with lasso offset") ;  nested.compare0(xgb.perf.cv[,6] ,  orf.perf.cv[,3],pow=pow)   
    }
  }
  
  if ((doxgb == 1) & (doann == 1)) {
    cat ("\n XGBoost (tuned) - ANN                       ") ;  nested.compare0(xgb.perf.cv[,4] ,  ann.perf.cv[,1],pow=pow)   
    if (ensemble[6]) { 
      cat ("\n XGBoost lasso feature - ANN, l lasso feature ") ;  nested.compare0(xgb.perf.cv[,5] ,  ann.perf.cv[,6],pow=pow)   
    } else if (ensemble[2]) { 
      cat ("\n XGBoost lasso feature - ANN lasso feature   ") ;  nested.compare0(xgb.perf.cv[,5] ,  ann.perf.cv[,2],pow=pow)   
    } 
    if (family == "gaussian") {
      if (ensemble[8]) { 
        cat ("\n XGBoost lasso offset-ANN l lasso offset(upated) ") ;  nested.compare0(xgb.perf.cv[,6] ,  ann.perf.cv[,8],pow=pow)   
      } else if (ensemble[4]) { 
        cat ("\n XGBoost lasso offset - ANN, lasso offset   ") ;  nested.compare0(xgb.perf.cv[,6] ,  ann.perf.cv[,4],pow=pow)   
      } else if (ensemble[7]) { 
        cat ("\n XGBoost lasso offset - ANN l lasso offset  ") ;  nested.compare0(xgb.perf.cv[,6] ,  ann.perf.cv[,7],pow=pow)   
      } else if (ensemble[3]) { 
        cat ("\n XGBoost offset - ANN lasso offset          ") ;  nested.compare0(xgb.perf.cv[,6] ,  ann.perf.cv[,3],pow=pow)   
      }  
    }
  }
  
  if (doxgb) { cat("\n") }
  
  if ((dorf == 1) & (doorf == 1)) {
    cat ("\n RF - ORF                                    ") ;  nested.compare0(rf.perf.cv[,1] ,  orf.perf.cv[,1], pow=pow) 
    if (sum(ensemble[c(2,6)]) > 0) {
      cat ("\n RF lasso feature - ORF l lasso feature      " ) ;  nested.compare0(rf.perf.cv[,2] , orf.perf.cv[,2],pow=pow) 
    }  
    if (sum(ensemble[c(3,4,7,8)]) > 0) {
      cat ("\n RF lasso feature - ORF lasso feature        " ) ;  nested.compare0(rf.perf.cv[,3] , orf.perf.cv[,3],pow=pow)  
    }
    cat("\n")
  }
  
  if ((dorf == 1) & (doann == 1)) {
    cat ("\n RF - ANN                                    ") ;  nested.compare0(rf.perf.cv[,1] ,  ann.perf.cv[,1],pow=pow) 
    if (ensemble[6]) {
      cat ("\n RF lasso feature - ANN  l lasso feature     " ) ;  nested.compare0(rf.perf.cv[,2] ,  ann.perf.cv[,6],pow=pow) 
    } else if (ensemble[2]) {
      cat ("\n RF lasso feature - ANN lasso feature        " ) ;  nested.compare0(rf.perf.cv[,2] ,  ann.perf.cv[,2],pow=pow)  
    }
    if (ensemble[8]) {
      cat ("\n RF lasso offset - ANN l lasso offset (upated) " ) ;  nested.compare0(rf.perf.cv[,3] ,  ann.perf.cv[,8],pow=pow)   
    } else if (ensemble[4]) {
      cat ("\n RF lasso offset - ANN lasso offset           " ) ;  nested.compare0(rf.perf.cv[,3] ,  ann.perf.cv[,4],pow=pow)  
    } else if (ensemble[7]) {
      cat ("\n RF lasso offset - ANN, l lasso offset (upated) " ) ;  nested.compare0(rf.perf.cv[,3] ,  ann.perf.cv[,7],pow=pow) 
    } else if (ensemble[3]) {
      cat ("\n RF lasso offset - ANN, lasso offset          " ) ;  nested.compare0(rf.perf.cv[,3] ,  ann.perf.cv[,3],pow=pow)  
    }
    cat("\n")
  }
  
  if ((doorf == 1) & (doann == 1)) {
    cat ("\n ORF - ANN                                   ") ;  nested.compare0(orf.perf.cv[,1] ,  ann.perf.cv[,1],pow=pow) 
    if (ensemble[6]) {
      cat ("\n ORF lasso feature - ANN  l lasso feature    " ) ;  nested.compare0(orf.perf.cv[,2] ,  ann.perf.cv[,6],pow=pow) 
    } else if (ensemble[2]) {
      cat ("\n ORF lasso feature - ANN lasso feature       " ) ;  nested.compare0(orf.perf.cv[,2] ,  ann.perf.cv[,2],pow=pow)  
    }
    if (ensemble[8]) {
      cat ("\n ORF lasso offset - ANN l lasso offset (upated)" ) ;  nested.compare0(orf.perf.cv[,3] ,  ann.perf.cv[,8],pow=pow)   
    } else if (ensemble[4]) {
      cat ("\n ORF lasso offset - ANN lasso offset          " ) ;  nested.compare0(orf.perf.cv[,3] ,  ann.perf.cv[,4],pow=pow)  
    } else if (ensemble[7]) {
      cat ("\n ORF lasso offset - ANN, l lasso offset (upated)" ) ;  nested.compare0(orf.perf.cv[,3] ,  ann.perf.cv[,7],pow=pow) 
    } else if (ensemble[3]) {
      cat ("\n ORF lasso offset - ANN, lasso offset         " ) ;  nested.compare0(orf.perf.cv[,3] ,  ann.perf.cv[,3],pow=pow)  
    }
    cat("\n")
  }
  
  cat("\n")
}

###############################################################################################################
###############################################################################################################

#' A redirect to nested.compare
#'
#' @description 
#' See nested.compare(), as glmnetr() is depricated  
#'
#' @param object A nested.glmnetr output object.
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4 
#' @param pow the power to which the average of correlations is to be raised.   
#' @param type determines what type of nested cross validation performance measures are 
#' compared.  Possible values are "devrat" to compare the deviance ratios, i.e. the 
#' fractional reduction in deviance relative to the null model deviance, 
#' "agree" to compare agreement, "lincal" to compare the linear calibration 
#' slope coefficients, "intcal" to compare the linear calibration intercept 
#' coefficients, from the nested cross validation. 
#'  
#' @return A printout to the R console. 
#' 
#' @export
#' 

glmnetr.compcv = function(object, digits=4, type="devrat", pow=1) {
  nested.compare(object, digits, type, pow) 
}
  
###############################################################################################################
###############################################################################################################
