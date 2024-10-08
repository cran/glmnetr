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

nested.compare0 = function(a, b, pow=1, bootstrap=0, digits=4, txt=0) { 
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
  if (txt==1) {
    cat ( paste0(  " estimate (95% CI): ", round(meandiff, digits=digits), " (", round(lo, digits=digits), ", ", 
                   round(up, digits=digits), ") , p=", round(p_, digits=digits) ) )
  } else {
    cat ( paste0( round(meandiff, digits=digits), " (", round(lo, digits=digits), ", ", 
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
nested.compare = function(object, type="devrat", digits=4, pow=1) {
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
  
  if (dostep == 1) {
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
    doxgb = 0 ; dorf = 0 ; doorf = 0 ; doann = 0 ; dorpart = 0 ; dostep = 0 ; doaic = 0 ; 
  }
  
  cat ("  Model performance comparison in terms of", pm, "\n\n" )   
  cat ("  Comparison                                estimate   (95% CI)         p\n") 
  
  if (dolasso == 1) {
    cat ("\n lasso.minR  - lasso.min                     ") ;  nested.compare0(lasso.perf.rep[,4] , lasso.perf.rep[,2], pow, bootstrap, digits) 
    cat ("\n lasso.minR  - lasso.minR0                   ") ;  nested.compare0(lasso.perf.rep[,4] , lasso.perf.rep[,6], pow, bootstrap, digits)   
    cat ("\n lasso.min   - lasso.minR0                   ") ;  nested.compare0(lasso.perf.rep[,2] , lasso.perf.rep[,6], pow, bootstrap, digits)   
    cat("\n")
  }
  
  #  print(xgb.perf.rep)
  
  if (doxgb == 1) {
    cat ("\n XGBoost (tuned) - XGBoost (simple)          ") ;  nested.compare0(xgb.perf.rep[,4] , xgb.perf.rep[,1], pow, bootstrap, digits)  
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n XGBoost (tuned) lasso feature - no feature  ") ;  nested.compare0(xgb.perf.rep[,5] , xgb.perf.rep[,4], pow, bootstrap, digits)   
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      cat ("\n XGBoost (tuned) lasso offset - no offset    ") ;  nested.compare0(xgb.perf.rep[,6] , xgb.perf.rep[,4], pow, bootstrap, digits)   
    }
    cat("\n")
  }
  
  if (dorf == 1) {
    lr = 0 
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n RF with lasso feature - no feature          ") ;  nested.compare0(rf.perf.rep[,2] , rf.perf.rep[,1], pow, bootstrap, digits)   
      lr = 1 
    }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) {
      cat ("\n RF with lasso offset - no offset            ") ;  nested.compare0(rf.perf.rep[,3] , rf.perf.rep[,1], pow, bootstrap, digits)   
      lr = 1 
    }
    if (lr == 1) { cat("\n") } 
  }
  
  if (doorf == 1) {
    lr = 0 
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n ORF with lasso feature - no feature          ") ;  nested.compare0(orf.perf.rep[,2] , orf.perf.rep[,1], pow, bootstrap, digits)   
      lr = 1 
    }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) {
      cat ("\n ORF with lasso offset - no offset            ") ;  nested.compare0(orf.perf.rep[,3] , orf.perf.rep[,1], pow, bootstrap, digits)   
      lr = 1 
    }
    if (lr == 1) { cat("\n") } 
  }
  
  if (doann == 1) {
    lr = 0 
    if (sum(ensemble[6])> 0) {
      cat ("\n ANN with with lasso feature - no feature    ") ;  nested.compare0(ann.perf.rep[,6] , ann.perf.rep[,1], pow, bootstrap, digits) ; lr = 1 ;
    } else if (sum(ensemble[2])> 0) {
      cat ("\n ANN with with lasso feature - no feature    ") ;  nested.compare0(ann.perf.rep[,2] , ann.perf.rep[,1], pow, bootstrap, digits) ; lr = 1 ;
      
    } 
    if (sum(ensemble[8])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  nested.compare0(ann.perf.rep[,8] , ann.perf.rep[,1], pow, bootstrap, digits) ; lr = 1 ;
    } else if (sum(ensemble[7])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  nested.compare0(ann.perf.rep[,7] , ann.perf.rep[,1], pow, bootstrap, digits) ; lr = 1 ;
    } else     if (sum(ensemble[4])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  nested.compare0(ann.perf.rep[,4] , ann.perf.rep[,1], pow, bootstrap, digits) ; lr = 1 ;
    } else     if (sum(ensemble[3])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  nested.compare0(ann.perf.rep[,3] , ann.perf.rep[,1], pow, bootstrap, digits) ; lr = 1 ;
    } 
    if (lr == 1) { cat("\n") } 
  }
  
  if (dostep == 1) {
    cat ("\n step (df) - step (p)                        ") ;  nested.compare0(step.perf.rep[,1]      , step.perf.rep[,2], pow, bootstrap, digits) ;  cat("\n")
  }
  
  #  cat("\n")
  
  if ((dolasso == 1) & (doxgb == 1)) {
    cat ("\n lasso.minR - XGB (tuned)                    ") ;  nested.compare0(lasso.perf.rep[,4] , xgb.perf.rep[,4], pow, bootstrap, digits) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - XGB with lasso feature         ") ;  nested.compare0(lasso.perf.rep[,4] , xgb.perf.rep[,5], pow, bootstrap, digits) 
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      cat ("\n lasso.minR - XGB with lasso offset          ") ;  nested.compare0(lasso.perf.rep[,4] , xgb.perf.rep[,6], pow, bootstrap, digits)   
    }
  }
  
  if ((dolasso == 1) & (dorf == 1)) {
    cat ("\n lasso.minR - Random Forest                  ") ;  nested.compare0(lasso.perf.rep[,4] , rf.perf.rep[,1], pow, bootstrap, digits) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - RF with lasso feature          ") ;  nested.compare0(lasso.perf.rep[,4] , rf.perf.rep[,2], pow, bootstrap, digits) 
    }
    if ( (sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian") ) {
      cat ("\n lasso.minR - RF with lasso offset           ") ;  nested.compare0(lasso.perf.rep[,4] , rf.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if ((dolasso == 1) & (doorf == 1)) {
    cat ("\n lasso.minR - Oblique Random Forest          ") ;  nested.compare0(lasso.perf.rep[,4] , orf.perf.rep[,1], pow, bootstrap, digits) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - ORF with lasso feature         ") ;  nested.compare0(lasso.perf.rep[,4] , orf.perf.rep[,2], pow, bootstrap, digits) 
    }
    if ( (sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian") ) {
      cat ("\n lasso.minR - ORF with lasso offset          ") ;  nested.compare0(lasso.perf.rep[,4] , orf.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if ((dolasso == 1) & (doann == 1)) {
    cat ("\n lasso.minR - ANN                            ") ;  nested.compare0(lasso.perf.rep[,4] , ann.perf.rep[,1], pow, bootstrap, digits) 
    if (ensemble[6]) { 
      cat ("\n lasso.minR - ANN l lasso feature            ") ;  nested.compare0(lasso.perf.rep[,4] , ann.perf.rep[,6], pow, bootstrap, digits)   
    } else if (ensemble[2]) { 
      cat ("\n lasso.minR - ANN lasso feature              ") ;  nested.compare0(lasso.perf.rep[,4] ,  ann.perf.rep[,2], pow, bootstrap, digits)   
    }  
    if (ensemble[8]) { 
      cat ("\n lasso.minR - ANN l lasso offset (upated)    ") ;  nested.compare0(lasso.perf.rep[,4] , ann.perf.rep[,8], pow, bootstrap, digits)   
    } else if (ensemble[4]) { 
      cat ("\n lasso.minR - ANN lasso offset               ") ;  nested.compare0(lasso.perf.rep[,4] , ann.perf.rep[,4], pow, bootstrap, digits)  
    } else if (ensemble[7]) { 
      cat ("\n lasso.minR - ANN l lasso offset (upated)    ") ;  nested.compare0(lasso.perf.rep[,4] , ann.perf.rep[,7], pow, bootstrap, digits)   
    }  else if (ensemble[3]) { 
      cat ("\n lasso.minR - ANN lasso offset               ") ;  nested.compare0(lasso.perf.rep[,4] , ann.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if (dolasso) { cat("\n") } 
  
  if ((doxgb == 1) & (dorf == 1)) {
    cat ("\n XGBoost (tuned) - RF                        ") ;  nested.compare0(xgb.perf.rep[,4] ,  rf.perf.rep[,1], pow, bootstrap, digits)   
    if (sum(ensemble[c(2,6)]) > 0) {
      cat ("\n XGBoost lasso feature-RF with lasso feature ") ;  nested.compare0(xgb.perf.rep[,5] ,  rf.perf.rep[,2], pow, bootstrap, digits)   
    } 
    if ( (sum(ensemble[c(3,4,7,8)]) > 0)  & (family == "gaussian") ) {
      cat ("\n XGBoost lasso offset-  RF with lasso offset ") ;  nested.compare0(xgb.perf.rep[,6] ,  rf.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if ((doxgb == 1) & (doorf == 1)) {
    cat ("\n XGBoost (tuned) - ORF                       ") ;  nested.compare0(xgb.perf.rep[,4] ,  orf.perf.rep[,1], pow, bootstrap, digits)   
    if (sum(ensemble[c(2,6)]) > 0) {
      cat ("\n XGBoost lasso feature-ORF with lasso feature") ;  nested.compare0(xgb.perf.rep[,5] ,  orf.perf.rep[,2], pow, bootstrap, digits)   
    } 
    if ( (sum(ensemble[c(3,4,7,8)]) > 0)  & (family == "gaussian") ) {
      cat ("\n XGBoost lasso offset-  ORF with lasso offset") ;  nested.compare0(xgb.perf.rep[,6] ,  orf.perf.rep[,3], pow, bootstrap, digits)   
    }
  }
  
  if ((doxgb == 1) & (doann == 1)) {
    cat ("\n XGBoost (tuned) - ANN                       ") ;  nested.compare0(xgb.perf.rep[,4] ,  ann.perf.rep[,1], pow, bootstrap, digits)   
    if (ensemble[6]) { 
      cat ("\n XGBoost lasso feature - ANN, l lasso feature ") ;  nested.compare0(xgb.perf.rep[,5] ,  ann.perf.rep[,6], pow, bootstrap, digits)   
    } else if (ensemble[2]) { 
      cat ("\n XGBoost lasso feature - ANN lasso feature   ") ;  nested.compare0(xgb.perf.rep[,5] ,  ann.perf.rep[,2], pow, bootstrap, digits)   
    } 
    if (family == "gaussian") {
      if (ensemble[8]) { 
        cat ("\n XGBoost lasso offset-ANN l lasso offset(upated) ") ;  nested.compare0(xgb.perf.rep[,6] ,  ann.perf.rep[,8], pow, bootstrap, digits)   
      } else if (ensemble[4]) { 
        cat ("\n XGBoost lasso offset - ANN, lasso offset   ") ;  nested.compare0(xgb.perf.rep[,6] ,  ann.perf.rep[,4], pow, bootstrap, digits)   
      } else if (ensemble[7]) { 
        cat ("\n XGBoost lasso offset - ANN l lasso offset  ") ;  nested.compare0(xgb.perf.rep[,6] ,  ann.perf.rep[,7], pow, bootstrap, digits)   
      } else if (ensemble[3]) { 
        cat ("\n XGBoost offset - ANN lasso offset          ") ;  nested.compare0(xgb.perf.rep[,6] ,  ann.perf.rep[,3], pow, bootstrap, digits)   
      }  
    }
  }
  
  if (doxgb) { cat("\n") }
  
  if ((dorf == 1) & (doorf == 1)) {
    cat ("\n RF - ORF                                    ") ;  nested.compare0(rf.perf.rep[,1] ,  orf.perf.rep[,1], pow, bootstrap, digits) 
    if (sum(ensemble[c(2,6)]) > 0) {
      cat ("\n RF lasso feature - ORF l lasso feature      " ) ;  nested.compare0(rf.perf.rep[,2] , orf.perf.rep[,2], pow, bootstrap, digits) 
    }  
    if (sum(ensemble[c(3,4,7,8)]) > 0) {
      cat ("\n RF lasso feature - ORF lasso feature        " ) ;  nested.compare0(rf.perf.rep[,3] , orf.perf.rep[,3], pow, bootstrap, digits)  
    }
    cat("\n")
  }
  
  if ((dorf == 1) & (doann == 1)) {
    cat ("\n RF - ANN                                    ") ;  nested.compare0(rf.perf.rep[,1] ,  ann.perf.rep[,1], pow, bootstrap, digits) 
    if (ensemble[6]) {
      cat ("\n RF lasso feature - ANN  l lasso feature     " ) ;  nested.compare0(rf.perf.rep[,2] ,  ann.perf.rep[,6], pow, bootstrap, digits) 
    } else if (ensemble[2]) {
      cat ("\n RF lasso feature - ANN lasso feature        " ) ;  nested.compare0(rf.perf.rep[,2] ,  ann.perf.rep[,2], pow, bootstrap, digits)  
    }
    if (ensemble[8]) {
      cat ("\n RF lasso offset - ANN l lasso offset (upated) " ) ;  nested.compare0(rf.perf.rep[,3] ,  ann.perf.rep[,8], pow, bootstrap, digits)   
    } else if (ensemble[4]) {
      cat ("\n RF lasso offset - ANN lasso offset           " ) ;  nested.compare0(rf.perf.rep[,3] ,  ann.perf.rep[,4], pow, bootstrap, digits)  
    } else if (ensemble[7]) {
      cat ("\n RF lasso offset - ANN, l lasso offset (upated) " ) ;  nested.compare0(rf.perf.rep[,3] ,  ann.perf.rep[,7], pow, bootstrap, digits) 
    } else if (ensemble[3]) {
      cat ("\n RF lasso offset - ANN, lasso offset          " ) ;  nested.compare0(rf.perf.rep[,3] ,  ann.perf.rep[,3], pow, bootstrap, digits)  
    }
    cat("\n")
  }
  
  if ((doorf == 1) & (doann == 1)) {
    cat ("\n ORF - ANN                                   ") ;  nested.compare0(orf.perf.rep[,1] ,  ann.perf.rep[,1], pow, bootstrap, digits) 
    if (ensemble[6]) {
      cat ("\n ORF lasso feature - ANN  l lasso feature    " ) ;  nested.compare0(orf.perf.rep[,2] ,  ann.perf.rep[,6], pow, bootstrap, digits) 
    } else if (ensemble[2]) {
      cat ("\n ORF lasso feature - ANN lasso feature       " ) ;  nested.compare0(orf.perf.rep[,2] ,  ann.perf.rep[,2], pow, bootstrap, digits)  
    }
    if (ensemble[8]) {
      cat ("\n ORF lasso offset - ANN l lasso offset (upated)" ) ;  nested.compare0(orf.perf.rep[,3] ,  ann.perf.rep[,8], pow, bootstrap, digits)   
    } else if (ensemble[4]) {
      cat ("\n ORF lasso offset - ANN lasso offset          " ) ;  nested.compare0(orf.perf.rep[,3] ,  ann.perf.rep[,4], pow, bootstrap, digits)  
    } else if (ensemble[7]) {
      cat ("\n ORF lasso offset - ANN, l lasso offset (upated)" ) ;  nested.compare0(orf.perf.rep[,3] ,  ann.perf.rep[,7], pow, bootstrap, digits) 
    } else if (ensemble[3]) {
      cat ("\n ORF lasso offset - ANN, lasso offset         " ) ;  nested.compare0(orf.perf.rep[,3] ,  ann.perf.rep[,3], pow, bootstrap, digits)  
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
#' See nested.compare(), as glmnetr.compcv() is depricated  
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
#' @seealso
#'   \code{\link{nested.compare}}
#' 
#' @export
#' 

glmnetr.compcv = function(object, digits=4, type="devrat", pow=1) {
  nested.compare(object, type, digits, pow) 
}
  
###############################################################################################################
###############################################################################################################
