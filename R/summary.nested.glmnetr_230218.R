####################################################################################################################################
####################################################################################################################################
#' Summarize a a nested.glmnetr() output object
#'
#' @description 
#' Summarize the model fit from a nested.glmnetr() output object, i.e. the fit of 
#' a cross-validation informed relaxed lasso model fit, inferred by nested cross 
#' validation.  Else summarize the cross-validated model fit.    
#'
#' @param object a nested.glmnetr() output object.  
#' @param cvfit  default of FALSE to summarize fit of a cross validation informed 
#' relaxed lasso model fit, inferred by nested cross validation.  Option of TRUE 
#' will describe the cross validation informed relaxed lasso model itself. 
#' @param printg1 TRUE to also print out the fully penalized lasso beta, else to suppress.  
#' Only applies to cvfit=TRUE.
#' @param ... Additional arguments passed to the summary function.  
#' 
#' @return - a nested cross validation fit summary, or a cross validation model summary.  
#' 
#' @seealso
#'   \code{\link{glmnetr.compcv}} , \code{\link{summary.cv.stepreg}} , \code{\link{nested.glmnetr}} 
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
#' summary(fit3)
#' }
#' 
summary.nested.glmnetr = function(object, cvfit=FALSE, printg1=FALSE, ...) {
  
  if (cvfit==TRUE) {
    cv.glmnet.fit = object$cv.glmnet.fit  
    summary(cv.glmnet.fit,printg1=printg1)
  } else {
    
    tuning        = object$tuning      
#    names (object$sample)[1] = "family" 
    family        = object$sample[1] 
    if (tuning[4]==1) { 
      lassonzerocv  = object$lassonzerocv
      lassodeviancv = object$lassodeviancv
      lassolincalcv = object$lassolincalcv
      lassoagreecv  = object$lassoagreecv
      lasso.naive.agree = object$lasso.naive.agree
      lassoAveNZero  = colMeans(lassonzerocv)  
      lassoAveDevian = colMeans(lassodeviancv)  
      lassoAveLincal = colMeans(lassolincalcv)  
      lassoAveAgree  = colMeans(lassoagreecv) 
    } 
    if (tuning[5]==1) { 
      rpartnzerocv  = object$rpartnzerocv
      rpartdeviancv = object$rpartdeviancv
      rpartlincalcv = object$rpartlincalcv
      rpartagreecv  = object$rpartagreecv
      rpart.nzero   = object$rpart.nzero
      rpart.agree.naive = object$rpart.agree.naive
      rpartAveNZero  = colMeans(rpartnzerocv)  
      rpartAveDevian = colMeans(rpartdeviancv)  
      rpartAveLincal = colMeans(rpartlincalcv)  
      rpartAveAgree  = colMeans(rpartagreecv) 
    } 
    if ((tuning[6] == 1) | (tuning[7]==1)) {
      stepdeviancv  = object$stepdeviancv
      stepagreecv   = object$stepagreecv 
      step_dfcv     = object$step_dfcv
      step_pcv      = object$step_pcv
      StepAveDevian = colMeans( stepdeviancv)
      StepAveAgree  = colMeans( stepagreecv )
      StepAve_df    = colMeans( step_dfcv   )
      StepAve_p     = colMeans( step_pcv    )
    }
    if (tuning[6]==1) { func.fit.aic      = object$func.fit.aic      }
    if (tuning[7]==1) { cv.stepreg.fit    = object$cv.stepreg.fit    }  
  
    cat(paste0("\n"  , " Sample information including number of records, "))     
    if (family %in% c("cox","binomial")) { cat(paste0("events,")) }
    cat(paste0( "number of columns in", "\n", " design (predictor, X) matrix, and df (rank) of design matrix: ", "\n") )
    if (family %in% c("cox","binomial")) { print(object$sample) 
    } else { print(object$sample[-3]) }

#    cat(paste0("\n"  , " Number of columns in design (predictor, X) matrix, and df (rank) of design matrix: ", "\n") )
#    print( object$sample[3:4] ) 
    
    if ((tuning[6]==0) & (tuning[7]==0)) { 
      cat(paste0("\n"  , " Tuning parameters for lasso/rpart model: ", "\n") )
      print(object$tuning[c(2,8)]) 
    } else { 
      cat(paste0("\n"  , " Tuning parameters for models: ", "\n") )
      print(object$tuning) 
    }
    
#    if (tuning[7]==1) {
#      cat(paste0("\n", " Average deviance for null model", "\n") )    ## pull other data applicable to all models 
#      print( round(StepAveDevian[1], digits = 4 ) ) 
#    }
    
    if (family == "cox") { perunit = "deviance per event " } else { perunit = "deviance per record " }
    
    if (tuning[4] == 1) {
      cat(paste0("\n" , " Nested Cross Validation averages for LASSO (1se and min), Relaxed LASSO, and gamma=0 LASSO : ", "\n") )
      cat(paste0("\n" , "      ", perunit, ": ", "\n") )
      print( round( lassoAveDevian , digits = 4) )
      
      cat(paste0("\n" , "      number of nonzero model terms : ", "\n") )
      print( round( lassoAveNZero , digits = 2) )    
      
      cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
      print( round( lassoAveLincal , digits = 4) ) 
      
      # In summary.nested.glmnetr state only "concordance" or "R-square" in output
      
      if (family %in% c("cox","binomial")) { 
        cat(paste0("\n", "      agreement (concordance):         ", "\n") )
      } else { 
        cat(paste0("\n", "      agreement (r-square):            ", "\n") )
      } 
      
      print( round( lassoAveAgree , digits = 4) ) 
      
      cat(paste0("\n", " Naive agreement for cross validation informed lasso model : ",  "\n") )
      names(lasso.naive.agree) = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0") 
      print( round( lasso.naive.agree , digits=4) )
    }                                                                             ## how to get rid of [1,]  ??

    if (tuning[5] == 1) { 
      cat(paste0("\n" , " Nested Cross Validation averages for RPART : ", "\n") )

      cat(paste0("\n", "    average number of terms used in cv informed models : ", "\n") )
      names(rpartAveNZero) = c("cp=0.00", "cp=0.01", "cp=0.02") 
      print(round( rpartAveNZero, digits=1)) 
      
      if (family %in% c("cox","binomial")) { 
        cat(paste0("\n", "    average agreement (concordance):         ", "\n") )
      } else { 
        cat(paste0("\n", "    agreement (r-square):            ", "\n") )
      } 
      
      names(rpartAveAgree) = c("cp=0.00", "cp=0.01", "cp=0.02") 
      print( round( rpartAveAgree , digits = 4) ) 
      
      cat(paste0("\n", " Cross validation informed RPART model : ",  "\n") )
      
      cat(paste0("\n", "   number of terms used in model : ", "\n") )
      names(rpart.nzero) = c("cp=0.00", "cp=0.01", "cp=0.02") 
      print(round( rpart.nzero, digits=1)) 
      
      cat(paste0("\n", "   naive agreement : ",  "\n") )
      names(rpart.agree.naive) = c("cp=0.00", "cp=0.01", "cp=0.02") 
      print( round( rpart.agree.naive , digits=4) )
    }                                                                             ## how to get rid of [1,]  ??
        
    if (tuning[7]==1) {
#      cv.stepreg.fit = object$cv.stepreg.fit
#      cv.stepreg.fit.df = cv.stepreg.fit$best.df    #func.fit.df
      cat(paste0("\n", " Nested Cross Validation stepwise regression model (df): ", "\n") )   
      cat(paste0("      Average deviance : ", round(StepAveDevian[1]   ,  digits = 4), "\n") )    
      cat(paste0("      Average model df : ", round( StepAve_df[1], digits=2 ), "\n") )
      if (family %in% c("cox","binomial")) {
        cat(paste0("      Concordance      : ", round(StepAveAgree[1],  digits = 4), "\n") )        
        cat(paste0(" Naive concordance based upon the same (all) data as model derivation (df): ",
                   round(cv.stepreg.fit$cvfit.df[6], digits=4) , "\n") )
      } else {
        cat(paste0("      R-square         : ", round(StepAveAgree[1],  digits = 4), "\n") )        
        cat(paste0(" Naive R-square based upon the same (all) data as model derivation (df): ",
                    round(cv.stepreg.fit$cvfit.df[6], digits=4) , "\n") )
      }
      
      cat(paste0("\n", " Nested Cross Validation stepwise regression model (p): ", "\n") )   
      cat(paste0("      Average deviance : ", round(StepAveDevian[2]   ,  digits = 4), "\n") )    
      cat(paste0("      Average model p  : ", round( StepAve_p [1], digits=3 ), "\n") )
      cat(paste0("      Average model df : ", round( StepAve_df[2], digits=2 ), "\n") )
      if (family %in% c("cox","binomial")) {
        cat(paste0("      Concordance      : ", round(StepAveAgree[2],  digits = 4), "\n") )        
        cat(paste0(" Naive concordance based upon the same (all) data as model derivation (df): ",
                   round(cv.stepreg.fit$cvfit.p[6], digits=4) , "\n") )
      } else {
        cat(paste0("      R-square         : ", round(StepAveAgree[2],  digits = 4), "\n") )        
        cat(paste0(" Naive R-square based upon the same (all) data as model derivation (p): ",
                   round(cv.stepreg.fit$cvfit.p[6], digits=4) , "\n") )
      }
    }
    
    if (tuning[6]==1) {
      cat(paste0("\n", " Cross Validation results for stepwise regression model: (AIC)", "\n") )    
      cat(paste0("      Average deviance : ", round(StepAveDevian[3],  digits = 4), "\n") )    
      cat(paste0("      Average model df : ", round( StepAve_df[3], digits=2 ), "\n") )
      cat(paste0("      Concordance      : ", round( StepAveAgree[3], digits = 4), "\n") )   
      if (family %in% c("cox","binomial")) {
        cat(paste0(" Naive concordance based upon the same (all) data as model derivation (AIC) : ", 
                   round(object$func.fit.aic$aic.fit.n[6], digits=4) , "\n") )
      } else {
        cat(paste0(" Naive R-square based upon the same (all) data as model derivation (AIC) : ",         
                   round( object$func.fit.aic$aic.fit.n[6], digits=4 ), "\n") ) 
        #        round( 1-var(all.fit.aic$residuals)/var(all.fit.aic$y), digits=4), "\n") ) 
      }
    }
    
    cat("\n")  
  }
} 

###############################################################################################################
###############################################################################################################

#' Give predicteds based upon the cv.glmnet output object contained in the nested.glmnetr output object.
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
#'   \code{\link{predict.cv.glmnetr}}  
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
  if (class(object) %in% c("nested.glmnetr","nested.stepreg")) { object = object$cv.glmnet.fit.f }
  predict.cv.glmnetr( object=object, xs_new=xs_new, lam=lam, gam=gam, comment=comment) 
}

################################################################################
################################################################################

#' A glmnetr specifc paired t-test
#'
#' @description 
#' Perform a paired t-test as called from glmnetr.compcv().  
#'
#' @param a One term
#' @param b A second term 
#'
#' @return  A t-test
#' 
#' @importFrom stats t.test 
#'
glmnetr.compcv0 = function(a,b) { 
  tdiff = t.test(a-b)
  cat ( paste0(  " estimate (95% CI): ", round(tdiff$estimate, digits=4), " (", round(tdiff$conf.int[1],digits=4), ", ", 
                 round(tdiff$conf.int[2],digits=4), ") , p=", round(tdiff$p.value,digits=4) ) )
}

################################################################################
################################################################################

#' Compare cross validation fits from a nested.glmnetr output.
#'
#' @description 
#' Compare cross-validation model fits in terms of average concordance from the 
#' nested cross validaiton fits.  
#'
#' @param object A nested.glmnetr output object.
#'
#' @return A printout to the R console. 
#' 
#' @seealso
#'   \code{\link{summary.nested.glmnetr}} 
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
#' glmnetr.compcv(fit3)
#' }
#' 
glmnetr.compcv = function(object) {
  tuning        = object$tuning
  dolasso = tuning[4]
  dorpart = tuning[5]
  doaic   = tuning[6]
  dostep  = tuning[7]
  lassoagreecv  = object$lassoagreecv
  stepagreecv   = object$stepagreecv
  if (dolasso == 1) {
    cat ("\n", "C lasso.min   - lasso.minR , ") ;  glmnetr.compcv0(lassoagreecv[,2] , lassoagreecv[,4]) 
    cat ("\n", "C lasso.min   - lasso.minR0, ") ;  glmnetr.compcv0(lassoagreecv[,2] , lassoagreecv[,6]) ;  
    cat ("\n", "C lasso.minR  - lasso.minR0, ") ;  glmnetr.compcv0(lassoagreecv[,4] , lassoagreecv[,6]) ;  cat("\n")
  }
  if ((dolasso == 1) & (dostep == 1)) {
    cat ("\n", "C lasso.min   - step (df),   ") ;  glmnetr.compcv0(lassoagreecv[,2] , stepagreecv[,1]) 
    cat ("\n", "C lasso.minR  - step (df),   ") ;  glmnetr.compcv0(lassoagreecv[,4] , stepagreecv[,1]) 
    cat ("\n", "C lasso.minR0 - step (df),   ") ;  glmnetr.compcv0(lassoagreecv[,6] , stepagreecv[,1])     ;  cat("\n")
    
    cat ("\n", "C lasso.min   - step (p),   ") ;  glmnetr.compcv0(lassoagreecv[,2] , stepagreecv[,2]) 
    cat ("\n", "C lasso.minR  - step (p),   ") ;  glmnetr.compcv0(lassoagreecv[,4] , stepagreecv[,2]) 
    cat ("\n", "C lasso.minR0 - step (p),   ") ;  glmnetr.compcv0(lassoagreecv[,6] , stepagreecv[,2])    ;  cat("\n")
    
    cat ("\n", "C step (df) - step (p),     ") ;  glmnetr.compcv0(stepagreecv[,1]      , stepagreecv[,2])    ;  cat("\n")
  }
  
  if ((dolasso == 1) & (doaic == 1)) {
    cat ("\n", "C lasso.min   - step (AIC),  ") ;  glmnetr.compcv0(lassoagreecv[,2] , stepagreecv[,3]) 
    cat ("\n", "C relax.minR  - step (AIC),  ") ;  glmnetr.compcv0(lassoagreecv[,4] , stepagreecv[,3]) 
    cat ("\n", "C lasso.minR0 - step (AIC),  ") ;  glmnetr.compcv0(lassoagreecv[,6] , stepagreecv[,3]) ;  cat("\n")    
  }
  
  cat("\n")
}

####################################################################################################################################
####################################################################################################################################

