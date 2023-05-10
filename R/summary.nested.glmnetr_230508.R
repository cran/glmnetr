####################################################################################################################################
####################################################################################################################################
#' Print an abbreviated summary of a nested.glmnetr() output object
#'
#' @param x a nested.glmnetr() output object.  
#' @param ... additional inputs for the print function.  This is not used here.   
#' 
#' @return - a nested cross validation fit summary, or a cross validation model summary.  
#' 
#' @seealso
#'    \code{\link{nested.glmnetr}}, \code{\link{glmnetr.compcv}} , \code{\link{summary.nested.glmnetr}} 
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
#' print(fit3)
#' }
#' 
print.nested.glmnetr = function(x, ...) {
  summary(x,short=1) 
}
#
################################################################################  
#' Summarize a nested.glmnetr() output object
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
#' @param short optionally print just the CV agreement summaries (short=1)
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
summary.nested.glmnetr = function(object, cvfit=FALSE, printg1=FALSE, short=0, ...) {
  
  if (cvfit==TRUE) {
    cv_glmnet_fit = object$cv_glmnet_fit  
    summary(cv_glmnet_fit,printg1=printg1)
  } else {
    tuning  = object$tuning
#    tuning  = tuning[-c(10:17)]
    dolasso = tuning[4]
    doann   = tuning[5]
    doxgb   = tuning[6]
    dorpart = tuning[7]
    dostep  = tuning[8]
    doaic   = tuning[9]
    doann_  = (doann[1] > 0)*1   
    ensemble = object$ensemble
    do_ncv  = object$do_ncv
#    temp_ = c(1,3,5,2,4)
#    ensemble = rep(0,6)
#    ensemble = temp_
    ensemble2 = ensemble[c(1,4,2,3,5,6)]  ## from complexity order to ensemble input order    
    ensemble 
    ensemble2 
    
    family = object$sample[1]
    sample = object$sample 
    sample[6] = round(as.numeric(sample[6]), digits=3)
    if (family=="cox") { names(sample)[6]="null.dev/events" 
    } else { names(sample)[6] = "null.dev/obs"  }

    if (dolasso==1) { 
      lasso.nzero.cv  = object$lasso.nzero.cv
      lasso.devian.cv = object$lasso.devian.cv
      lasso.cal.devian.cv = object$lasso.cal.devian.cv
      lasso.lincal.cv = object$lasso.lincal.cv
      lasso.agree.cv  = object$lasso.agree.cv
      lasso.agree.naive = object$lasso.agree.naive
      lassoAveNZero  = colMeans(lasso.nzero.cv)  
      lassoAveDevian = colMeans(lasso.devian.cv)  
      lassoCalAveDevian = colMeans(lasso.cal.devian.cv)  
      lassoAveLincal = colMeans(lasso.lincal.cv)  
      lassoAveAgree  = colMeans(lasso.agree.cv) 
    } 
    if (doxgb==1) { 
      xgb.devian.cv   = object$xgb.devian.cv
      xgb.lincal.cv   = object$xgb.lincal.cv
      xgb.agree.cv    = object$xgb.agree.cv
      xgb.agree.naive = object$xgb.agree.naive 
      xgbAveDevian = colMeans(xgb.devian.cv)   
      xgbAveLincal = colMeans(xgb.lincal.cv)   
      xgbAveAgree  = colMeans(xgb.agree.cv)    
    } 
    if (doann_==1) { 
      nms = c("Uninformed", "lasso terms", "lasso feat.", "init w's", "update w's", "offset") 
      ann.devian.cv = object$ann.devian.cv
      ann.lincal.cv = object$ann.lincal.cv
      ann.agree.cv  = object$ann.agree.cv 
      ann.agree.naive = object$ann.agree.naive[(ensemble2==1)]
      annAveDevian = colMeans(ann.devian.cv)[(ensemble2==1)]  
      annAveLincal = colMeans(ann.lincal.cv)[(ensemble2==1)]
      annAveAgree  = colMeans(ann.agree.cv )[(ensemble2==1)]
    } 
    if (dorpart==1) { 
      rpart.nzero.cv  = object$rpart.nzero.cv
      rpart.devian.cv = object$rpart.devian.cv
      rpart.lincal.cv = object$rpart.lincal.cv
      rpart.agree.cv  = object$rpart.agree.cv
      rpart.nzero     = object$rpart.nzero         
      rpart.agree.naive = object$rpart.agree.naive 
      rpartAveNZero  = colMeans(rpart.nzero.cv[,c(3,2,1,6,5,4,9,8,7)])    
      rpartAveDevian = colMeans(rpart.devian.cv[,c(3,2,1,6,5,4,9,8,7)])   
      rpartAveLincal = colMeans(rpart.lincal.cv[,c(3,2,1,6,5,4,9,8,7)])   
      rpartAveAgree  = colMeans(rpart.agree.cv[,c(3,2,1,6,5,4,9,8,7)])   
      rpartAveNZero  = rpartAveNZero [ c(1:6) ]
      rpartAveDevian = rpartAveDevian[ c(1:6) ] 
      rpartAveLincal = rpartAveLincal[ c(1:6) ]
      rpartAveAgree = rpartAveAgree  [ c(1:6) ]
      rpartAveNZero = rpartAveNZero  [ c(1:6) ]
    } 
    if ((doaic == 1) | (dostep==1)) {
      step.devian.cv  = object$step.devian.cv
      step.lincal.cv  = object$step.lincal.cv
      step.agree.cv   = object$step.agree.cv 
      step_df_cv    = object$step_df_cv
      step_p_cv      = object$step_p_cv
      StepAveDevian = colMeans( step.devian.cv)
      StepAveLincal = colMeans( step.lincal.cv)
      StepAveAgree  = colMeans( step.agree.cv )
      StepAve_df    = colMeans( step_df_cv   )
      StepAve_p     = colMeans( step_p_cv    )
    }
    if (doann == 1) { 
      if        (ensemble2[4]==1) { ann_cv = object$ann_fit_4 ; whichann = nms[4] ;
      } else if (ensemble2[5]==1) { ann_cv = object$ann_fit_5 ; whichann = nms[5] ;
      } else if (ensemble2[3]==1) { ann_cv = object$ann_fit_3 ; whichann = nms[3] ;
      } else if (ensemble2[2]==1) { ann_cv = object$ann_fit_2 ; whichann = nms[2] ;
      } else if (ensemble2[1]==1) { ann_cv = object$ann_fit_1 ; whichann = nms[1] ;
      }
    }
    if (dostep==1) { cv.stepreg.fit    = object$cv.stepreg.fit }  
    if (doaic ==1) { func.fit.aic      = object$func.fit.aic   }
  
    cat(paste0("\n"  , " Sample information including number of records, "))     
    if (family %in% c("cox","binomial")) { cat(paste0("events, ")) }
    cat(paste0( "number of columns in", "\n", " design (predictor, X) matrix, and df (rank) of design matrix: ", "\n") )
    if (family %in% c("cox","binomial")) { print(sample) 
    } else { print(sample[-3]) }

    if ((doaic==0) & (dostep==0)) { 
      cat(paste0("\n"  , " Tuning parameters for models : ", "\n") )
      if (doann_==1) {
        print(object$tuning[c(2,10,11)]) 
      } else {
        print(object$tuning[c(2,10)]) 
      }
    } else { 
      cat(paste0("\n"  , " Tuning parameters for models : ", "\n") )
      if (doann_==1) {
        print(object$tuning[c(1,2,10,11)]) 
      } else {
        print(object$tuning[c(1,2,10)]) 
      }
    }
    
    if (doann_ == 1) {
      cat(paste0("\n"  , " Tuning parameters for  ", whichann, "  ANN model : ", "\n") )
      print(ann_cv$modelsum[c(1:9)]) 
    }
    
#    if (doaic==1) {
#      cat(paste0("\n", " Average deviance for null model", "\n") )    ## pull other data applicable to all models 
#      print( round(StepAveDevian[1], digits = 4 ) ) 
#    }
    
    if (family == "cox") { perunit = "deviance per event " } else { perunit = "deviance per record " }
    
    ## LASSO ###################################################################
    
    if (dolasso == 1) {
      if (do_ncv == 1) {
        if (short==1) {
          if (family %in% c("cox","binomial")) { 
            cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed LASSO : ",  "\n") )
          } else { 
            cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed LASSO : ",  "\n") )
         } 
          print( round( lassoAveAgree , digits = 4) ) 
        } else {
          cat(paste0("\n\n" , " Nested Cross Validation averages for LASSO (1se and min), Relaxed LASSO, and gamma=0 LASSO : ", "\n") )
          cat(paste0("\n" , "      ", perunit, ": ", "\n") )
          print( round( lassoAveDevian , digits = 4) )
      
#          cat(paste0("\n\n" , " Nested Cross Validation averages for Linear Calibrated LASSO : ", "\n") )
          cat(paste0("\n" , "      ", perunit, "(linerly calibrated) : ", "\n") )
          print( round( lassoCalAveDevian , digits = 4) )
        
          cat(paste0("\n" , "      number of nonzero model terms : ", "\n") )
          print( round( lassoAveNZero , digits = 2) )    
      
          cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
          print( round( lassoAveLincal , digits = 4) ) 
      
          # In summary.nested.glmnetr state only "concordance" or "R-square" in output
      
          if (family %in% c("cox","binomial")) { 
            cat(paste0("\n", "      agreement (concordance) :         ", "\n") )
          } else { 
            cat(paste0("\n", "      agreement (r-square) :            ", "\n") )
          } 
      
          print( round( lassoAveAgree , digits = 4) ) 
        }
      }
      if ((short != 1) | (do_ncv == 0)) {
        cat(paste0("\n", " Naive agreement for cross validation informed LASSO : ",  "\n") )
        names(lasso.agree.naive) = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0", "ridge" )
        names(lasso.agree.naive) = c("1se", "min", "1seR", "minR", "1seR.G0", "minR.G0", "ridge" )
        print( round( lasso.agree.naive , digits=4) )
      }
    }        
    
    ## XGB #####################################################################
    
    if (doxgb == 1) { 
      if (do_ncv == 1) {
      if (short==1) {
        if (family %in% c("cox","binomial")) { 
          cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed XGBoost : ",  "\n") )
        } else { 
          cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed XGBoost : ",  "\n") )
        } 
        print( round( xgbAveAgree , digits=4) )
        
      } else {
        cat(paste0("\n\n" , " Nested Cross Validation averages for XGBoost model : ", "\n") )
        cat(paste0("\n" , "      ", perunit, ": ", "\n") )
        print( round( xgbAveDevian , digits = 4) )
      
        cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
        print( round( xgbAveLincal , digits = 4) ) 
      
        if (family %in% c("cox","binomial")) { 
          cat(paste0("\n", "      average agreement (concordance) :         ", "\n") )
        } else { 
          cat(paste0("\n", "      agreement (r-square) :            ", "\n") )
        } 
      
        print( round( xgbAveAgree , digits = 4) ) 
      
#        cat(paste0("\n", " Cross validation informed XGBoost model : ",  "\n") )
      }
      }
      if ((short != 1) |  (do_ncv == 0)) {
        cat(paste0("\n", " Naive agreement for cross validation informed XGBoost model : ",  "\n") )
        print( round( xgb.agree.naive , digits=4) )
      }
    }                                                                             ## how to get rid of [1,]  ??
    
    ##### ANN ##################################################################
    
    if (doann == 1) { 
      if (do_ncv == 1) {
        if (short==1) {
          if (family %in% c("cox","binomial")) { 
            cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed Neural Network : ",  "\n") )
          } else { 
            cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed Neural Network : ",  "\n") )
          } 
          print( round( annAveAgree , digits = 4) ) 
        } else {
          cat(paste0("\n\n" , " Nested Cross Validation averages for Neural Network : ", "\n") )
        
          cat(paste0("\n" , "      ", perunit, ": ", "\n") )
          print( round( annAveDevian , digits = 4) )
        
          cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
          print( round( annAveLincal , digits = 4) ) 
        
          if (family %in% c("cox","binomial")) { 
            cat(paste0("\n", "      average agreement (concordance) :         ", "\n") )
          } else { 
            cat(paste0("\n", "      agreement (r-square) :            ", "\n") )
          } 
        
          print( round( annAveAgree , digits = 4) ) 
        
          cat(paste0("\n", " Cross validation informed Neural Network : ",  "\n") )
        }
      }
      if ((short != 1) | (do_ncv == 0)) {
        cat(paste0("\n", "      naive agreement : ",  "\n") )
        print( round( ann.agree.naive, digits=4) )
      }
    }                                                                           
    
    ##### RPART ################################################################
    
    if (dorpart == 1) { 
      if ((short==1) & (do_ncv == 1)) {
        if (family %in% c("cox","binomial")) { 
          cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed RPART : ",  "\n") )
        } else { 
          cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed RPART : ",  "\n") )
        } 
        print( round( rpartAveAgree , digits = 4) ) 
      }
      if ((short == 0) & (do_ncv == 1)) {
        cat(paste0("\n\n" , " Nested Cross Validation averages for RPART model : ", "\n") )
      
        cat(paste0("\n" , "      ", perunit, ": ", "\n") )
        print( round( rpartAveDevian , digits = 4) )
      
        cat(paste0("\n", "      average number of terms used in cv informed models : ", "\n") )
        print(round( rpartAveNZero, digits=1) ) 
      
        cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
        print( round( rpartAveLincal , digits = 4) ) 
      
        if (family %in% c("cox","binomial")) { 
          cat(paste0("\n", "      average agreement (concordance) :         ", "\n") )
        } else { 
          cat(paste0("\n", "      agreement (r-square) :            ", "\n") )
        } 
      }      
      if ((short == 0) | (do_ncv==0)) {
        print( round( rpartAveAgree , digits = 4) ) 
      
        cat(paste0("\n", " Cross validation informed RPART model : ",  "\n") )
        cat(paste0("\n", "      number of terms used in model : ", "\n") )
        print(round( rpart.nzero, digits=1)[c(3,2,1,6,5,4)] ) 

        cat(paste0("\n", "      naive agreement : ",  "\n") )
        print( round( rpart.agree.naive , digits=4)[c(3,2,1,6,5,4)] )
      }
    }      
    
    ##### STEP #################################################################
      
    if ( (dostep==1) & (short==1) & (do_ncv == 1)) {
      if (family %in% c("cox","binomial")) { 
        cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed STEPWISE : ",  "\n") )
      } else { 
        cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed STEPWISE : ",  "\n") )
      } 
      print( round(StepAveAgree[c(1,2)], digits = 4) ) 
    }
    
    if ( (doaic==1) & (short==1) & (do_ncv == 1) ) {
      if (family %in% c("cox","binomial")) { 
        cat(paste0("\n", " Cross validation agreement (concordance) for the STEPWISE based AIC : ",  "\n") )
      } else { 
        cat(paste0("\n", " Cross validation agreement (r-square) for the STEPWISE based AIC : ",  "\n") )
      } 
       print( round(StepAveAgree[c(3)], digits = 4) ) 
    }
    
    if ((dostep==1) & (short!=1) & (do_ncv == 1)) {
#      cv.stepreg.fit = object$cv.stepreg.fit
#      cv.stepreg.fit.df = cv.stepreg.fit$best.df    #func.fit.df
      cat(paste0("\n\n", " Nested Cross Validation STEPWISE regression model (df): ", "\n") )   
      cat(paste0("      Average linear calibration coefficient: ", round(StepAveLincal[1] ,  digits = 3), "\n") )    
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
      
      cat(paste0("\n", " Nested Cross Validation STEPWISE regression model (p): ", "\n") )   
      cat(paste0("      Average linear calibration coefficient p : ", round(StepAveLincal[1] ,  digits = 3), "\n") )    
      cat(paste0("      Average deviance : ", round(StepAveDevian[2] ,  digits = 4), "\n") )    
      cat(paste0("      Average model p  : ", round( StepAve_p [1], digits=3 ), "\n") )
      cat(paste0("      Average model df : ", round( StepAve_df[2], digits=2 ), "\n") )
      if (family %in% c("cox","binomial")) {
        cat(paste0("      Concordance      : ", round(StepAveAgree[2],  digits = 4), "\n") )        
        cat(paste0(" Naive concordance based upon the same (all) data as model derivation (p): ",
                   round(cv.stepreg.fit$cvfit.p[6], digits=4) , "\n") )
#        cat(paste0(" Naive concordance based upon the same (all) data as model derivation (df): ",
#                   round(cv.stepreg.fit$cvfit.df[6], digits=4) , "\n") )
      } else {
        cat(paste0("      R-square         : ", round(StepAveAgree[2],  digits = 4), "\n") )  
        cat(paste0(" Naive R-square based upon the same (all) data as model derivation (p): ",
                   round(cv.stepreg.fit$cvfit.p[6], digits=4) , "\n") )
#        cat(paste0(" Naive R-square based upon the same (all) data as model derivation (df): ",
#                   round(cv.stepreg.fit$cvfit.df[6], digits=4) , "\n") )
      }
    }  
    
#    if ((dostep == 1) & (short == 1) & (do_ncv == 0)) {
#    if ((dostep == 1) & (short==1) & (do_ncv==0)) {
      if ((dostep==1) & ((short!=1) | (do_ncv == 0)))  {
      if (family %in% c("cox","binomial")) {
        cat(paste0("\n Naive concordance based upon the same (all) data as model derivation (df): ",
                   round(cv.stepreg.fit$cvfit.df[6], digits=4) , "\n") )
        cat(paste0("\n Naive concordance based upon the same (all) data as model derivation (p): ",
                   round(cv.stepreg.fit$cvfit.p[6], digits=4) , "\n") )
      } else {
        cat(paste0("\n Naive R-square based upon the same (all) data as model derivation (df): ",
                   round(cv.stepreg.fit$cvfit.df[6], digits=4) , "\n") )
        cat(paste0("\n Naive R-square based upon the same (all) data as model derivation (p): ",
                   round(cv.stepreg.fit$cvfit.p[6], digits=4) , "\n") )
      }
    }
    
    if ((doaic==1) & (short!=1) & (do_ncv==1)) {
      cat(paste0("\n Cross Validation results for STEPWISE regression model: (AIC)", "\n") ) 
      cat(paste0("      Average linear calibration coefficient : ", round(StepAveLincal[3] ,  digits = 3), "\n") )    
      cat(paste0("      Average deviance : ", round(StepAveDevian[3],  digits = 4), "\n") )    
      cat(paste0("      Average model df : ", round( StepAve_df[3], digits=2 ), "\n") )
      cat(paste0("      Concordance      : ", round( StepAveAgree[3], digits = 4), "\n") )   
    }
    
    if ((doaic==1) & ((short!=1) | (do_ncv == 0)))  {
      if (do_ncv == 0) { cat("\n") }
      if (family %in% c("cox","binomial")) {
        cat(paste0(" Naive concordance based upon the same (all) data as model derivation (AIC) : ", 
                     round(object$func.fit.aic$aic.fit.n[6], digits=4) , "\n") )
      } else {
        cat(paste0(" Naive R-square based upon the same (all) data as model derivation (AIC) : ",         
                   round( object$func.fit.aic$aic.fit.n[6], digits=4 ), "\n") ) 
        #        round( 1-var(all.fit.aic$residuals)/var(all.fit.aic$y), digits=4), "\n") ) 
      }
    }
    
    ##### END STEP #################################################################
    cat("\n")  
  } #### end of if summmary 
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
  if (class(object) %in% c("nested.glmnetr","nested.stepreg")) { object = object$cv_glmnet_fit }
  retobj = predict.cv.glmnetr( object=object, xs_new=xs_new, lam=lam, gam=gam, comment=comment) 
  return(retobj)
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
  tuning  = object$tuning
  dolasso = tuning[4]
  doann   = tuning[5]
  doxgb   = tuning[6]
  dorpart = tuning[7]
  dostep  = tuning[8]
  doaic   = tuning[9]
  lasso.agree.cv  = object$lasso.agree.cv
  xgb.agree.cv    = object$xgb.agree.cv
  rpart.agree.cv  = object$rpart.agree.cv
  step.agree.cv   = object$step.agree.cv
  if (dolasso == 1) {
    cat ("\n", "C lasso.min   - lasso.minR , ") ;  glmnetr.compcv0(lasso.agree.cv[,2] , lasso.agree.cv[,4]) 
    cat ("\n", "C lasso.min   - lasso.minR0, ") ;  glmnetr.compcv0(lasso.agree.cv[,2] , lasso.agree.cv[,6]) ;  
    cat ("\n", "C lasso.minR  - lasso.minR0, ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , lasso.agree.cv[,6]) ;  cat("\n")
  }

  if (doxgb == 1) {
    cat ("\n", "C XGBoost (tuned) - XGBoost (simple), ") ;  glmnetr.compcv0(xgb.agree.cv[,4] , xgb.agree.cv[,1]) ;  cat("\n")
  }
  
  if ((dolasso == 1) & (doxgb == 1)) {
    cat ("\n", "C lasso.min   - XGBoost (tuned), ") ;  glmnetr.compcv0(lasso.agree.cv[,2] , xgb.agree.cv[,4]) 
    cat ("\n", "C lasso.minR  - XGBoost (tuned), ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , xgb.agree.cv[,4]) 
    cat ("\n", "C lasso.minR0 - XGBoost (tuned), ") ;  glmnetr.compcv0(lasso.agree.cv[,6] , xgb.agree.cv[,4])  ;  cat("\n")
    cat ("\n", "C lasso.min   - XGBoost (simple),") ;  glmnetr.compcv0(lasso.agree.cv[,2] , xgb.agree.cv[,1]) 
    cat ("\n", "C lasso.minR  - XGBoost (simple),") ;  glmnetr.compcv0(lasso.agree.cv[,4] , xgb.agree.cv[,1]) 
    cat ("\n", "C lasso.minR0 - XGBoost (simple),") ;  glmnetr.compcv0(lasso.agree.cv[,6] , xgb.agree.cv[,1])  ;  cat("\n")
    
  }
  
  if ((dolasso == 1) & (dorpart == 1)) {
    cat ("\n", "C lasso.min   - RPART (df), ") ;  glmnetr.compcv0(lasso.agree.cv[,2] , rpart.agree.cv[,1]) 
    cat ("\n", "C lasso.minR  - RPART (df), ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , rpart.agree.cv[,1]) 
    cat ("\n", "C lasso.minR0 - RPART (df), ") ;  glmnetr.compcv0(lasso.agree.cv[,6] , rpart.agree.cv[,1])  ;  cat("\n")
  }
  
  if ((dolasso == 1) & (dostep == 1)) {
    cat ("\n", "C lasso.min   - step (df),   ") ;  glmnetr.compcv0(lasso.agree.cv[,2] , step.agree.cv[,1]) 
    cat ("\n", "C lasso.minR  - step (df),   ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , step.agree.cv[,1]) 
    cat ("\n", "C lasso.minR0 - step (df),   ") ;  glmnetr.compcv0(lasso.agree.cv[,6] , step.agree.cv[,1])     ;  cat("\n")
    
    cat ("\n", "C lasso.min   - step (p),   ") ;  glmnetr.compcv0(lasso.agree.cv[,2] , step.agree.cv[,2]) 
    cat ("\n", "C lasso.minR  - step (p),   ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , step.agree.cv[,2]) 
    cat ("\n", "C lasso.minR0 - step (p),   ") ;  glmnetr.compcv0(lasso.agree.cv[,6] , step.agree.cv[,2])    ;  cat("\n")
    
    cat ("\n", "C step (df) - step (p),     ") ;  glmnetr.compcv0(step.agree.cv[,1]      , step.agree.cv[,2])    ;  cat("\n")
  }
  
  if ((dolasso == 1) & (doaic == 1)) {
    cat ("\n", "C lasso.min   - step (AIC),  ") ;  glmnetr.compcv0(lasso.agree.cv[,2] , step.agree.cv[,3]) 
    cat ("\n", "C relax.minR  - step (AIC),  ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , step.agree.cv[,3]) 
    cat ("\n", "C lasso.minR0 - step (AIC),  ") ;  glmnetr.compcv0(lasso.agree.cv[,6] , step.agree.cv[,3]) ;  cat("\n")    
  }
  
  if ((doxgb == 1) & (dorpart == 1)) {
    cat ("\n", "C XGBoost (tuned) - RPART (cp=0), ") ;  glmnetr.compcv0(xgb.agree.cv[,4] , rpart.agree.cv[,1]) ;  cat("\n")
  }
  
  if ((doxgb == 1) & (dostep == 1)) {
    cat ("\n", "C XGBoost (tuned) - step (df), ") ;  glmnetr.compcv0(xgb.agree.cv[,4] ,  step.agree.cv[,1]) ;  cat("\n")
    cat ("\n", "C XGBoost (tuned) - step ( p), ") ;  glmnetr.compcv0(xgb.agree.cv[,4] ,  step.agree.cv[,2]) ;  cat("\n")
  }
  
  if ((doxgb == 1) & (doaic == 1)) {
    cat ("\n", "C XGBoost (tuned) - step (AIC), ") ;  glmnetr.compcv0(xgb.agree.cv[,4] ,  step.agree.cv[,3]) ;  cat("\n")
  }

  cat("\n")
}

####################################################################################################################################
####################################################################################################################################

