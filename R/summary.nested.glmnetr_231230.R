####################################################################################################################################
####################################################################################################################################
#' Print an abbreviated summary of a nested.glmnetr() output object
#'
#' @param x a nested.glmnetr() output object.  
#' @param ... additional pass through inputs for the print function.
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
#' @param digits digits for printing of deviances, linear calibration coefficients 
#' and agreement (concordances and R-squares).
#' @param Call 1 to print call used in generation of the object, 0 or NULL to not print 
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
# object = nestedcv ; cvfit = FALSE ; printg1 = FALSE ; short = 0 ; digits = 3 ; 
#'
summary.nested.glmnetr = function(object, cvfit=FALSE, printg1=FALSE, short=0, digits=3, Call=NULL, ...) {
  
  if (cvfit==TRUE) {
    cv_glmnet_fit = object$cv_glmnet_fit  
    summary(cv_glmnet_fit,printg1=printg1)
  } else {
    if (!is.null(Call)) { 
      if (Call != 0) { Call = object$Call 
      } else { Call = NULL }
    }
    sample  = object$sample 
    tuning  = object$tuning
    fits    = object$fits 
    dolasso = fits[1]
    doxgb   = fits[2]
    dorf    = fits[3]
    dorpart = fits[4]
    doann   = fits[5]
    dostep  = fits[6]
    doaic   = fits[7]
    ensemble = object$ensemble
    do_ncv  = object$do_ncv

    family = sample[1]
    sample[6] = round(as.numeric(sample[6]), digits=digits)
    if (family=="cox") { names(sample)[6]="null.dev/events" 
    } else { names(sample)[6] = "null.dev/obs"  }

    if (dolasso == 1) { 
      lasso.nzero.cv  = object$lasso.nzero.cv
      lasso.devian.cv = object$lasso.devian.cv
      lasso.cal.devian.cv = object$lasso.cal.devian.cv
      lasso.lincal.cv = object$lasso.lincal.cv
      lasso.agree.cv  = object$lasso.agree.cv
      lasso.devian.naive = object$lasso.devian.naive 
      lasso.agree.naive = object$lasso.agree.naive
      lassoAveNZero  = colMeans(lasso.nzero.cv)  
      lassoAveDevian = colMeans(lasso.devian.cv)  
      lassoCalAveDevian = colMeans(lasso.cal.devian.cv)  
      lassoAveLincal = colMeans(lasso.lincal.cv)  
      lassoAveAgree  = colMeans(lasso.agree.cv)
      if (family == "gaussian") { 
        lasso.agree.naive = lasso.agree.naive ^2
        lassoAveAgree = lassoAveAgree ^2
      }
      ## sqrt( apply(lasso.agree.cv,2,var) ) 
    } 
    if (doxgb==1) { 
      xgb.devian.cv   = object$xgb.devian.cv
      xgb.lincal.cv   = object$xgb.lincal.cv
      xgb.agree.cv    = object$xgb.agree.cv
      xgb.agree.naive = object$xgb.agree.naive 
      xgbAveDevian = colMeans(xgb.devian.cv)   
      xgbAveLincal = colMeans(xgb.lincal.cv)   
      xgbAveAgree  = colMeans(xgb.agree.cv)    
      if (family == "gaussian") { 
        xgb.agree.naive = xgb.agree.naive ^2 
        xgbAveAgree = xgbAveAgree ^2 
      }
    } 
    if (dorf == 1) { 
      rf.devian.cv   = object$rf.devian.cv
      rf.lincal.cv   = object$rf.lincal.cv
      rf.agree.cv    = object$rf.agree.cv
      rf.mtry.cv     = object$rf.mtry.cv
      rf.mtry        = object$rf.mtry
      rf.agree.naive = object$rf.agree.naive 
      rfAveDevian = colMeans(rf.devian.cv)   
      rfAveLincal = colMeans(rf.lincal.cv)   
      rfAveAgree  = colMeans(rf.agree.cv)    
      rfAveMtry   = colMeans(rf.mtry.cv)   
      if (family == "gaussian") { 
        rf.agree.naive = rf.agree.naive ^2 
        rfAveAgree  = rfAveAgree ^2 
      }
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
      rpartAveNZero  = rpartAveNZero [ c(1:9) ]
      rpartAveDevian = rpartAveDevian[ c(1:9) ] 
      rpartAveLincal = rpartAveLincal[ c(1:9) ]
      rpartAveAgree = rpartAveAgree  [ c(1:9) ]
      rpartAveNZero = rpartAveNZero  [ c(1:9) ]
      if (family == "gaussian") { 
        rpart.agree.naive = rpart.agree.naive ^2 
        rpartAveAgree = rpartAveAgree ^2 
      }
    }
    if (doann ==1) { 
#      nms = c("Uninformed", "lasso terms", "lasso feat.", "init w's", "update w's", "offset") 
      nms = c("Uninformed", "lasso feat", "lasso w's", "lasso update", "lasso terms", "l/lasso feat", "l/lasso w's", "l/lasso update") 
      ann.devian.cv = object$ann.devian.cv
      ann.lincal.cv = object$ann.lincal.cv
      ann.agree.cv  = object$ann.agree.cv 
      ann.agree.naive = object$ann.agree.naive[(ensemble==1)]
      annAveDevian = colMeans(ann.devian.cv)[(ensemble==1)]  
      annAveLincal = colMeans(ann.lincal.cv)[(ensemble==1)]
      annAveAgree  = colMeans(ann.agree.cv )[(ensemble==1)]
      if (family == "gaussian") { 
        ann.agree.naive = ann.agree.naive ^2
        annAveAgree = annAveAgree ^2 
      }
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
      if (family == "gaussian") { 
        StepAveAgree = StepAveAgree ^2 
      }
    }
    if (doann == 1) { 
      if        (ensemble[4]==1) { ann_cv = object$ann_fit_4 ; whichann = nms[4] ;
      } else if (ensemble[8]==1) { ann_cv = object$ann_fit_8 ; whichann = nms[8] ;
      } else if (ensemble[3]==1) { ann_cv = object$ann_fit_3 ; whichann = nms[3] ;
      } else if (ensemble[7]==1) { ann_cv = object$ann_fit_7 ; whichann = nms[7] ;
      } else if (ensemble[2]==1) { ann_cv = object$ann_fit_2 ; whichann = nms[2] ;
      } else if (ensemble[6]==1) { ann_cv = object$ann_fit_6 ; whichann = nms[6] ;
      } else if (ensemble[1]==1) { ann_cv = object$ann_fit_1 ; whichann = nms[1] ;
      } else if (ensemble[5]==1) { ann_cv = object$ann_fit_5 ; whichann = nms[5] ;
      }
    }
    if (dostep==1) { cv.stepreg.fit    = object$cv.stepreg.fit }  
    if (doaic ==1) { func.fit.aic      = object$func.fit.aic   }
  
    if (!is.null(Call)) { 
      cat(paste0("\n     function call :\n\n"))  
          print(Call) 
#          cat(paste0("\n")) 
    } 
    
    cat(paste0("\n"  , " Sample information including number of records, "))     
    if (family %in% c("cox","binomial")) { cat(paste0("events, ")) }
    cat(paste0( "number of columns in", "\n", " design (predictor, X) matrix, and df (rank) of design matrix: ", "\n") )
    if (family %in% c("cox","binomial")) { print(sample) 
    } else { print(sample[-3]) }
    
    if (!is.null(object$dep_names)) { 
      cat(paste0("\n"  , " Dependent Variable(s) : ", "\n") )
      print(object$dep_names) 
    }

    if ((dolasso ==0) & (dostep==0) & (doaic==0)) { 
      cat(paste0("\n"  , " Tuning parameters for models : ", "\n") )
      print(object$tuning[c(1,2)]) 
    } else if ((dostep==0) & (doaic==0)) { 
      cat(paste0("\n"  , " Tuning parameters for models : ", "\n") )
      print(object$tuning[c(1:5)]) 
    } else { 
      cat(paste0("\n"  , " Tuning parameters for models : ", "\n") )
      print(object$tuning) 
    }
    
    if (doann  == 1) {
      cat(paste0("\n"  , " Tuning parameters for  ", whichann, "  ANN model : ", "\n") )
      print(ann_cv$modelsum[c(1:9,11:12)]) 
    }
    
#    if (doaic==1) {
#      cat(paste0("\n", " Average deviance for null model", "\n") )    ## pull other data applicable to all models 
#      print( round(StepAveDevian[1], digits = digits ) ) 
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
          print( round( lassoAveAgree , digits = digits) ) 
        } else {

#          cat(paste0("\n" , " Original Sample Null Deviance : ", "\n") )
#          cat(paste0( "       ", round( object$null.devian0.naive , digits = digits) ,"\n") ) 
#          object$cv_glmnet_fit
#          object$lasso.devain.naive
          
#          object$null.devian0.cv
          cat(paste0("\n" , " Average Null Deviance for leave out folds in outer loop : ", "\n") )
          cat(paste0( "       ", round( mean(object$null.devian0.cv) , digits = digits),"\n") )

          cat(paste0("\n\n" , " Nested Cross Validation averages for LASSO (1se and min), Relaxed LASSO, and gamma=0 LASSO : ", "\n") )
          cat(paste0("\n" , "      ", perunit, ": ", "\n") )
          print( round( lassoAveDevian , digits = digits) )
          
#          cat(paste0("\n\n" , " Nested Cross Validation averages for Linear Calibrated LASSO : ", "\n") )
          cat(paste0("\n" , "      ", perunit, "(linerly calibrated) : ", "\n") )
          print( round( lassoCalAveDevian , digits = digits) )
        
          cat(paste0("\n" , "      number of nonzero model terms : ", "\n") )
          print( round( lassoAveNZero , digits = max((digits-2),1) ) )
      
          cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
          print( round( lassoAveLincal , digits = digits) ) 
      
          # In summary.nested.glmnetr state only "concordance" or "R-square" in output
      
          if (family %in% c("cox","binomial")) { 
            cat(paste0("\n", "      agreement (concordance) :         ", "\n") )
          } else { 
            cat(paste0("\n", "      agreement (r-square) :            ", "\n") )
          } 
          print( round( lassoAveAgree , digits = digits) ) 
        }
      }
      if ((short != 1) | (do_ncv == 0)) {
        cat(paste0("\n", " Naive deviance for cross validation informed LASSO : ",  "\n") )
        print(  lasso.devian.naive ) 
        
        cat(paste0("\n", " Number of non-zero terms in cross validation informed LASSO : \n") )
        lassoNZero = c(rep(0,7))
        lassoNZero[1] = object$cv_glmnet_fit$nzero [ object$cv_glmnet_fit$index[2] ]
        lassoNZero[2] = object$cv_glmnet_fit$nzero [ object$cv_glmnet_fit$index[1] ]
        lassoNZero[3] = object$cv_glmnet_fit$relaxed$nzero.1se
        lassoNZero[4] = object$cv_glmnet_fit$relaxed$nzero.min
        lassoNZero[5] = object$cv_glmnet_fit$nzero [ object$cv_glmnet_fit$relaxed$index.g0[2] ] 
        lassoNZero[6] = object$cv_glmnet_fit$nzero [ object$cv_glmnet_fit$relaxed$index.g0[1] ] 
        lassoNZero[7] = object$cv_ridge_fit$nzero[object$cv_ridge_fit$index][1]
        names(lassoNZero) = c("1se", "min", "1seR", "minR", "1seR.G0", "minR.G0", "ridge" )
        print( round( lassoNZero , digits=digits) )
        
        cat(paste0("\n", " Naive agreement for cross validation informed LASSO : ",  "\n") )
        names(lasso.agree.naive) = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0", "ridge" )
        names(lasso.agree.naive) = c("1se", "min", "1seR", "minR", "1seR.G0", "minR.G0", "ridge" )
        print( round( lasso.agree.naive , digits=digits) )
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
        print( round( xgbAveAgree , digits = digits) )
        
      } else {
        cat(paste0("\n\n" , " Nested Cross Validation averages for XGBoost model : ", "\n") )
        cat(paste0("\n" , "      ", perunit, ": ", "\n") )
        print( round( xgbAveDevian , digits = digits) )
      
        cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
        print( round( xgbAveLincal , digits = digits) ) 
      
        if (family %in% c("cox","binomial")) { 
          cat(paste0("\n", "      average agreement (concordance) :         ", "\n") )
        } else { 
          cat(paste0("\n", "      agreement (r-square) :            ", "\n") )
        } 
      
        print( round( xgbAveAgree , digits = digits) ) 
      
#        cat(paste0("\n", " Cross validation informed XGBoost model : ",  "\n") )
      }
      }
      if ((short != 1) |  (do_ncv == 0)) {
        cat(paste0("\n", " Naive agreement for cross validation informed XGBoost model : ",  "\n") )
        print( round( xgb.agree.naive , digits=digits) )
      }
    }                                                                             ## how to get rid of [1,]  ??
    
    ##### Random Forest ########################################################
    
    if (dorf == 1) { 
      if ((short==1) & (do_ncv == 1)) { 
        if (family %in% c("cox","binomial")) { 
          if (sum(ensemble[2:8])==0) { 
            cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed Random Forest :", 
                       round( rfAveAgree[1], digits = digits))) 
          } else {
            cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed Random Forest : \n") )
            print( round( rfAveAgree[1:2] , digits = digits)) 
          }
        } else { 
          if (sum(ensemble[2:8])==0) { 
            cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed Random Forest : ", 
                       round( rfAveAgree[1] , digits = digits), "\n") )
          } else {
            cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed Random Forest : \n") )
            print( round( rfAveAgree[1:3] , digits = digits)) 
          }
        } 
      }
      if ((short == 0) & (do_ncv == 1)) {
        if (sum(ensemble[2:8])==0) {
          cat(paste0("\n\n" , " Nested Cross Validation averages for Random Forest : ", "\n") )
          cat(paste0("\n" , "      ", perunit, ": ", round( rfAveDevian[1] , digits = digits), "\n") )
          cat(paste0("\n", "      average number of variables random selected for the RF : ", round( rfAveMtry[1], digits=1) , "\n") )
          cat(paste0("\n", "      linear calibration coefficient :   ", round( rfAveLincal[1] , digits = digits) , "\n") )
          if (family %in% c("cox","binomial")) { 
            cat(paste0("\n", "      average agreement (concordance) :  ", round( rfAveAgree[1] , digits = digits), "\n") )
          } else { 
            cat(paste0("\n", "      agreement (r-square) :             ", round( rfAveAgree[1] , digits = digits), "\n") )
          } 
        } else {
          cat(paste0("\n\n" , " Nested Cross Validation averages for Random Forest : ", "\n") )
          cat(paste0("\n" , "      ", perunit, ": \n")) 
          if (family %in% c("cox","binomial")) { 
            print( round( rfAveDevian[1:2] , digits = digits) )
            cat(paste0("\n", "      average number of variables randomly selected for the RF : \n" ) )
            print( round( rfAveMtry  [1:2], digits=1) ) 
            cat(paste0("\n", "      linear calibration coefficient :   \n") )
            print( round( rfAveLincal[1:2] , digits = digits) )
            cat(paste0("\n", "      average agreement (concordance) :  \n") )
            print( round( rfAveAgree [1:2] , digits = digits) )
          } else { 
            print( round( rfAveDevian[1:3] , digits = digits) )
            cat(paste0("\n", "      average number of variables randomly selected for the RF : \n" ) )
            print( round( rfAveMtry  [1:3], digits=1) ) 
            cat(paste0("\n", "      linear calibration coefficient :   \n") )
            print( round( rfAveLincal[1:3] , digits = digits) )
            cat(paste0("\n", "      agreement (r-square) :             \n") )
            print( round( rfAveAgree [1:3] , digits = digits) )
          } 
        }
      } 
      if (((short==1) & (do_ncv == 0)) | ((short==0) & (do_ncv == 1))) {
        if (family %in% c("cox","binomial")) { 
          cat(paste0("\n", " Naive Random Forest agreement (concordance) :  \n") )
          print( round( rf.agree.naive[1:2] , digits = digits) )
        } else if (family %in% c("gaussian")) { 
          cat(paste0("\n", " Naive Random Forest agreement (r-square) :      \n") )
          print( round( rf.agree.naive[1:3] , digits = digits) )
        } 
      }
    }
#      if ((short == 0) | (do_ncv==0)) {
#        cat(paste0("\n", " Cross validation informed Random Forest : ",  "\n") )
#        cat(paste0("\n", "      number of terms used in model : ", "\n") )
#        print(round( rf.mtry[3], digits=1)) 
#        cat(paste0("\n", "      naive agreement : ",  "\n") )
#        print( round( rf.agree.naive[3] , digits=digits)[c(3,2,1,6,5,4,9,8,7)] )
#      }
#    }      
    
    ##### RPART ################################################################
    
    if (dorpart == 1) { 
      if ((short==1) & (do_ncv == 1)) {
        if (family %in% c("cox","binomial")) { 
          cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed RPART : ",  "\n") )
        } else { 
          cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed RPART : ",  "\n") )
        } 
        print( round( rpartAveAgree , digits = digits) ) 
      }
      if ((short == 0) & (do_ncv == 1)) {
        cat(paste0("\n\n" , " Nested Cross Validation averages for RPART : ", "\n") )
        
        cat(paste0("\n" , "      ", perunit, ": ", "\n") )
        print( round( rpartAveDevian , digits = digits) )
        
        cat(paste0("\n", "      average number of terms used in cv informed models : ", "\n") )
        print(round( rpartAveNZero, digits=1) ) 
        
        cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
        print( round( rpartAveLincal , digits = digits) ) 
        
        if (family %in% c("cox","binomial")) { 
          cat(paste0("\n", "      average agreement (concordance) :         ", "\n") )
        } else { 
          cat(paste0("\n", "      agreement (r-square) :            ", "\n") )
        } 
        print( round( rpartAveAgree , digits = digits) ) 
      }      
      if ((short == 0) | (do_ncv==0)) {
        cat(paste0("\n", " Cross validation informed RPART : ",  "\n") )
        cat(paste0("\n", "      number of terms used in model : ", "\n") )
        print(round( rpart.nzero, digits=1)[c(3,2,1,6,5,4,9,8,7)] ) 
        
        cat(paste0("\n", "      naive agreement : ",  "\n") )
        print( round( rpart.agree.naive , digits=digits)[c(3,2,1,6,5,4,9,8,7)] )
      }
    }      
    
    ##### ANN ##################################################################
    
    if (doann == 1) { 
      if (do_ncv == 1) {
        if (short==1) {
          if (family %in% c("cox","binomial")) { 
            cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed Neural Network : ",  "\n") )
          } else { 
            cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed Neural Network : ",  "\n") )
          } 
          print( round( annAveAgree , digits = digits) ) 
        } else {
          cat(paste0("\n\n" , " Nested Cross Validation averages for neural network : ", "\n") )
        
          cat(paste0("\n" , "      ", perunit, ": ", "\n") )
          print( round( annAveDevian , digits = digits) )
        
          cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
          print( round( annAveLincal , digits = digits) ) 
        
          if (family %in% c("cox","binomial")) { 
            cat(paste0("\n", "      average agreement (concordance) :         ", "\n") )
          } else { 
            cat(paste0("\n", "      agreement (r-square) :            ", "\n") )
          } 
        
          print( round( annAveAgree , digits = digits) ) 
        
        }
      }
      if ((short != 1) | (do_ncv == 0)) {
        cat(paste0("\n", " Cross validation informed neural network : ",  "\n") )
        cat(paste0("\n", "      naive agreement : ",  "\n") )
        print( round( ann.agree.naive, digits=digits) )
        cat(paste0("\n") )
      }
    }                                                                           
    
    ##### STEP #################################################################
      
    if ( (dostep==1) & (short==1) & (do_ncv == 1)) {
      if (family %in% c("cox","binomial")) { 
        cat(paste0("\n", " Nested cross validation agreement (concordance) for cross validation informed STEPWISE : ",  "\n") )
      } else { 
        cat(paste0("\n", " Nested cross validation agreement (r-square) for cross validation informed STEPWISE : ",  "\n") )
      } 
      print( round(StepAveAgree[c(1,2)], digits = digits) ) 
    }
    
    if ( (doaic==1) & (short==1) & (do_ncv == 1) ) {
      if (family %in% c("cox","binomial")) { 
        cat(paste0("\n", " Cross validation agreement (concordance) for the STEPWISE based AIC : ",  "\n") )
      } else { 
        cat(paste0("\n", " Cross validation agreement (r-square) for the STEPWISE based AIC : ",  "\n") )
      } 
       print( round(StepAveAgree[c(3)], digits = digits) ) 
    }
    
    if ((dostep==1) & (short!=1) & (do_ncv == 1)) {
#      cv.stepreg.fit.df = cv.stepreg.fit$best.df    #func.fit.df ; cv.stepreg.fit$best.p 
      cat(paste0("\n\n", " Nested Cross Validation STEPWISE regression model (df): ", "\n") )   
      cat(paste0("      Average linear calibration coefficient: ", round(StepAveLincal[1] ,  digits = 3), "\n") )    
      cat(paste0("      Average deviance : ", round(StepAveDevian[1]   ,  digits = digits), "\n") )    
      cat(paste0("      Average model df : ", round( StepAve_df[1], digits=2 ), "\n") )
      if (family %in% c("cox","binomial")) {
        cat(paste0("      Concordance      : ", round(StepAveAgree[1],  digits = digits), "\n") )        
        cat(paste0(" Naive concordance based upon the same (all) data as model derivation (df): ",
                   round(cv.stepreg.fit$cvfit.df[6], digits=digits) , "\n") )
      } else {
        cat(paste0("      R-square         : ", round(StepAveAgree[1],  digits = digits), "\n") )        
        cat(paste0(" Naive R-square based upon the same (all) data as model derivation (df): ",
                   round(cv.stepreg.fit$cvfit.df[6], digits=digits) , "\n") )
      } 
      cat(paste0("    Model df ", length(cv.stepreg.fit$func.fit.df$coefficients), "\n")) 
      
      cat(paste0("\n", " Nested Cross Validation STEPWISE regression model (p): ", "\n") )   
      cat(paste0("      Average linear calibration coefficient p : ", round(StepAveLincal[1] ,  digits = 3), "\n") )    
      cat(paste0("      Average deviance : ", round(StepAveDevian[2] ,  digits = digits), "\n") )    
      cat(paste0("      Average model p  : ", round( StepAve_p [1], digits=3 ), "\n") )
      cat(paste0("      Average model df : ", round( StepAve_df[2], digits=2 ), "\n") )
      if (family %in% c("cox","binomial")) {
        cat(paste0("      Concordance      : ", round(StepAveAgree[2],  digits = digits), "\n") )        
        cat(paste0(" Naive concordance based upon the same (all) data as model derivation (p): ",
                   round(cv.stepreg.fit$cvfit.p[6], digits=digits) , "\n") )
      } else {
        cat(paste0("      R-square         : ", round(StepAveAgree[2]^2,  digits = digits), "\n") )  
        cat(paste0(" Naive R-square based upon the same (all) data as model derivation (p): ",
                   round(cv.stepreg.fit$cvfit.p[6], digits=digits) , "\n") )
      }
      cat(paste0("    Model df ", length(cv.stepreg.fit$func.fit.p$coefficients), "\n")) 
    }  
  
#    if ((dostep == 1) & (short == 1) & (do_ncv == 0)) {
#    if ((dostep == 1) & (short==1) & (do_ncv==0)) {
#    if ((dostep==1) & ((short!=1) | (do_ncv == 0)))  {
    
    if ((dostep==1) | (doaic==1)) {
      if (do_ncv == 0) { 
        cat(paste0("\n Cross validation informed STEPWISE regression models : \n"))
      }
    }
    
    if (dostep==1) { 
      if (do_ncv == 0) { 
        if (family %in% c("cox","binomial")) {
          cat(paste0("\n Naive concordance based upon the same (all) data as model derivation (df): ",
                     round(cv.stepreg.fit$cvfit.df[6], digits=digits) , "\n") )
          cat(paste0("    Model df ", length(cv.stepreg.fit$func.fit.df$coefficients), "\n")) 
          cat(paste0("\n Naive concordance based upon the same (all) data as model derivation (p): ",
                     round(cv.stepreg.fit$cvfit.p[6], digits=digits) , "\n") )
          cat(paste0("    Model df ", length(cv.stepreg.fit$func.fit.p$coefficients), "\n")) 
        } else {
          cat(paste0("\n Naive R-square based upon the same (all) data as model derivation (df): ",
                     round(cv.stepreg.fit$cvfit.df[6]^2, digits=digits) , "\n") )
          cat(paste0("    Model df ", length(cv.stepreg.fit$func.fit.df$coefficients), "\n")) 
          cat(paste0("\n Naive R-square based upon the same (all) data as model derivation (p): ",
                     round(cv.stepreg.fit$cvfit.p[6]^2, digits=digits) , "\n") )
          cat(paste0("    Model df ", length(cv.stepreg.fit$func.fit.p$coefficients), "\n")) 
        }
      }
    }
    
    if (doaic==1) {
      if ((short!=1) & (do_ncv==1)) {
        cat(paste0("\n Cross Validation results for STEPWISE regression model: (AIC)", "\n") ) 
        cat(paste0("      Average linear calibration coefficient : ", round(StepAveLincal[3] ,  digits = 3), "\n") )    
        cat(paste0("      Average deviance : ", round(StepAveDevian[3],  digits = digits), "\n") )    
        cat(paste0("      Average model df : ", round( StepAve_df[3], digits=2 ), "\n") )
        if (family %in% c('cox','binomial')) {
          cat(paste0("      Concordance      : ", round( StepAveAgree[3], digits = digits), "\n") )   
        } else {
          cat(paste0("      R-square         : ", round( StepAveAgree[3]^2, digits = digits), "\n") )   
        }
      }
    
      if ((short!=1) | (do_ncv == 0))  {
        if (do_ncv == 0) { cat("\n") }
        if (family %in% c("cox","binomial")) {
          cat(paste0(" Naive concordance based upon the same (all) data as model derivation (AIC) : ", 
                       round(object$func.fit.aic$aic.fit.n[6], digits=digits) , "\n") )
        } else {
          cat(paste0(" Naive R-square based upon the same (all) data as model derivation (AIC) : ",         
                     round( object$func.fit.aic$aic.fit.n[6], digits=digits ), "\n") ) 
          #        round( 1-var(all.fit.aic$residuals)/var(all.fit.aic$y), digits=digits), "\n") ) 
        }
        cat(paste0("    Model df ", length(object$func.fit.aic$coefficients), "\n")) 
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
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4
#'
#' @return  A t-test
#' 
#' @importFrom stats t.test 
#'
glmnetr.compcv0 = function(a,b,digits=4) { 
  tdiff = t.test(a-b)
  cat ( paste0(  " estimate (95% CI): ", round(tdiff$estimate, digits=digits), " (", round(tdiff$conf.int[1],digits=digits), ", ", 
                 round(tdiff$conf.int[2],digits=digits), ") , p=", round(tdiff$p.value,digits=digits) ) )
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
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4 
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
glmnetr.compcv = function(object, digits=4) {
  tuning  = object$tuning
  fits    = object$fits 
  dolasso = fits[1]
  doxgb   = fits[2]
  dorf    = fits[3]
  dorpart = fits[4]
  doann   = fits[5]
  dostep  = fits[6]
  doaic   = fits[7]
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

