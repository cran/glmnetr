################################################################################  
#' Summarize a nested.glmnetr() output objects version 0.4-2
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
#' @param pow the power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  Pow is ignored for the family of "cox" and "binomial".   
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
#' @noRd
#'
summary.nested.glmnetr_0_4_2 = function(object, cvfit=FALSE, pow=2, printg1=FALSE, short=0, digits=3, Call=NULL, ...) {
  
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
        lasso.agree.naive = lasso.agree.naive ^pow
        lassoAveAgree = lassoAveAgree ^pow
      }
      ## sqrt( apply(lasso.agree.cv,2,var) ) 
    } 
    if (doxgb==1) { 
      en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
      en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
      en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
      enx = c(en1, en2, en3, en1, en2, en3) 
      xgb.devian.cv   = object$xgb.devian.cv
      xgb.lincal.cv   = object$xgb.lincal.cv
      xgb.agree.cv    = object$xgb.agree.cv
      xgb.agree.naive = object$xgb.agree.naive[enx==1] 
      xgbAveDevian = colMeans(xgb.devian.cv)[enx==1] 
      xgbAveLincal = colMeans(xgb.lincal.cv)[enx==1] 
      xgbAveAgree  = colMeans(xgb.agree.cv)[enx==1] 
      if (family == "gaussian") { 
        xgb.agree.naive = xgb.agree.naive ^pow 
        xgbAveAgree = xgbAveAgree ^pow 
      }
    } 
    if (dorf == 1) { 
      en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
      en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
      en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
      if (family != "gaussian") { en3 = 0 }
      enx = c(en1, en2, en3) 
      rf.devian.cv   = object$rf.devian.cv
      rf.lincal.cv   = object$rf.lincal.cv
      rf.agree.cv    = object$rf.agree.cv
      rf.mtry.cv     = object$rf.mtry.cv
      rf.mtry        = object$rf.mtry
      rf.agree.naive = object$rf.agree.naive[enx==1] 
      rfAveDevian = colMeans(rf.devian.cv)[enx==1]   
      rfAveLincal = colMeans(rf.lincal.cv)[enx==1]   
      rfAveAgree  = colMeans(rf.agree.cv)[enx==1]    
      rfAveMtry   = colMeans(rf.mtry.cv)[enx==1]   
      if (family == "gaussian") { 
        rf.agree.naive = rf.agree.naive ^pow 
        rfAveAgree  = rfAveAgree ^pow 
      }
    } 
    if (dorpart==1) { 
      en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
      en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
      en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
      if (family == "binomial") { en3 = 0 }
      enx = c(en1,en1,en1, en2,en2,en2, en3,en3,en3) 
      rpart.nzero.cv  = object$rpart.nzero.cv
      rpart.devian.cv = object$rpart.devian.cv
      rpart.lincal.cv = object$rpart.lincal.cv
      rpart.agree.cv  = object$rpart.agree.cv
      rpart.nzero     = object$rpart.nzero [c(3,2,1,6,5,4,9,8,7)]          [enx==1]
      rpart.agree.naive = object$rpart.agree.naive [c(3,2,1,6,5,4,9,8,7)]  [enx==1]
      rpartAveNZero  = colMeans(rpart.nzero.cv[,c(3,2,1,6,5,4,9,8,7)])    
      rpartAveDevian = colMeans(rpart.devian.cv[,c(3,2,1,6,5,4,9,8,7)])   
      rpartAveLincal = colMeans(rpart.lincal.cv[,c(3,2,1,6,5,4,9,8,7)])   
      rpartAveAgree  = colMeans(rpart.agree.cv[,c(3,2,1,6,5,4,9,8,7)])   
      rpartAveNZero  = rpartAveNZero [ c(1:9) ] [enx==1]
      rpartAveDevian = rpartAveDevian[ c(1:9) ] [enx==1]
      rpartAveLincal = rpartAveLincal[ c(1:9) ] [enx==1]
      rpartAveAgree = rpartAveAgree  [ c(1:9) ] [enx==1]
      rpartAveNZero = rpartAveNZero  [ c(1:9) ] [enx==1]
      if (family == "gaussian") { 
        rpart.agree.naive = rpart.agree.naive ^pow 
        rpartAveAgree = rpartAveAgree ^pow 
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
        ann.agree.naive = ann.agree.naive ^pow
        annAveAgree = annAveAgree ^pow 
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
        StepAveAgree = StepAveAgree ^pow 
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
    
    if ( pow == 2 ) { gagree = "R-square" } else { gagree = "Correlation" }
    if (family %in% c("cox","binomial")) { gagree = "Concordance" ; pow = 1 } 
    
    ## LASSO ###################################################################
    
    if (dolasso == 1) {
      if (do_ncv == 1) {
        if (short==1) {
          cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed LASSO : ",  "\n") )
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
      
          cat(paste0("\n", "      agreement (", gagree, ") :            ", "\n") )
          print( round( lassoAveAgree , digits = digits) ) 
        }
      }
      if ((short != 1) | (do_ncv == 0)) {
        cat(paste0("\n", " Naive deviance for cross validation informed LASSO : ",  "\n") )
        print( round( lasso.devian.naive , digits=digits ) ) 
        
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
        
        cat(paste0("\n", " Naive agreement (", gagree, ") for cross validation informed LASSO : ",  "\n") )
        names(lasso.agree.naive) = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0", "ridge" )
        names(lasso.agree.naive) = c("1se", "min", "1seR", "minR", "1seR.G0", "minR.G0", "ridge" )
        print( round( lasso.agree.naive , digits=digits) )
      }
    }        
    
    ## XGB #####################################################################
    
    if (doxgb == 1) { 
      if (do_ncv == 1) {
      if (short==1) {
        cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed XGBoost : ",  "\n") )
        print( round( xgbAveAgree , digits = digits) )
      } else {
        cat(paste0("\n\n" , " Nested Cross Validation averages for XGBoost model : ", "\n") )
        cat(paste0("\n" , "      ", perunit, ": ", "\n") )
        print( round( xgbAveDevian , digits = digits) )
      
        cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
        print( round( xgbAveLincal , digits = digits) ) 
      
        cat(paste0("\n", "      agreement (", gagree, ") :            ", "\n") )
        print( round( xgbAveAgree , digits = digits) ) 
      
#        cat(paste0("\n", " Cross validation informed XGBoost model : ",  "\n") )
      }
      }
      if ((short != 1) |  (do_ncv == 0)) {
        cat(paste0("\n", " Naive agreement (", gagree, ") for cross validation informed XGBoost model : ",  "\n") )
        print( round( xgb.agree.naive , digits=digits) )
      }
    }                                                                             ## how to get rid of [1,]  ??
    
    ##### Random Forest ########################################################
    
    if (dorf == 1) { 
      if ((short==1) & (do_ncv == 1)) { 
        if (family %in% c("cox","binomial")) { 
          if (sum(ensemble[2:8])==0) { 
            cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed Random Forest :", 
                       round( rfAveAgree[1], digits = digits),"\n")) 
          } else {
            cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed Random Forest : \n") )
            print( round( rfAveAgree[1:2] , digits = digits)) 
          }
        } else { 
          if (sum(ensemble[2:8])==0) { 
            cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed Random Forest : ", 
                       round( rfAveAgree[1] , digits = digits), "\n") )
          } else {
            cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed Random Forest : \n") )
            print( round( rfAveAgree[1:3] , digits = digits)) 
          }
        } 
      }
      if ((short == 0) & (do_ncv == 1)) {
        if (sum(ensemble[2:8])==0) {
          cat(paste0("\n\n" , " Nested Cross Validation averages for Random Forest : ", "\n") )
          cat(paste0("      ", perunit, ": ", round( rfAveDevian[1] , digits = digits), "\n") )
          cat(paste0("      average number of variables random selected for the RF : ", round( rfAveMtry[1], digits=1) , "\n") )
          cat(paste0("      linear calibration coefficient :   ", round( rfAveLincal[1] , digits = digits) , "\n") )
          cat(paste0("      average agreement (", gagree, ") :  ", round( rfAveAgree[1] , digits = digits), "\n") )
        } else {
          cat(paste0("\n\n" , " Nested Cross Validation averages for Random Forest : ", "\n") )
          cat(paste0("\n" , "      ", perunit, ": \n")) 
          print( round( rfAveDevian , digits = digits) )
          cat(paste0("\n", "      average number of variables randomly selected for the RF : \n" ) )
          print( round( rfAveMtry  , digits=1) ) 
          cat(paste0("\n", "      linear calibration coefficient :   \n") )
          print( round( rfAveLincal , digits = digits) )
          cat(paste0("\n", "      average agreement (", gagree, ") :  \n") )
          print( round( rfAveAgree  , digits = digits) )
        }
      } 
      if (((short==1) & (do_ncv == 0)) | ((short==0) & (do_ncv == 1))) {
        cat(paste0("\n", " Naive Random Forest agreement (", gagree, ") :      \n") )
        print( round( rf.agree.naive , digits = digits) )
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
        cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed RPART : ",  "\n") )
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
        
        cat(paste0("\n", "      average agreement (", gagree, ") :         ", "\n") )
        print( round( rpartAveAgree , digits = digits) ) 
      }      
      if ((short == 0) | (do_ncv==0)) {
        cat(paste0("\n", " Cross validation informed RPART : ",  "\n") )
        cat(paste0("\n", "      number of terms used in model : ", "\n") )
        print(round( rpart.nzero, digits=1) ) 
        
        cat(paste0("\n", "      naive agreement (", gagree, ") : ",  "\n") )
        print( round( rpart.agree.naive , digits=digits) )
      }
    }      
    
    ##### ANN ##################################################################
    
    if (doann == 1) { 
      if (do_ncv == 1) {
        if (short==1) {
          cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed Neural Network : ",  "\n") )
          print( round( annAveAgree , digits = digits) ) 
        } else {
          cat(paste0("\n\n" , " Nested Cross Validation averages for neural network : ", "\n") )
        
          cat(paste0("\n" , "      ", perunit, ": ", "\n") )
          print( round( annAveDevian , digits = digits) )
        
          cat(paste0("\n", "      linear calibration coefficient : ", "\n") )
          print( round( annAveLincal , digits = digits) ) 
        
          cat(paste0("\n", "      average agreement (", gagree, ") :            ", "\n") )
          print( round( annAveAgree , digits = digits) ) 
        }
      }
      if ((short != 1) | (do_ncv == 0)) {
        cat(paste0("\n", " Cross validation informed neural network : ",  "\n") )
        cat(paste0("\n", "      naive agreement (", gagree, ") : ",  "\n") )
        print( round( ann.agree.naive, digits=digits) )
        cat(paste0("\n") )
      }
    }                                                                           
    
    ##### STEP #################################################################
      
    if ( (short==1) & (do_ncv == 1) ) {
      if ( (dostep==1) & (doaic==1) ) {
        cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed STEPWISE and ",  "\n") )
        cat(paste0( "   STEPWISE based AIC : ",  "\n") )
        print( round(StepAveAgree[c(1,2,3)], digits = digits) ) 
      } else if (dostep==1) {
        cat(paste0("\n", " Nested cross validation agreement (", gagree, ") for cross validation informed STEPWISE : ",  "\n") )
        print( round(StepAveAgree[c(1,2)], digits = digits) ) 
      } else if (doaic==1) {
        cat(paste0("\n", " Cross validation agreement (", gagree, ") for the STEPWISE based AIC : ",  "\n") )
        print( round(StepAveAgree[c(3)], digits = digits) )  
      }
    }
    
    if ( (short!=1) & (do_ncv == 1) & (dostep==1) ) {
#      cv.stepreg.fit.df = cv.stepreg.fit$best.df    #func.fit.df ; cv.stepreg.fit$best.p 
      cat(paste0("\n", " Nested Cross Validation STEPWISE regression model (df): ", "\n") )   
      cat(paste0("      Average linear calibration coefficient: ", round(StepAveLincal[1] ,  digits = 3), "\n") )    
      cat(paste0("      Average deviance : ", round(StepAveDevian[1]   ,  digits = digits), "\n") )    
      cat(paste0("      Average model df : ", round( StepAve_df[1], digits=2 ), "\n") )
      cat(paste0("      Average ", gagree, " : ", round(StepAveAgree[1],  digits = digits), "\n") )        
      cat(paste0(" Naive ", gagree, " based upon the same (all) data as model derivation (df): ",
                   round(cv.stepreg.fit$cvfit.df[6], digits=digits) , "\n") )
      cat(paste0("    Model df ", length(cv.stepreg.fit$func.fit.df$coefficients), "\n")) 
      
      cat(paste0("\n", " Nested Cross Validation STEPWISE regression model (p): ", "\n") )   
      cat(paste0("      Average linear calibration coefficient : ", round(StepAveLincal[1] ,  digits = 3), "\n") )    
      cat(paste0("      Average deviance : ", round(StepAveDevian[2] ,  digits = digits), "\n") )    
      cat(paste0("      Average model p  : ", round( StepAve_p [1], digits=3 ), "\n") )
      cat(paste0("      Average model df : ", round( StepAve_df[2], digits=2 ), "\n") )
      cat(paste0("      Average ", gagree, " : ", round(StepAveAgree[2],  digits = digits), "\n") )  
      cat(paste0(" Naive ", gagree, " based upon the same (all) data as model derivation (p): ",
                   round(cv.stepreg.fit$cvfit.p[6], digits=digits) , "\n") )
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
        cat(paste0("\n Naive ", gagree, " based upon the same (all) data as model derivation (df): ",
                   round(cv.stepreg.fit$cvfit.df[6]^pow, digits=digits) , "\n") )
        cat(paste0("    Model df ", length(cv.stepreg.fit$func.fit.df$coefficients), "\n")) 
        cat(paste0("\n Naive ", gagree, " based upon the same (all) data as model derivation (p): ",
                   round(cv.stepreg.fit$cvfit.p[6]^pow, digits=digits) , "\n") )
        cat(paste0("    Model df ", length(cv.stepreg.fit$func.fit.p$coefficients), "\n")) 
      }
    }
    
    if (doaic==1) {
      if ((short!=1) & (do_ncv==1)) {
        cat(paste0("\n Cross Validation results for STEPWISE regression model: (AIC)", "\n") ) 
        cat(paste0("      Average linear calibration coefficient : ", round(StepAveLincal[3] ,  digits = 3), "\n") )    
        cat(paste0("      Average deviance : ", round(StepAveDevian[3],  digits = digits), "\n") )    
        cat(paste0("      Average model df : ", round( StepAve_df[3], digits=2 ), "\n") )
        cat(paste0("      Average ", gagree, " : ", round( StepAveAgree[3], digits = digits), "\n") )   
      }
    
      if ((short!=1) | (do_ncv == 0))  {
        if (do_ncv == 0) { cat("\n") }
          cat(paste0(" Naive ", gagree, "based upon the same (all) data as model derivation (AIC) : ",         
                     round( object$func.fit.aic$aic.fit.n[6], digits=digits ), "\n") ) 
          #        round( 1-var(all.fit.aic$residuals)/var(all.fit.aic$y), digits=digits), "\n") ) 
        cat(paste0("    Model df ", length(object$func.fit.aic$coefficients), "\n")) 
      }
    }
    
    ##### END STEP #################################################################
    cat("\n")  
  } #### end of if summmary 
} 

###############################################################################################################
###############################################################################################################
################################################################################
################################################################################

#' Calculate agreement differences with CI and p
#' 
#' @description 
#' Perform a paired t-test as called from glmnetr.compcv().  
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
#' 
glmnetr.compcv0_0_4_2 = function(a, b, digits=4, txt=0, pow=1) { 
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
#' @param pow the power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  pow is ignored for the family of "cox" and "binomial".   
#' @return A printout to the R console. 
#' 
#' @seealso
#'   \code{\link{summary.nested.glmnetr}} 
#' 
#' @noRd
#' 
glmnetr.compcv_0_4_2 = function(object, digits=4, pow=1) {
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
  ensemble = object$ensemble 
  lasso.agree.cv  = object$lasso.agree.cv
  xgb.agree.cv    = object$xgb.agree.cv
  rf.agree.cv     = object$rf.agree.cv
  ann.agree.cv    = object$ann.agree.cv
  rpart.agree.cv  = object$rpart.agree.cv
  step.agree.cv   = object$step.agree.cv
  
#  pow = 1 
  
  if (family == "gaussian") { 
    if (pow == 2) { 
      pm = "R-square"
    } else { 
      pow = 1 
      pm = "Correlation" 
    }
  } else if (family %in% c("cox","binomial")) { 
    pow = 1 
    pm = "Concordance" 
  }
  
#  if ( pow == 2 ) { gagree = "R-square" } else { gagree = "Correlation" }
#  if (family %in% c("cox","binomial")) { gagree = "Concordance" ; pow = 1 }
  
  if (sum(ensemble[c(2:8)]) > 0 ) {
    cat ("\n  Ensemble paramter used when fitting models : \n         ensemble\n" ) 
    cat(paste0("    (", ensemble[1],",", ensemble[2],",", ensemble[3],",", ensemble[4],", ",
                 ensemble[5],",", ensemble[6],",", ensemble[7],",", ensemble[8],")\n\n")) 
  } 
  
  if ( ensemble[c(1)] == 0 ) {
    cat ("\n Simple models with informaiton form loass fot not run.  Output is abbreviated. \n\n" ) 
    doxgb = 0 ; dorf = 0 ; doann = 0 ; dorpart = 0 ; dostep = 0 ; doaic = 0 ; 
  }
  
  cat ("  Model performance comparison in terms of ", pm, "\n\n" )   
  cat ("  Comparison                                estimate   (95% CI)         p\n") 


  if (dolasso == 1) {
    cat ("\n lasso.minR  - lasso.min                     ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , lasso.agree.cv[,2],pow=pow) 
    cat ("\n lasso.minR  - lasso.minR0                   ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , lasso.agree.cv[,6],pow=pow)   
    cat ("\n lasso.min   - lasso.minR0                   ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,2] , lasso.agree.cv[,6],pow=pow)   
    cat("\n")
  }
  
#  print(xgb.agree.cv)

  if (doxgb == 1) {
    cat ("\n XGBoost (tuned) - XGBoost (simple)          ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,4] , xgb.agree.cv[,1],pow=pow) ; 
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n XGBoost (tuned) lasso feature - no feature  ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,5] , xgb.agree.cv[,4],pow=pow) ;  
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      cat ("\n XGBoost (tuned) lasso offset - no offset    ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,6] , xgb.agree.cv[,4],pow=pow) ;  
    }
    cat("\n")
  }
  
  if (dorf == 1) {
    lr = 0 
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n RF with lasso feature - no feature          ") ;  glmnetr.compcv0_0_4_2(rf.agree.cv[,2] , rf.agree.cv[,1],pow=pow) ;  
      lr = 1 
    }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) {
      cat ("\n RF with lasso offset - no offset            ") ;  glmnetr.compcv0_0_4_2(rf.agree.cv[,3] , rf.agree.cv[,1],pow=pow) ;  
      lr = 1 
    }
    if (lr == 1) { cat("\n") } 
  }
  
  if (doann == 1) {
    lr = 0 
    if (sum(ensemble[6])> 0) {
      cat ("\n ANN with with lasso feature - no feature    ") ;  glmnetr.compcv0_0_4_2(ann.agree.cv[,6] , ann.agree.cv[,1],pow=pow) ; lr = 1 
    } else if (sum(ensemble[2])> 0) {
      cat ("\n ANN with with lasso feature - no feature    ") ;  glmnetr.compcv0_0_4_2(ann.agree.cv[,2] , ann.agree.cv[,1],pow=pow) ; lr = 1 

    } 
    if (sum(ensemble[8])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  glmnetr.compcv0_0_4_2(ann.agree.cv[,8] , ann.agree.cv[,1],pow=pow) ; lr = 1 
    } else if (sum(ensemble[7])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  glmnetr.compcv0_0_4_2(ann.agree.cv[,7] , ann.agree.cv[,1],pow=pow) ; lr = 1  
    } else     if (sum(ensemble[4])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  glmnetr.compcv0_0_4_2(ann.agree.cv[,4] , ann.agree.cv[,1],pow=pow) ; lr = 1 
    } else     if (sum(ensemble[3])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  glmnetr.compcv0_0_4_2(ann.agree.cv[,3] , ann.agree.cv[,1],pow=pow) ; lr = 1 
    } 
    if (lr == 1) { cat("\n") } 
  }
  
  if (dostep == 1) {
    cat ("\n step (df) - step (p)                        ") ;  glmnetr.compcv0_0_4_2(step.agree.cv[,1]      , step.agree.cv[,2],pow=pow)    ;  cat("\n")
  }
  
  cat("\n")

  if ((dolasso == 1) & (doxgb == 1)) {
    cat ("\n lasso.minR - XGB (tuned)                    ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , xgb.agree.cv[,4],pow=pow) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - XGB with lasso feature         ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , xgb.agree.cv[,5],pow=pow) 
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      cat ("\n lasso.minR - XGB with lasso offset          ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , xgb.agree.cv[,6],pow=pow)   
    }
  }

  if ((dolasso == 1) & (dorf == 1)) {
    cat ("\n lasso.minR - Random Forest                  ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , rf.agree.cv[,1],pow=pow) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - RF with lasso feature          ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , rf.agree.cv[,2],pow=pow) 
    }
    if ( (sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian") ) {
      cat ("\n lasso.minR - RF with lasso offset           ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , xgb.agree.cv[,3],pow=pow)   
    }
  }

  if ((dolasso == 1) & (doann == 1)) {
    cat ("\n lasso.minR - ANN                            ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , ann.agree.cv[,1],pow=pow) 
    if (ensemble[6]) { 
      cat ("\n lasso.minR - ANN l lasso feature            ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , ann.agree.cv[,6],pow=pow)   
    } else if (ensemble[2]) { 
      cat ("\n lasso.minR - ANN lasso feature              ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] ,  ann.agree.cv[,2],pow=pow)   
    }  
    if (ensemble[8]) { 
      cat ("\n lasso.minR - ANN l lasso offset (upated)    ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , ann.agree.cv[,8],pow=pow)   
    } else if (ensemble[4]) { 
      cat ("\n lasso.minR - ANN lasso offset               ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , ann.agree.cv[,4],pow=pow)  
    } else if (ensemble[7]) { 
      cat ("\n lasso.minR - ANN l lasso offset (upated)    ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , ann.agree.cv[,7],pow=pow)   
    }  else if (ensemble[3]) { 
      cat ("\n lasso.minR - ANN lasso offset               ") ;  glmnetr.compcv0_0_4_2(lasso.agree.cv[,4] , ann.agree.cv[,3],pow=pow)   
    }
  }

  if (dolasso) { cat("\n") } 

  if ((doxgb == 1) & (dorf == 1)) {
    cat ("\n XGBoost (tuned) - RF                        ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,4] ,  rf.agree.cv[,1],pow=pow)   
    if (sum(ensemble[c(2,6)]) > 0) {
      cat ("\n XGBoost - RF with lasso feature             ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,5] ,  rf.agree.cv[,2],pow=pow)   
    } 
    if ( (sum(ensemble[c(3,4,7,8)]) > 0)  & (family == "gaussian") ) {
      cat ("\n XGBoost - RF with lasso offset              ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,6] ,  rf.agree.cv[,3],pow=pow)   
    }
  }

  if ((doxgb == 1) & (doann == 1)) {
    cat ("\n XGBoost (tuned) - ANN                       ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,4] ,  ann.agree.cv[,1],pow=pow)   
    if (ensemble[6]) { 
      cat ("\n XGBoost - ANN, lasso feature                ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,5] ,  ann.agree.cv[,6],pow=pow)   
    } else if (ensemble[2]) { 
      cat ("\n XGBoost - ANN, lasso feature                ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,5] ,  ann.agree.cv[,2],pow=pow)   
    } 
    if (family == "gaussian") {
      if (ensemble[8]) { 
        cat ("\n XGBoost - ANN, l lasso offset (upated)      ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,6] ,  ann.agree.cv[,8],pow=pow)   
      } else if (ensemble[4]) { 
        cat ("\n XGBoost - ANN, lasso offset                 ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,6] ,  ann.agree.cv[,4],pow=pow)   
      } else if (ensemble[7]) { 
        cat ("\n XGBoost - ANN, l lasso offset (updated)     ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,6] ,  ann.agree.cv[,7],pow=pow)   
      } else if (ensemble[3]) { 
        cat ("\n XGBoost - ANN, lasso offset                 ") ;  glmnetr.compcv0_0_4_2(xgb.agree.cv[,6] ,  ann.agree.cv[,3],pow=pow)   
      }  
    }
  }
  
  if (doxgb) { cat("\n") }
  
  if ((dorf == 1) & (doann == 1)) {
    cat ("\n RF - ANN,                                   ") ;  glmnetr.compcv0_0_4_2(rf.agree.cv[,1] ,  ann.agree.cv[,1],pow=pow) 
    if (ensemble[6]) {
      cat ("\n RF - ANN, l lasso feature                   " ) ;  glmnetr.compcv0_0_4_2(rf.agree.cv[,2] ,  ann.agree.cv[,6],pow=pow) 
    } else if (ensemble[2]) {
      cat ("\n RF - ANN, lasso feature                     " ) ;  glmnetr.compcv0_0_4_2(rf.agree.cv[,2] ,  ann.agree.cv[,2],pow=pow)  
    }
    if (ensemble[8]) {
      cat ("\n RF - ANN, l lasso offset (upated)           " ) ;  glmnetr.compcv0_0_4_2(rf.agree.cv[,3] ,  ann.agree.cv[,8],pow=pow)   
    } else if (ensemble[4]) {
      cat ("\n RF - ANN, lasso offset                      " ) ;  glmnetr.compcv0_0_4_2(rf.agree.cv[,3] ,  ann.agree.cv[,4],pow=pow)  
    } else if (ensemble[7]) {
      cat ("\n RF - ANN, l lasso offset (upated)           " ) ;  glmnetr.compcv0_0_4_2(rf.agree.cv[,3] ,  ann.agree.cv[,7],pow=pow) 
    } else if (ensemble[3]) {
      cat ("\n RF - ANN, lasso offset                      " ) ;  glmnetr.compcv0_0_4_2(rf.agree.cv[,3] ,  ann.agree.cv[,3],pow=pow)  
    }
    cat("\n")
  }
  
  cat("\n")
}

####################################################################################################################################
####################################################################################################################################

