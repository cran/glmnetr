###### 240215 ######################################################################################################################
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
  summary(x) 
}
#
################################################################################
#' round elements of a summary.glmnetr() output
#'
#' @param summdf a summary data frame from summary.nested.glmnetr() obtained using 
#' the option table=0 
#' @param digits the minimum number of decimals to display the elements of the data
#' frame 
#' @param do_ncv 1 (default) if the summdf object is a summary for an analysis including
#' nested cross validation, 0 if only the full data models were fit.
#'
#' @return a data frame with same form as the input but with rounding for easier 
#' display
#' 
#' @export 
#'
roundperf = function(summdf, digits=3, do_ncv=1) {
  if (do_ncv == 1) {
    set1 = c(1,2,5,6,8) 
    set2 = c(3) 
    set3 = c(4:7)
  } else {
    set1 = c(1,3) 
    set2 = NULL 
    set3 = c(2)
  }
  for ( j_ in set1 ) {
    digitst = digits 
    for (i_ in c(0:3)) {
      if (min(abs(summdf[j_])) < 2*10^(-digits-i_)) { digitst = digits + i_ + 1 }
    }
    #        print(c(j_,digitst))
    summdf[j_] = round(summdf[j_], digits=digitst)
  }
  for ( j_ in set2 ) {
    digitst = digits 
    for (i_ in c(0:3)) {
      #          print( summdf[j_] )
      #          print( summdf[j_] - 1)
      if (min(abs(summdf[j_][!is.na(summdf[j_])]-1)) < 2*10^(-digits-i_)) { digitst = digits + i_ + 1 }
      #          print(c(i_,digitst))
    }
    
    summdf[j_] = round(summdf[j_], digits=digitst)
  }
  for ( j_ in set3 ) {
    digitst = digits 
    for (i_ in c(0:3)) {
      if (max(abs(summdf[j_])) > (1 -2*10^(-digits-i_))) { digitst = digits + i_ + 1 }
    }
    summdf[j_] = round(summdf[j_], digits=digitst)
  }
  return( summdf )
}

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
#' @param pow the power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  Pow is ignored for the family of "cox" and "binomial".   
#' @param printg1 TRUE to also print out the fully penalized lasso beta, else to suppress.  
#' Only applies to cvfit=TRUE.
#' @param digits digits for printing of deviances, linear calibration coefficients 
#' and agreement (concordances and R-squares).
#' @param Call 1 to print call used in generation of the object, 0 or NULL to not print 
#' @param onese 0 (default) to not include summary for 1se lasso fits in tables, 1 to include 
#' @param table 1 to print table to console, 0 to output the tabled information to a data frame
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
# cvfit = FALSE ; pow=2 ; printg1 = FALSE ; digits = 3 ; Call=NULL ; onese = 0 ; table = 1 ; 
#'
summary.nested.glmnetr = function(object, cvfit=FALSE, pow=2, printg1=FALSE, 
                                  digits=4, Call=NULL, onese=0, table=1, ...) {
  
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
#    sample[6] = round(as.numeric(sample[6]), digits=digits)
#    if (family=="cox") { names(sample)[6]="null.dev/events" 
#    } else { names(sample)[6] = "null.dev/obs"  }
    
    if ( pow == 2 ) { gagree = "R-square" } else { gagree = "Correlation" }
    if (family %in% c("cox","binomial")) { gagree = "Concordance" ; pow = 1 } 

    annrownames = c("ANN Uninformed", "ANN lasso features", "ANN lasso weights", "ANN lasso weights, reset",
                    "ANN only lasso terms", "ANN lasso terms, lasso features", "ANN lasso terms, lasso weights", "ANN lasso terms, weights reset") 

    colnames1 = c("Ave DevRat", "Ave Int", "Ave Slope", paste("Ave", gagree), "Ave Non Zero", "Naive Devian", paste("Naive", gagree), "Non Zero" )
    colnames0 = c( "Naive Devian", paste("Naive", gagree), "Non Zero" )

    if (doann == 1) { 
      if        (ensemble[4]==1) { ann_cv = object$ann_fit_4 ; whichann = annrownames[4] ;
      } else if (ensemble[8]==1) { ann_cv = object$ann_fit_8 ; whichann = annrownames[8] ;
      } else if (ensemble[3]==1) { ann_cv = object$ann_fit_3 ; whichann = annrownames[3] ;
      } else if (ensemble[7]==1) { ann_cv = object$ann_fit_7 ; whichann = annrownames[7] ;
      } else if (ensemble[2]==1) { ann_cv = object$ann_fit_2 ; whichann = annrownames[2] ;
      } else if (ensemble[6]==1) { ann_cv = object$ann_fit_6 ; whichann = annrownames[6] ;
      } else if (ensemble[1]==1) { ann_cv = object$ann_fit_1 ; whichann = annrownames[1] ;
      } else if (ensemble[5]==1) { ann_cv = object$ann_fit_5 ; whichann = annrownames[5] ;
      }
    }
    
    if (dostep==1) { cv.stepreg.fit    = object$cv.stepreg.fit }  
    if (doaic ==1) { func.fit.aic      = object$func.fit.aic   }

    if ( table %in% c(1,2) ) { 
      if (!is.null(Call)) { 
        cat(paste0("\n     function call :\n\n"))  
        print(Call) 
        #          cat(paste0("\n")) 
      } 
      
      cat(paste0("\n"  , " Sample information including number of records, "))     
      if (family %in% c("cox","binomial")) { cat(paste0("events, ")) }
      cat(paste0( "number of columns in", "\n", " design (predictor, X) matrix, and df (rank) of design matrix: ", "\n") )
#      if ((family %in% c("cox","binomial"))==0) { sample = sample[c(-3)] }

      x_1 = as.numeric(sample[6:7])
      if (min(x_1) > 0.2) { x_1 = round(x_1,digits=2)
      } else if (min(x_1) > 0.02 ) { x_1 = round(x_1,digits=3)
      } else if (min(x_1) > 0.002) { x_1 = round(x_1,digits=4)
      }
      sample[6:7] = x_1 
#      print(sample)
      x_1 = as.numeric(sample[8])
      if ((x_1 > 0)==1) {
        if (x_1 > 0.2) { x_1 = round(x_1,digits=2)
        } else if (x_1 > 0.02 ) { x_1 = round(x_1,digits=3)
        } else if (x_1 > 0.002) { x_1 = round(x_1,digits=4)
        }
      }
      sample[8] = x_1 
            
      if ((family %in% c("cox","binomial"))==1) { print(sample,quote=FALSE)
      } else { print(sample[c(-3)],quote=FALSE)  }
      
      x_ = sample
      nms_ = names(x_)
      x_1 = matrix(as.numeric(x_[-1]), nrow=length(x_)-1,ncol=1 )
#      if (min(x_1) > 0.1) { x_1 = round(x_1,digits=2)
#      } else if (min(x_1) > 0.01 ) { x_1 = round(x_1,digits=3)
#      } else if (min(x_1) > 0.001) { x_1 = round(x_1,digits=4)
#      }
      rownames(x_1) = paste0("  ", nms_[-1])
      colnames(x_1) = ""
#      cat(paste0("\n  ", paste(nms_[1], "               ", x_[1],"")))
#      print(x_1)
      
      if (!is.null(object$dep_names)) { 
        cat(paste0("\n"  , " Dependent Variable(s) : ", "\n") )
        print(object$dep_names) 
      }
      
      if ((dolasso ==0) & (dostep==0) & (doaic==0)) { 
        cat(paste0("\n"  , " Tuning parameters for models : ", "\n") )
        print(object$tuning[c(1,2)], quote=FALSE) 
      } else if ((dostep==0) & (doaic==0)) { 
        cat(paste0("\n"  , " Tuning parameters for models : ", "\n") )
        print(object$tuning[c(1:5)], quote=FALSE) 
      } else { 
        cat(paste0("\n"  , " Tuning parameters for models : ", "\n") )
        print(object$tuning, quote=FALSE) 
      }
      
      if (doann  == 1) {
        cat(paste0("\n"  , " Tuning parameters for  ", whichann, "  model : ", "\n") )
#        print(ann_cv$modelsum[c(1:9,11:12)]) 
        class( ann_cv$modelsum[c(1:9,11:12)] )
        x = ann_cv$modelsum[c(1:9,11:12)]
        namesx = names(x)
        mx = matrix(x,nrow=1,ncol=11) 
        colnames(mx) = namesx
        print(data.frame(mx)[1,c(1:11)], row.names=FALSE)
      }
    }
    
#    if (doaic==1) {
#      cat(paste0("\n", " Average deviance for null model", "\n") )    ## pull other data applicable to all models 
#      print( round(StepAveDevian[1], digits = digits ) ) 
#    }
    
#    lassot = lasso 
#    roundperf(lasso,3)

    if (family == "cox") { perunit = "deviance per event " } else { perunit = "deviance per record " }
    
    get.DevRat = function( devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv ) {
      ncoldevian = dim(devian.cv)[2] 
      AveDevRat = rep(1,ncoldevian)
      for (j_ in c(1:ncoldevian)) {
        AveDevRat[j_] = devrat_(devian.cv[,j_], null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )[[2]]
      }
      return( AveDevRat )
    }
    
    null.m2LogLik.cv = object$null.m2LogLik.cv
    sat.m2LogLik.cv = object$sat.m2LogLik.cv
    n.cv = object$n.cv 

    ## PERFORMANCE SUMMARIES ###################################################
    ## PERFORMANCE SUMMARIES ###################################################
    
    lasso = NULL 
    if (dolasso == 1) {
      lassoAveDevRat = get.DevRat( object$lasso.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
      lassoAveDevian = colMeans(object$lasso.devian.cv, na.rm=TRUE)  
      lassoAveIntcal = colMeans(object$lasso.intcal.cv, na.rm=TRUE)  
      lassoAveLincal = colMeans(object$lasso.lincal.cv, na.rm=TRUE)  
      lassoAveAgree  = colMeans(object$lasso.agree.cv, na.rm=TRUE)
      lassoAveNzero  = colMeans(object$lasso.nzero.cv, na.rm=TRUE)  
      lasso.agree.naive = object$lasso.agree.naive 
      if (family == "gaussian") { 
        lasso.agree.naive = lasso.agree.naive ^pow
        lassoAveAgree = lassoAveAgree ^pow
      }
      ## sqrt( apply(lasso.agree.cv,2,var) ) 

      if (do_ncv == 1) {
        if ( table %in% c(1,2) ) { cat( "\n  LASSO: Ave is for (nested) CV model performance summary, else naive summary for \n",
                                   "         fit on all data \n" ) }
        lasso = data.frame( lassoAveDevRat , lassoAveIntcal, lassoAveLincal , lassoAveAgree, lassoAveNzero, 
                            object$lasso.devian.naive, object$lasso.agree.naive, object$lasso.nzero)
        names(lasso) = c("Ave DevRat" , "Ave Int", "Ave Slope" , paste0("Ave ", gagree) , "Ave Non Zero", "Naive Devian", paste0("Naive ", gagree) , "Non Zero" )
      } else {
        if ( table %in% c(1,2) ) { cat( "\n  LASSO: Naive summary for fit on all data \n" ) }
        
        lasso = data.frame( object$lasso.devian.naive, object$lasso.agree.naive, object$lasso.nzero)
        names(lasso) = c( "Naive Devian", paste0("Naive ", gagree) , "Non Zero" )
      }
      rownames = paste0(c(rep("LASSO ",6),""), row.names(lasso)) 
      rownames[7] = "Ridge" 
      row.names(lasso) = rownames 
      
      if (onese == 0) { lasso = lasso[c(2,4,6,7),] } 
      lassor = roundperf(lasso, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { lassor = lassor[,-2] }      
      if ( table %in% c(1,2) ) { print( lassor ) }
    }        
      
    ## XGB #####################################################################
    
    xgb = NULL 
    if (doxgb == 1) { 
      en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
      en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
      en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
      enx = c(en1, en2, en3, en1, en2, en3) 
      xgbAveDevRat = get.DevRat( object$xgb.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
      xgbAveDevian = colMeans( object$xgb.devian.cv, na.rm=TRUE )
      xgbAveIntcal = colMeans( object$xgb.intcal.cv, na.rm=TRUE )
      xgbAveLincal = colMeans( object$xgb.lincal.cv, na.rm=TRUE )
      xgbAveAgree  = colMeans( object$xgb.agree.cv, na.rm=TRUE ) 
      xgbAveNzero  = colMeans( object$xgb.nzero.cv, na.rm=TRUE ) 
      xgb.agree.naive = object$xgb.agree.naive
      if (family == "gaussian") { 
        xgb.agree.naive = xgb.agree.naive ^pow 
        xgbAveAgree = xgbAveAgree ^pow 
      }
        
      if (do_ncv == 1) {
        if ( table %in% c(1,2) ) { cat( "\n  XGBoost: Ave is for (nested) CV model performance summary, else naive summary for\n",
                                   "           fit on all data \n" ) }
        xgb = data.frame( xgbAveDevRat , xgbAveIntcal, xgbAveLincal , xgbAveAgree, xgbAveNzero, object$xgb.devian.naive, xgb.agree.naive, object$xgb.nzero )
        names( xgb ) = colnames1 
        xgb
      } else {
        if ( table %in% c(1,2) ) { cat( "\n  XGBoost: Naive model performance summary for fit on all data\n" ) }
        xgb = data.frame( object$xgb.devian.naive, xgb.agree.naive , object$xgb.nzero )
        names(xgb) = c( "Naive Devian", paste("Naive", gagree), "Non Zero" )
      }
      row.names(xgb) = c("XGB Simple", "XGB Feature", "XGB Offset", "XGB Tuned", "XGB Tuned Feature", "XGB Tuned Offset" )
      xgb = xgb[enx==1,]
      xgbr = roundperf(xgb, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { xgbr = xgbr[,-2] }      
      if ( table %in% c(1,2) ) { print( xgbr ) }
    }   
    
    ##### Random Forest ########################################################
    
    rf = NULL 
    if (dorf == 1) { 
      en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
      en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
      en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
      if (family != "gaussian") { en3 = 0 }
      enx = c(en1, en2, en3) 
      rfAveDevRat = get.DevRat( object$rf.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
      rfAveDevian = colMeans(object$rf.devian.cv, na.rm=TRUE) # [enx==1]   
      rfAveIntcal = colMeans(object$rf.intcal.cv, na.rm=TRUE)
      rfAveLincal = colMeans(object$rf.lincal.cv, na.rm=TRUE)
      rfAveAgree  = colMeans(object$rf.agree.cv, na.rm=TRUE)
      rfAveMtry   = colMeans(object$rf.mtry.cv, na.rm=TRUE)
      rf.agree.naive = object$rf.agree.naive
      if (family == "gaussian") { 
        rf.agree.naive = rf.agree.naive ^pow 
        rfAveAgree  = rfAveAgree ^pow 
      }
      
      if (do_ncv == 1) {
        if ( table %in% c(1,2) ) { cat( "\n  Random Forest: Ave is for CV model performance summary, else naive\n",
                                 "                   summary for fit on all data \n" ) }
        
        rf = data.frame( rfAveDevRat , rfAveIntcal , rfAveLincal , rfAveAgree, rfAveMtry, object$rf.devian.naive, rf.agree.naive, object$rf.mtry  )
        names( rf ) = colnames1 
        rf
      } else {
        if ( table %in% c(1,2) ) { cat( "\n  Random Forest: Naive model performance summary for fit on all data\n" ) }
        rf = data.frame( object$rf.devian.naive, rf.agree.naive , object$rf.mtry )
        names(rf) = c( "Naive Devian", paste("Naive", gagree), "Non Zero" )
      }
      row.names(rf) = c("RF Simple", "RF lasso Feature", "RF lasso Offset" )
      rf = rf[enx==1,]
      rfr = roundperf(rf, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { rfr = rfr[,-2] }      
      if ( table %in% c(1,2) ) { print( rfr ) }
    }

    ##### ANN ##################################################################
    
    ann = NULL 
    if (doann == 1) { 
      annAveDevRat = get.DevRat( object$ann.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
      annAveDevian = colMeans(object$ann.devian.cv, na.rm=TRUE)               # [(ensemble==1)]
      annAveIntcal = colMeans(object$ann.intcal.cv, na.rm=TRUE)
      annAveLincal = colMeans(object$ann.lincal.cv, na.rm=TRUE)
      annAveAgree  = colMeans(object$ann.agree.cv , na.rm=TRUE)
      annAveNzero  = colMeans(object$ann.nzero.cv , na.rm=TRUE)
      ann.agree.naive = object$ann.agree.naive
      if (family == "gaussian") { 
        ann.agree.naive = ann.agree.naive ^pow
        annAveAgree = annAveAgree ^pow 
      }
      
      if (do_ncv == 1) {
        if ( table %in% c(1,2) ) { cat( "\n  Artificial Neural Network: Ave is for (nested) CV model performance summary,\n",
                                        "       else naive summary for fit on all data \n" ) }
        ann = data.frame( annAveDevRat , annAveIntcal , annAveLincal , annAveAgree, annAveNzero, object$ann.devian.naive, ann.agree.naive, object$ann.nzero )
        names( ann ) = colnames1 
        ann
      } else {
        if ( table %in% c(1,2) ) { cat( "\n  Artificial Neural Network: Naive model performance summary for fit on all data\n" ) }
        ann = data.frame( object$ann.devian.naive, ann.agree.naive , object$ann.nzero )
        names(ann) = c( "Naive Devian", paste("Naive", gagree), "Non Zero" )
        ann
      }
      row.names(ann)  = annrownames 
      ann = ann[ensemble[c(1:8)]==1,]
      annr = roundperf(ann, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { annr = annr[,-2] }      
      if ( table %in% c(1,2) ) { print( annr) }
    }                                                                           
    
    ##### RPART ################################################################
    
    rpart = NULL 
    if (dorpart == 1) { 
        en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
        en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
        en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
        if (family == "binomial") { en3 = 0 }
        enx = c(en1,en1,en1, en2,en2,en2, en3,en3,en3) 
        rpartAveDevRat = get.DevRat( object$rpart.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
        rpartAveDevian = colMeans(object$rpart.devian.cv, na.rm=TRUE) 
        rpartAveIntcal = colMeans(object$rpart.intcal.cv, na.rm=TRUE) 
        rpartAveLincal = colMeans(object$rpart.lincal.cv, na.rm=TRUE) 
        rpartAveAgree  = colMeans(object$rpart.agree.cv, na.rm=TRUE) 
        rpartAveNzero  = colMeans(object$rpart.nzero.cv, na.rm=TRUE)                          # [c(3,2,1,6,5,4,9,8,7)]          # [enx==1]
        rpart.agree.naive = object$rpart.agree.naive
        if (family == "gaussian") { 
          rpart.agree.naive = rpart.agree.naive ^pow 
          rpartAveAgree = rpartAveAgree ^pow 
        }
        if (do_ncv == 1) {
          if ( table %in% c(1,2) ) { cat( "\n  Recursive Partitioning: Ave is for CV model performance summary, else\n", 
                                     "                          naive summary for fit on all data \n" ) }
          rpart = data.frame( rpartAveDevRat , rpartAveIntcal , rpartAveLincal , rpartAveAgree, rpartAveNzero, object$rpart.devian.naive, rpart.agree.naive, object$rpart.nzero )
          #        rpart = rpart[c(3,2,1,6,5,4,9,8,7),]
          names( rpart ) = colnames1 
          rpart
        } else {
          if ( table %in% c(1,2) ) { cat( "\n  Recursive Partitioning: Naive model performance summary for fit on all data\n" ) }
          rpart = data.frame( object$rpart.devian.naive, rpart.agree.naive , object$rpart.nzero )
          names(rpart) = c( "Naive Devian", paste("Naive", gagree), "Non Zero" )
        }
        rownames = row.names(rpart)
        rownames[c(4:9)] = substr(rownames[c(4:9)], 3,10)
#        rownames = paste(c("RPART", "RPART", "RPART", "RPART lasso feature", "RPART lasso feature", "RPART lasso feature",
#                           "RPART lasso offset", "RPART lasso offset", "RPART lasso offset"), rownames)
        rownames = paste(c("RPART", "RPART", "RPART", "RPART l. feat", "RPART l. feat", "RPART l. feat",
                           "RPART l. offs", "RPART l. offs", "RPART l. offs"), rownames)
        if ( table %in% c(1,2) ) { 
          if ((en2 == 1) | (en3 == 1)) { cat( "          l. feat for lasso included as feature, l. offs for lasso included as offset\n" )
          } else if (en2 == 1) { cat( "                  l. feat for lasso included as feature\n" )
          } else if (en3 == 1) { cat( "                  l. offs for lasso included as offset\n" ) }
        }
        row.names(rpart) = rownames 
        rpart = rpart[enx==1,]
        rpartr = roundperf(rpart, digits, do_ncv) 
        if ((family == "cox") & (do_ncv==1)) { rpartr = rpartr[,-2] }      
        if ( table %in% c(1,2) ) { print( rpartr ) } 
    }      
    
    ##### STEPWISE p or df tuned & AIC #########################################

    step = NULL 
    if ((dostep==1) | (doaic==1)) {
      StepAveDevRat = get.DevRat( object$step.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
      StepAveDevian = colMeans( object$step.devian.cv, na.rm=TRUE)
      StepAveIntcal = colMeans( object$step.intcal.cv, na.rm=TRUE)
      StepAveLincal = colMeans( object$step.lincal.cv, na.rm=TRUE)
      StepAveAgree  = colMeans( object$step.agree.cv , na.rm=TRUE)
      StepAve.nzero = colMeans( object$step.nzero.cv , na.rm=TRUE)
      StepAve.p     = colMeans( object$step.p.cv     , na.rm=TRUE)
      step.agree.naive = object$step.agree.naive 
      if (family == "gaussian") { 
        StepAveAgree = StepAveAgree ^pow 
      }
      
      if (do_ncv == 1) {
        if ( table %in% c(1,2) ) {
          if        ((dostep == 1) & (doaic == 1)) { 
            cat( "\n  Stepwise tuned and AIC: Ave is for (nested) CV model performance summary, else\n")
          } else if ((dostep == 1) & (doaic == 0)) { 
            cat( "\n  Stepwise tuned: Ave is for (nested) CV model performance summary, else\n")
          } else if ((dostep == 0) & (doaic == 1)) { 
            cat( "\n  Stepwise AIC: Ave is for (nested) CV model performance summary, else\n")
          }
          cat( "       naive summary for fit on all data \n" ) 
        }
        step = data.frame( StepAveDevRat , StepAveIntcal , StepAveLincal , StepAveAgree, StepAve.nzero, object$step.devian.naive, step.agree.naive, object$step.nzero )
         names( step ) = colnames1 
#        names( step ) = c("Deviance", "Inter", "Slope", "Concor", "Non Zero", "Deviancd", "ConroR", "Non Zero" )
      } else {
        if ( table %in% c(1,2) ) { cat( "\n  Stepwise tuned and AIC: Naive model performance summary for fit on all data\n" ) }
        step = data.frame( object$step.devian.naive, step.agree.naive , object$nzero )
        names(step) = colnames0 
      }
      row.names(step) = c("Stepwise df tuned", "Stepwise p tuned", "Stepwise AIC selected" )
      if        ((dostep == 1) & (doaic == 0)) { 
        step = step[c(1,2),] 
      } else if ((dostep == 0) & (doaic == 1)) { 
        step = step[c(3),]  
      }
      stepr = roundperf(step, digits, do_ncv)
      if ((family == "cox") & (do_ncv==1)) { stepr = stepr[,-2] }      
      if ( table %in% c(1,2) ) { print( stepr ) }
    }
    cat("\n")  
#    (names(lasso) == names(ann))*1 
#    if ( table %in% c(0) ) { return( rbind(lasso, xgb, rf, rpart, ann, step) ) }
    if ( table %in% c(0,2) ) { return( rbind(lasso, xgb, rf, ann, rpart, step) ) }
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
glmnetr.compcv0 = function(a, b, digits=4, txt=0, pow=1) { 
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
glmnetr.compcv = function(object, digits=4, pow=1) {
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
      pm = " R-square  "
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
    cat ("\n lasso.minR  - lasso.min                     ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , lasso.agree.cv[,2],pow=pow) 
    cat ("\n lasso.minR  - lasso.minR0                   ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , lasso.agree.cv[,6],pow=pow)   
    cat ("\n lasso.min   - lasso.minR0                   ") ;  glmnetr.compcv0(lasso.agree.cv[,2] , lasso.agree.cv[,6],pow=pow)   
    cat("\n")
  }
  
#  print(xgb.agree.cv)

  if (doxgb == 1) {
    cat ("\n XGBoost (tuned) - XGBoost (simple)          ") ;  glmnetr.compcv0(xgb.agree.cv[,4] , xgb.agree.cv[,1],pow=pow) ; 
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n XGBoost (tuned) lasso feature - no feature  ") ;  glmnetr.compcv0(xgb.agree.cv[,5] , xgb.agree.cv[,4],pow=pow) ;  
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      cat ("\n XGBoost (tuned) lasso offset - no offset    ") ;  glmnetr.compcv0(xgb.agree.cv[,6] , xgb.agree.cv[,4],pow=pow) ;  
    }
    cat("\n")
  }
  
  if (dorf == 1) {
    lr = 0 
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n RF with lasso feature - no feature          ") ;  glmnetr.compcv0(rf.agree.cv[,2] , rf.agree.cv[,1],pow=pow) ;  
      lr = 1 
    }
    if ((sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian")) {
      cat ("\n RF with lasso offset - no offset            ") ;  glmnetr.compcv0(rf.agree.cv[,3] , rf.agree.cv[,1],pow=pow) ;  
      lr = 1 
    }
    if (lr == 1) { cat("\n") } 
  }
  
  if (doann == 1) {
    lr = 0 
    if (sum(ensemble[6])> 0) {
      cat ("\n ANN with with lasso feature - no feature    ") ;  glmnetr.compcv0(ann.agree.cv[,6] , ann.agree.cv[,1],pow=pow) ; lr = 1 
    } else if (sum(ensemble[2])> 0) {
      cat ("\n ANN with with lasso feature - no feature    ") ;  glmnetr.compcv0(ann.agree.cv[,2] , ann.agree.cv[,1],pow=pow) ; lr = 1 

    } 
    if (sum(ensemble[8])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  glmnetr.compcv0(ann.agree.cv[,8] , ann.agree.cv[,1],pow=pow) ; lr = 1 
    } else if (sum(ensemble[7])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  glmnetr.compcv0(ann.agree.cv[,7] , ann.agree.cv[,1],pow=pow) ; lr = 1  
    } else     if (sum(ensemble[4])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  glmnetr.compcv0(ann.agree.cv[,4] , ann.agree.cv[,1],pow=pow) ; lr = 1 
    } else     if (sum(ensemble[3])> 0) { 
      cat ("\n ANN with with lasso offset - no offset      ") ;  glmnetr.compcv0(ann.agree.cv[,3] , ann.agree.cv[,1],pow=pow) ; lr = 1 
    } 
    if (lr == 1) { cat("\n") } 
  }
  
  if (dostep == 1) {
    cat ("\n step (df) - step (p)                        ") ;  glmnetr.compcv0(step.agree.cv[,1]      , step.agree.cv[,2],pow=pow)    ;  cat("\n")
  }
  
  cat("\n")

  if ((dolasso == 1) & (doxgb == 1)) {
    cat ("\n lasso.minR - XGB (tuned)                    ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , xgb.agree.cv[,4],pow=pow) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - XGB with lasso feature         ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , xgb.agree.cv[,5],pow=pow) 
    }
    if (sum(ensemble[c(3,4,7,8)])> 0) {
      cat ("\n lasso.minR - XGB with lasso offset          ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , xgb.agree.cv[,6],pow=pow)   
    }
  }

  if ((dolasso == 1) & (dorf == 1)) {
    cat ("\n lasso.minR - Random Forest                  ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , rf.agree.cv[,1],pow=pow) 
    
    if (sum(ensemble[c(2,6)])> 0) {
      cat ("\n lasso.minR - RF with lasso feature          ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , rf.agree.cv[,2],pow=pow) 
    }
    if ( (sum(ensemble[c(3,4,7,8)])> 0) & (family == "gaussian") ) {
      cat ("\n lasso.minR - RF with lasso offset           ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , xgb.agree.cv[,3],pow=pow)   
    }
  }

  if ((dolasso == 1) & (doann == 1)) {
    cat ("\n lasso.minR - ANN                            ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , ann.agree.cv[,1],pow=pow) 
    if (ensemble[6]) { 
      cat ("\n lasso.minR - ANN l lasso feature            ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , ann.agree.cv[,6],pow=pow)   
    } else if (ensemble[2]) { 
      cat ("\n lasso.minR - ANN lasso feature              ") ;  glmnetr.compcv0(lasso.agree.cv[,4] ,  ann.agree.cv[,2],pow=pow)   
    }  
    if (ensemble[8]) { 
      cat ("\n lasso.minR - ANN l lasso offset (upated)    ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , ann.agree.cv[,8],pow=pow)   
    } else if (ensemble[4]) { 
      cat ("\n lasso.minR - ANN lasso offset               ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , ann.agree.cv[,4],pow=pow)  
    } else if (ensemble[7]) { 
      cat ("\n lasso.minR - ANN l lasso offset (upated)    ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , ann.agree.cv[,7],pow=pow)   
    }  else if (ensemble[3]) { 
      cat ("\n lasso.minR - ANN lasso offset               ") ;  glmnetr.compcv0(lasso.agree.cv[,4] , ann.agree.cv[,3],pow=pow)   
    }
  }

  if (dolasso) { cat("\n") } 

  if ((doxgb == 1) & (dorf == 1)) {
    cat ("\n XGBoost (tuned) - RF                        ") ;  glmnetr.compcv0(xgb.agree.cv[,4] ,  rf.agree.cv[,1],pow=pow)   
    if (sum(ensemble[c(2,6)]) > 0) {
      cat ("\n XGBoost lasso feature-RF with lasso feature ") ;  glmnetr.compcv0(xgb.agree.cv[,5] ,  rf.agree.cv[,2],pow=pow)   
    } 
    if ( (sum(ensemble[c(3,4,7,8)]) > 0)  & (family == "gaussian") ) {
      cat ("\n XGBoost lasso offset-  RF with lasso offset ") ;  glmnetr.compcv0(xgb.agree.cv[,6] ,  rf.agree.cv[,3],pow=pow)   
    }
  }

  if ((doxgb == 1) & (doann == 1)) {
    cat ("\n XGBoost (tuned) - ANN                         ") ;  glmnetr.compcv0(xgb.agree.cv[,4] ,  ann.agree.cv[,1],pow=pow)   
    if (ensemble[6]) { 
      cat ("\n XGBoost lasso feature - ANN, l lasso feature ") ;  glmnetr.compcv0(xgb.agree.cv[,5] ,  ann.agree.cv[,6],pow=pow)   
    } else if (ensemble[2]) { 
      cat ("\n XGBoost lasso feature - ANN lasso feature   ") ;  glmnetr.compcv0(xgb.agree.cv[,5] ,  ann.agree.cv[,2],pow=pow)   
    } 
    if (family == "gaussian") {
      if (ensemble[8]) { 
        cat ("\n XGBoost lasso offset-ANN l lasso offset(upated) ") ;  glmnetr.compcv0(xgb.agree.cv[,6] ,  ann.agree.cv[,8],pow=pow)   
      } else if (ensemble[4]) { 
        cat ("\n XGBoost lasso offset - ANN, lasso offset   ") ;  glmnetr.compcv0(xgb.agree.cv[,6] ,  ann.agree.cv[,4],pow=pow)   
      } else if (ensemble[7]) { 
        cat ("\n XGBoost lasso offset - ANN l lasso offset  ") ;  glmnetr.compcv0(xgb.agree.cv[,6] ,  ann.agree.cv[,7],pow=pow)   
      } else if (ensemble[3]) { 
        cat ("\n XGBoost offset - ANN lasso offset          ") ;  glmnetr.compcv0(xgb.agree.cv[,6] ,  ann.agree.cv[,3],pow=pow)   
      }  
    }
  }
  
  if (doxgb) { cat("\n") }
  
  if ((dorf == 1) & (doann == 1)) {
    cat ("\n RF - ANN                                    ") ;  glmnetr.compcv0(rf.agree.cv[,1] ,  ann.agree.cv[,1],pow=pow) 
    if (ensemble[6]) {
      cat ("\n RF lasso feature - ANN  l lasso feature     " ) ;  glmnetr.compcv0(rf.agree.cv[,2] ,  ann.agree.cv[,6],pow=pow) 
    } else if (ensemble[2]) {
      cat ("\n RF lasso feature - ANN lasso feature        " ) ;  glmnetr.compcv0(rf.agree.cv[,2] ,  ann.agree.cv[,2],pow=pow)  
    }
    if (ensemble[8]) {
      cat ("\n RF lasso offset - ANN l lasso offset (upated) " ) ;  glmnetr.compcv0(rf.agree.cv[,3] ,  ann.agree.cv[,8],pow=pow)   
    } else if (ensemble[4]) {
      cat ("\n RF lasso offset - ANN lasso offset           " ) ;  glmnetr.compcv0(rf.agree.cv[,3] ,  ann.agree.cv[,4],pow=pow)  
    } else if (ensemble[7]) {
      cat ("\n RF lasso offset - ANN, l lasso offset (upated) " ) ;  glmnetr.compcv0(rf.agree.cv[,3] ,  ann.agree.cv[,7],pow=pow) 
    } else if (ensemble[3]) {
      cat ("\n RF lasso offset - ANN, lasso offset          " ) ;  glmnetr.compcv0(rf.agree.cv[,3] ,  ann.agree.cv[,3],pow=pow)  
    }
    cat("\n")
  }
  
  cat("\n")
}

####################################################################################################################################
####################################################################################################################################

