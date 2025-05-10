################################################################################  
#' Summarize a nested.glmnetr() output objects version 0.5-1
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
#' @param tuning 1 to print tuning parameters, 0 (default) to not print
#' @param width character width of the text body preceding the performance 
#' measures which can be adjusted between 60 and 120.
#' @param ... Additional arguments passed to the summary function.  
#' 
#' @return - a nested cross validation fit summary, or a cross validation model summary.  
#' 
#' @seealso
#'   \code{\link{nested.compare}} , \code{\link{nested.cis}} , \code{\link{summary.cv.glmnetr}} , \code{\link{roundperf}} , 
#'   \code{\link{plot.nested.glmnetr}} , \code{\link{calplot}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
#' @noRd
#'
# cvfit = FALSE ; pow=2 ; printg1 = FALSE ; digits = 3 ; Call=NULL ; onese = 0 ; table = 1 ; 

summary.nested.glmnetr_0_5_1 = function(object, cvfit=FALSE, pow=2, printg1=FALSE, 
                                  digits=4, Call=NULL, onese=0, table=1, tuning=0, width=84, ...) {
# cvfit=FALSE ; pow=2 ; printg1=FALSE ; digits=4 ; Call=NULL ; onese=0 ; table=1 ; tuning=0 ; width=108 ;  
  
  get.DevRat = function( devian.rep, null.m2LogLik.rep, sat.m2LogLik.rep, n.rep, bootstrap=0 ) {
    if (bootstrap == 0) {
      AllDevRat = colSums( n.rep*( null.m2LogLik.rep - devian.rep ))/sum((null.m2LogLik.rep - sat.m2LogLik.rep)*n.rep )
    } else {
      AllDevRat = colMeans( ( null.m2LogLik.rep - devian.rep ) / (null.m2LogLik.rep - sat.m2LogLik.rep) ) 
    }
    return( AllDevRat )
  }
  
#  get.DevRat = function( devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv ) {
#    ncoldevian = dim(devian.cv)[2] 
#    AllDevRat = rep(1,ncoldevian)
#    for (j_ in c(1:ncoldevian)) {
#      AllDevRat[j_] = devrat_(devian.cv[,j_], null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )[[2]]
#    }
#    return( AllDevRat )
#  }
  
  mlength = 116 ## biggest for asthetics 
  mlength = 108 ## 3 lines  
  mlength = 95 ## 4 lines for all models 
  mlength = 77 ## 5 lines 
  mlength = 80 ## "standard" width
  mlength = 84 ## good for R markdown 
  if (is.null(width)) { width = 80 }
  mlength = width 
  
  if (cvfit==TRUE) {
    cv_glmnet_fit = object$cv_glmnet_fit  
    summary(cv_glmnet_fit,printg1=printg1)
  } else {
    if (!(table %in% c(0,1,2,3))) { table = 1 }
    if (!is.null(Call)) { 
      if (Call != 0) { Call = object$call 
      } else { Call = NULL }
    }
    sample  = object$sample 
#    tuning_  = object$tuning
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
    do_ncv  = object$do_ncv
    
    family = sample[1]
#    sample[6] = round(as.numeric(sample[6]), digits=digits)
#    if (family=="cox") { names(sample)[6]="null.dev/events" 
#    } else { names(sample)[6] = "null.dev/obs"  }
    
    if ( pow == 2 ) { gagree = "R-square" } else { gagree = "Correlation" }
    if (family %in% c("cox","binomial")) { gagree = "Concordance" ; pow = 1 } 

    annrownames = c("ANN Uninformed          ", "ANN lasso features             ", "ANN lasso weights             ", "ANN lasso weights, reset      ",
                    "ANN only lasso terms    ", "ANN lasso terms, lasso features", "ANN lasso terms, lasso weights", "ANN lasso terms, weights reset") 

    annrownames = c("ANN Uninformed          ", "ANN lasso features         ", "ANN lasso weights         ", "ANN lasso weights, reset  ",
                    "ANN only lasso terms    ", "ANN lasso terms, lasso feat", "ANN lasso terms, lasso wts", "ANN lasso terms, wts reset") 
    
    colnames1 = c("Ave DevRat", "Ave Int", "Ave Slope", paste("Ave", gagree), "Ave Non Zero", "Naive DevRat", paste("Naive", gagree), "Non Zero" )
    colnames0 = c( "Naive DevRat", paste("Naive", gagree), "Non Zero" )

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

    if (!is.null(Call)) { 
      cat(paste0("\n     function call :\n\n"))  
      print(Call) 
      cat(paste0("\n"))  
    }
    
    if ( table %in% c(1,2,3) ) { 
      cat(paste0("  Sample information including number of records, "))     
      if (family %in% c("cox","binomial")) { cat(paste0("events, ")) }
      cat(paste0( "number of columns in", "\n  design (predictor, X) matrix, and df (rank) of design matrix: ", "\n") )
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
            
      if (family == "cox") { 
        print(sample,quote=FALSE)
      } else if (family == "binomial") { 
        print(sample[-c(7,8)],quote=FALSE)
      } else { 
        print(sample[-c(3,7,8)],quote=FALSE)
      }
      
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
      
      if (tuning == 1) {
      
        ensemble = object$ensemble
        if (sum(ensemble[c(1,5)]) > 0) { base = 1 }
        if (sum(ensemble[c(2,6)]) > 0) { feat = 1 }
        if (sum(ensemble[c(3,4,7,8)]) > 0) { offs = 1 }
        
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
        
        if ( is.null(object$ver) ) { xx = 0 
        } else if ( substr(object$ver[2],1,5) == "0.4-3" ) { xx = 0
        } else { xx = 1 }
        
        if ( xx == 1 ) {
          
          if (doxgb == 1) {
            if (base == 1) { x=object$doxgb_tuned 
            } else if (feat == 1) { x=object$doxgb_tunedF 
            } else if (offs == 1) { x=object$doxgb_tunedO }
            cat(paste0("\n"  , " Tuning parameters for Gradient Boosting Machine models\n"))
            print( c( nfolds=x$nfold[1], nrounds=x$nrounds, early_stopping_rounds=x$early_stopping_rounds, 
               eta=x$eta, gamma=x$gamma, max_depth=x$max_depth, min_child_weight=x$min_child_weight, 
               colsample_bytree=x$colsample_bytree, lambda=x$lambda, alpha=x$alpha, 
               subsample=x$subsample, nrounds.final=x$nrounds.final) )
          }
          
          if (dorf == 1) {
            x = NULL
            if (base == 1) { x=object$dorf_base 
            } else if (feat == 1) { x=object$dorf_feat 
            } else if ((offs == 1) & (family == "gaussian")) { x=object$dorf_offs }
            if (!is.null(x)) { 
              cat(paste0("\n"  , " Tuning parameters for Random Forest models\n"))
              print( c(mtryc=x$mtryc, ntreec=x$ntreec, nsplitc=x$nsplitc) )
            } 
          }
          
          if (doorf == 2) {
            x = NULL
            if (base == 1) { x=object$doorf_base 
            } else if (feat == 1) { x=object$doorf_feat 
            } else if ((offs == 1) & (family == "gaussian")) { x=object$doorf_offs }
            if (!is.null(x)) { 
              cat(paste0("\n"  , " Tuning parameters for Oblique Random Forest models\n"))
              print( c(mtryc=x$mtryc, ntreec=x$ntreec, nsplitc=x$nsplitc) )
            } 
          }
          
        }
        
        if (doann  == 1) {
          cat(paste0("\n"  , " Tuning parameters for ", trimws(whichann), " model : ", "\n") )
          #        print(ann_cv$modelsum[c(1:9,11:12)]) 
          class( ann_cv$modelsum[c(1:9,11:12)] )
          x = ann_cv$modelsum[c(1:9,11:12)]
          namesx = names(x)
          mx = matrix(x,nrow=1,ncol=11) 
          colnames(mx) = namesx
          print(data.frame(mx)[1,c(1:11)], row.names=FALSE)
        }
        
      }
      
    }
    
    if (family == "cox") { perunit = "deviance per event " } else { perunit = "deviance per record " }
    
    null.m2LogLik.cv = object$null.m2LogLik.cv
    sat.m2LogLik.cv = object$sat.m2LogLik.cv
    n.cv = object$n.cv 
    
    last = "none"
    if (doaic  == 1) { last = "aic"  
    } else if (dostep == 1) { last = "step" 
    } else if (dorpart== 1) { last = "rpart" 
    } else if (doann  == 1) { last = "ann"   
    } else if (doorf  == 1) { last = "orf"    
    } else if (dorf   == 1) { last = "rf"    
    } else if (doxgb  == 1) { last = "xgb"   
    } else if (dolasso== 1) { last = "lasso" }

    if (!(table %in% c(0,3))) {
      hstring = "For" 
      i_ =  2
      if (dolasso== 1) { 
        hstring[2] = c("LASSO,") 
        i_ = i_ + 1 
      }
      if (doxgb  == 1) { 
        hstring[3] = "gradient boosting machine using XGBoost (XGB),"
        i_ = i_ + 1 
      }
      if (dorf   == 1) { 
        hstring[c(i_:(i_+2))] = c("Random", "Forest", "(RF),")
        i_ = i_ + 3 
      }
      if (doorf   == 1) { 
        hstring[c(i_:(i_+3))] = c("Oblique", "Random", "Forest", "(ORF),")
        i_ = i_ + 4 
      }
      if (doann  == 1) { 
        hstring[c(i_:(i_+3))] = c("Artificial", "Neural", "Networks", "(ANN),") 
        i_ = i_ + 4 
      }
      if (dorpart== 1) { 
        hstring[c(i_:(i_+4))] = c("Recursive", "And", "Partitioning", "Trees", "(RPART),") 
        i_ = i_ + 5 
      }
      if ((dostep== 1) & (doaic==1)) {
        hstring[c(i_:(i_+8))] = c( "Stepwise", "regression", "tuned", "by", "df", 
          "and", "p,", "and", "AIC,") 
        i_ = i_ + 9
      } else if (dostep==1) {
        hstring[c(i_:(i_+7))] = c( "and", "Stepwise", "regression", "tuned", "by", "df", 
                                   "and", "p,") 
        i_ = i_ + 8
      } else if (doaic ==1) {
        hstring[c(i_:(i_+4))] = c("and", "Akaike", "Information", "Criterion", "(AIC),") 
        i_ = i_ + 5 
      }
      if (do_ncv == 1) {
        hstring0 = c("average", "(Ave)", "model", "performance", "measures", "from",
                     "the", paste0(object$tuning[1], "-fold"), "(nested)", "cross", "validation", "are", "given", "together", "with", "naive",
                     "summaries", "calculated", "using", "all", "data", "without", "cross", "validation")
        tlen = length( hstring0 )
        hstring[c(i_:(i_+tlen-1))] = hstring0 
        i_ = i_ + tlen -1
      } else { 
        hstring0 = c("naive", "performance", "measures", "calculated", "using", 
            "all", "data", "without", "cross", "validation", "are", "given")
        tlen = length( hstring0 )
        hstring[c(i_:(i_+tlen-1))] = hstring0 
        i_ = i_ + tlen -1
        }
      lhstring = length(hstring)
#      print (lhstring) 
#      print(i_)
#      print (hstring) 
      
      mlength = min(mlength, 120)
      if (do_ncv == 0) { mlength = min(76,mlength) } 
      mlength = max(mlength,60) 
      
      cat("\n")
      clen = 0 
      for (i_ in c(1:lhstring)) { 
#        print(i_)
        if ( (clen + nchar(hstring[i_]) + 1) > mlength ) {
          cat(("\n")) 
          clen = 0 
        }
        if (hstring[i_] == "-fold") { cat(hstring[i_])  
          } else { cat(paste0(" ", hstring[i_])) }
        clen = clen + nchar(hstring[i_]) + 1 
      }
      cat("\n\n")
    }

    ## CALCULATE PERFORMANCE SUMMARIES #########################################
    ## CALCULATE PERFORMANCE SUMMARIES #########################################

    lasso = NULL 
    if (dolasso == 1) {
      if (do_ncv ==1) {
        lassoAllDevRat = get.DevRat( object$lasso.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
        lassoAveDevian = colMeans(object$lasso.devian.cv, na.rm=TRUE)  
        lassoAveIntcal = colMeans(object$lasso.intcal.cv, na.rm=TRUE)  
        lassoAveLincal = colMeans(object$lasso.lincal.cv, na.rm=TRUE)  
        lassoAveAgree  = colMeans(object$lasso.agree.cv, na.rm=TRUE)
        lassoAveNzero  = colMeans(object$lasso.nzero.cv, na.rm=TRUE)  
      }
      lasso.agree.naive = object$lasso.agree.naive 
      if (family == "gaussian") { 
        lasso.agree.naive = lasso.agree.naive ^pow
        if (do_ncv ==1) { lassoAveAgree = lassoAveAgree ^pow } 
      }
      ## sqrt( apply(lasso.agree.cv,2,var) ) 

      #        ((m2.ll.null - m2.ll.mod)/(m2.ll.null - m2.ll.sat ))
      lasso.devrat.naive = (as.numeric(object$sample[7]) - object$lasso.devian.naive) / 
        (as.numeric(object$sample[7]) - as.numeric(object$sample[8]))
  
      if (do_ncv == 1) {
        lasso = data.frame( lassoAllDevRat , lassoAveIntcal, lassoAveLincal , lassoAveAgree, lassoAveNzero, 
                            lasso.devrat.naive, object$lasso.agree.naive, object$lasso.nzero)
        names(lasso) = colnames1 
      } else {
        lasso = data.frame( lasso.devrat.naive, object$lasso.agree.naive, object$lasso.nzero) 
        names(lasso) = colnames0 
      }
      rownames = paste0(c(rep("LASSO ",6),""), row.names(lasso)) 
      rownames[7] = "Ridge                      " 
      row.names(lasso) = rownames 
      
      if (onese == 0) { lasso = lasso[c(2,4,6,7),] } 
      lassor = roundperf(lasso, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { lassor = lassor[,-2] }      
    }        
      
    ## XGB #####################################################################
    
    xgb = NULL 
    if (doxgb == 1) { 
      en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
      en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
      en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
      enx = c(en1, en2, en3, en1, en2, en3) 
      if (do_ncv ==1) {
        xgbAllDevRat = get.DevRat( object$xgb.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
        xgbAveDevian = colMeans( object$xgb.devian.cv, na.rm=TRUE )
        xgbAveIntcal = colMeans( object$xgb.intcal.cv, na.rm=TRUE )
        xgbAveLincal = colMeans( object$xgb.lincal.cv, na.rm=TRUE )
        xgbAveAgree  = colMeans( object$xgb.agree.cv, na.rm=TRUE ) 
        xgbAveNzero  = colMeans( object$xgb.nzero.cv, na.rm=TRUE ) 
      }
      xgb.agree.naive = object$xgb.agree.naive
      if (family == "gaussian") { 
        xgb.agree.naive = xgb.agree.naive ^pow 
        if (do_ncv ==1) { xgbAveAgree = xgbAveAgree ^pow }
      }
      xgb.devrat.naive = (as.numeric(object$sample[7]) - object$xgb.devian.naive) / 
                   (as.numeric(object$sample[7]) - as.numeric(object$sample[8]))
      if (do_ncv == 1) {
        xgb = data.frame( xgbAllDevRat , xgbAveIntcal, xgbAveLincal , xgbAveAgree, xgbAveNzero, xgb.devrat.naive, xgb.agree.naive, object$xgb.nzero )
        names( xgb ) = colnames1 
#        xgb
      } else {
        xgb = data.frame( xgb.devrat.naive, xgb.agree.naive , object$xgb.nzero )
        names(xgb) =colnames0
      }
      row.names(xgb) = c("XGB (not tuned)            ", "XGB lasso Feature          ", "XGB lasso Offset           ", 
                         "XGB Tuned                  ", "XGB Tuned lasso Feature    ", "XGB Tuned lasso Offset     " )
      xgb = xgb[enx==1,]
      xgbr = roundperf(xgb, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { xgbr = xgbr[,-2] }      
    }   
    
    ##### Random Forest ########################################################
    
    rf = NULL 
    if (dorf == 1) { 
      en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
      en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
      en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
      if (family != "gaussian") { en3 = 0 }
      enx = c(en1, en2, en3) 
      if (do_ncv ==1) {
        rfAllDevRat = get.DevRat( object$rf.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
        rfAveDevian = colMeans(object$rf.devian.cv, na.rm=TRUE) # [enx==1]   
        rfAveIntcal = colMeans(object$rf.intcal.cv, na.rm=TRUE)
        rfAveLincal = colMeans(object$rf.lincal.cv, na.rm=TRUE)
        rfAveAgree  = colMeans(object$rf.agree.cv, na.rm=TRUE)
        rfAveMtry   = colMeans(object$rf.mtry.cv, na.rm=TRUE)
      }
      rf.agree.naive = object$rf.agree.naive
      if (family == "gaussian") { 
        rf.agree.naive = rf.agree.naive ^pow 
        if (do_ncv ==1) { rfAveAgree  = rfAveAgree ^pow }
      }
      rf.devrat.naive = (as.numeric(object$sample[7]) - object$rf.devian.naive) / 
        (as.numeric(object$sample[7]) - as.numeric(object$sample[8]))
      if (do_ncv == 1) {
        rf = data.frame( rfAllDevRat , rfAveIntcal , rfAveLincal , rfAveAgree, rfAveMtry, rf.devrat.naive, rf.agree.naive, object$rf.mtry  )
        names( rf ) = colnames1 
        rf
      } else {
        rf = data.frame( rf.devrat.naive, rf.agree.naive , object$rf.mtry )
        names(rf) = colnames0 
      }
      row.names(rf) = c("RF Simple                  ", "RF lasso Feature           ", "RF lasso Offset            " )
      rf = rf[enx==1,]
      rfr = roundperf(rf, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { rfr = rfr[,-2] }      
    }

    ##### Oblique Random Forest ########################################################
    
    orf = NULL 
    if (doorf == 1) { 
      en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
      en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
      en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
      if (family != "gaussian") { en3 = 0 }
      enx = c(en1, en2, en3) 
      if (do_ncv ==1) {
        orfAllDevRat = get.DevRat( object$orf.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
        orfAveDevian = colMeans(object$orf.devian.cv, na.rm=TRUE) # [enx==1]   
        orfAveIntcal = colMeans(object$orf.intcal.cv, na.rm=TRUE)
        orfAveLincal = colMeans(object$orf.lincal.cv, na.rm=TRUE)
        orfAveAgree  = colMeans(object$orf.agree.cv, na.rm=TRUE)
        orfAveMtry   = colMeans(object$orf.mtry.cv, na.rm=TRUE)
      }
      orf.agree.naive = object$orf.agree.naive
      if (family == "gaussian") { 
        orf.agree.naive = orf.agree.naive ^pow 
        if (do_ncv ==1) { orfAveAgree  = orfAveAgree ^pow }
      }
      orf.devrat.naive = (as.numeric(object$sample[7]) - object$orf.devian.naive) / 
        (as.numeric(object$sample[7]) - as.numeric(object$sample[8]))
      if (do_ncv == 1) {
        orf = data.frame( orfAllDevRat , orfAveIntcal , orfAveLincal , orfAveAgree, orfAveMtry, orf.devrat.naive, orf.agree.naive, object$orf.mtry  )
        names( orf ) = colnames1 
        orf
      } else {
        orf = data.frame( orf.devrat.naive, orf.agree.naive , object$orf.mtry )
        names(orf) = colnames0 
      }
      row.names(orf) = c("ORF Simple                  ", "ORF lasso Feature           ", "ORF lasso Offset            " )
      orf = orf[enx==1,]
      orfr = roundperf(orf, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { orfr = orfr[,-2] }      
    }
    
    ##### ANN ##################################################################
    
    ann = NULL 
    if (doann == 1) { 
      if (do_ncv ==1) {
        annAllDevRat = get.DevRat( object$ann.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
        annAveDevian = colMeans(object$ann.devian.cv, na.rm=TRUE)               # [(ensemble==1)]
        annAveIntcal = colMeans(object$ann.intcal.cv, na.rm=TRUE)
        annAveLincal = colMeans(object$ann.lincal.cv, na.rm=TRUE)
        annAveAgree  = colMeans(object$ann.agree.cv , na.rm=TRUE)
        annAveNzero  = colMeans(object$ann.nzero.cv , na.rm=TRUE)
      }
      ann.agree.naive = object$ann.agree.naive
      if (family == "gaussian") { 
        ann.agree.naive = ann.agree.naive ^pow
        if (do_ncv ==1) { annAveAgree = annAveAgree ^pow }
      }
      ann.devrat.naive = (as.numeric(object$sample[7]) - object$ann.devian.naive) / 
                      (as.numeric(object$sample[7]) - as.numeric(object$sample[8]))
      if (do_ncv == 1) {
        ann = data.frame( annAllDevRat , annAveIntcal , annAveLincal , annAveAgree, annAveNzero, ann.devrat.naive, ann.agree.naive, object$ann.nzero )
        names( ann ) = colnames1 
      } else {
        ann = data.frame( ann.devrat.naive, ann.agree.naive , object$ann.nzero )
        names( ann ) = colnames0 
      }
      row.names(ann)  = annrownames 
      ann = ann[ensemble[c(1:8)]==1,]
      annr = roundperf(ann, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { annr = annr[,-2] }      
    }                                                                           
    
    ##### RPART ################################################################
    
    rpart = NULL 
    if (dorpart == 1) { 
      en1 = ifelse (sum(ensemble[c(1,5)])>=1, 1, 0)
      en2 = ifelse (sum(ensemble[c(2,6)])>=1, 1, 0)
      en3 = ifelse (sum(ensemble[c(3,4,7,8)])>=1, 1, 0)
      if (family == "binomial") { en3 = 0 }
      enx = c(en1,en1,en1, en2,en2,en2, en3,en3,en3) 
      if (do_ncv ==1) {
        rpartAllDevRat = get.DevRat( object$rpart.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
        rpartAveDevian = colMeans(object$rpart.devian.cv, na.rm=TRUE) 
        rpartAveIntcal = colMeans(object$rpart.intcal.cv, na.rm=TRUE) 
        rpartAveLincal = colMeans(object$rpart.lincal.cv, na.rm=TRUE) 
        rpartAveAgree  = colMeans(object$rpart.agree.cv, na.rm=TRUE) 
        rpartAveNzero  = colMeans(object$rpart.nzero.cv, na.rm=TRUE)                          # [c(3,2,1,6,5,4,9,8,7)]          # [enx==1]
      }
      rpart.agree.naive = object$rpart.agree.naive
      if (family == "gaussian") { 
        rpart.agree.naive = rpart.agree.naive ^pow 
        if (do_ncv ==1) { rpartAveAgree = rpartAveAgree ^pow }
      }
      rpart.devrat.naive = (as.numeric(object$sample[7]) - object$rpart.devian.naive) / 
        (as.numeric(object$sample[7]) - as.numeric(object$sample[8]))
      if (do_ncv == 1) {
        rpart = data.frame( rpartAllDevRat , rpartAveIntcal , rpartAveLincal , rpartAveAgree, rpartAveNzero, rpart.devrat.naive, rpart.agree.naive, object$rpart.nzero )
        #        rpart = rpart[c(3,2,1,6,5,4,9,8,7),]
        names( rpart ) = colnames1 
        rpart
      } else {
        rpart = data.frame( rpart.devrat.naive, rpart.agree.naive , object$rpart.nzero )
        names( rpart ) = colnames0  
      }
      rownames = row.names(rpart)
      rownames[c(4:9)] = substr(rownames[c(4:9)], 3,10)
      #        rownames = paste(c("RPART", "RPART", "RPART", "RPART lasso feature", "RPART lasso feature", "RPART lasso feature",
      #                           "RPART lasso offset", "RPART lasso offset", "RPART lasso offset"), rownames)
      rownames = paste(c("RPART", "RPART", "RPART", "RPART lasso feat", "RPART lasso feat", "RPART lasso feat",
                         "RPART lasso offs", "RPART lasso offs", "RPART lasso offs"), rownames,"       ")
      row.names(rpart) = rownames 
      rpart = rpart[enx==1,]
      rpartr = roundperf(rpart, digits, do_ncv) 
      if ((family == "cox") & (do_ncv==1)) { rpartr = rpartr[,-2] }      
    }      
    
    ##### STEPWISE p or df tuned & AIC #########################################

    step = NULL 
    if ((dostep==1) | (doaic==1)) {
      if (do_ncv == 1) {
      StepAllDevRat = get.DevRat( object$step.devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )
      object$step.devian.cv[ is.na(object$step.devian.cv) ] = max(object$step.devian.cv)
      StepAveDevian = colMeans( object$step.devian.cv, na.rm=TRUE)
      StepAveIntcal = colMeans( object$step.intcal.cv, na.rm=TRUE)
      StepAveLincal = colMeans( object$step.lincal.cv, na.rm=TRUE)
      StepAveAgree  = colMeans( object$step.agree.cv , na.rm=TRUE)
      StepAve.nzero = colMeans( object$step.nzero.cv , na.rm=TRUE)
      StepAve.p     = colMeans( object$step.p.cv     , na.rm=TRUE)
      }
      step.agree.naive = object$step.agree.naive 
      if (family == "gaussian") { 
        if (do_ncv == 1) { StepAveAgree = StepAveAgree ^pow }
      }
      step.devrat.naive = (as.numeric(object$sample[7]) - object$step.devian.naive) / 
        (as.numeric(object$sample[7]) - as.numeric(object$sample[8]))
      if (do_ncv == 1) {
        step = data.frame( StepAllDevRat , StepAveIntcal , StepAveLincal , StepAveAgree, StepAve.nzero, step.devrat.naive, step.agree.naive, object$step.nzero )
        names( step ) = colnames1 
#        names( step ) = c("Deviance", "Inter", "Slope", "Concor", "Non Zero", "Deviancd", "ConroR", "Non Zero" )
      } else {
        step = data.frame( step.devrat.naive, step.agree.naive , object$step.nzero )
        names(step) = colnames0 
      }
      row.names(step) = c("Stepwise df tuned          ", "Stepwise p tuned       ", "Stepwise AIC             " )
      if        ((dostep == 1) & (doaic == 0)) { 
        step = step[c(1,2),] 
      } else if ((dostep == 0) & (doaic == 1)) { 
        step = step[c(3),]  
      }
      stepr = roundperf(step, digits, do_ncv)
      if ((family == "cox") & (do_ncv==1)) { stepr = stepr[,-2] }      
    }
    
    ## PRINT ACTUAL PERFORMANCE SUMMARIES ######################################
    ## PRINT ACTUAL PERFORMANCE SUMMARIES ######################################
    
    allstat = rbind(lasso, xgb, rf, orf, ann, rpart, step)
    mlen = max( nchar( trimws(rownames(allstat)) ) ) 
    rnms = rownames(allstat)
    rnms = substring ( rnms, 1, mlen)
    rownames(allstat) = rnms
    
    if (dolasso == 1) {
      if ( table %in% c(3) ) {
        if (do_ncv == 1) {
          cat( "\n  LASSO: Ave is for (nested) CV model performance summary, else \n",
               "        naive summary for fit on all data \n" ) 
        } else {
          cat( "\n  LASSO: Naive summary for fit on all data \n" )         
        }
      }
      if ( table %in% c(1,2,3) ) { 
        rownames(lassor) = substring(rownames(lassor),1,mlen)
        print( lassor ) 
        cat("\n") 
      }
    }
    
    if (doxgb == 1) {
      if ( table %in% c(3) ) {
        if (do_ncv == 1) {
          cat( "\n  XGBoost: Ave is for (nested) CV model performance summary, else\n",
               "          naive summary for fit on all data \n" ) 
        } else {
          cat( "\n  XGBoost: Naive model performance summary for fit on all data\n" ) 
        }
      }
      if ( table %in% c(1,2,3) ) { 
        rownames(xgbr) = substring(rownames(xgbr),1,mlen)
        print( xgbr ) 
        cat("\n") 
      }
    }
    
    if (dorf == 1) { 
      if ( table %in% c(3) ) {
        if (do_ncv == 1) {
          cat( "\n  Random Forest: Ave is for CV model performance summary, else\n",
               "                   naive summary for fit on all data \n" ) 
        } else {
          cat( "\n  Random Forest: Naive model performance summary for fit on all data\n" ) 
        }
      }
      if ( table %in% c(1,2,3) ) { 
        rownames(rfr) = substring(rownames(rfr),1,mlen)
        print( rfr )
        cat("\n") 
      }
    }
    
    if (doorf == 1) { 
      if ( table %in% c(3) ) {
        if (do_ncv == 1) {
          cat( "\n  Random Forest: Ave is for CV model performance summary, else\n",
               "                   naive summary for fit on all data \n" ) 
        } else {
          cat( "\n  Random Forest: Naive model performance summary for fit on all data\n" ) 
        }
      }
      if ( table %in% c(1,2,3) ) { 
        rownames(orfr) = substring(rownames(orfr),1,mlen)
        print( orfr )
        cat("\n") 
      }
    }
      
    if (doann == 1) { 
      if ( table %in% c(3) ) {
        if (do_ncv == 1) {
          cat( "\n  Artificial Neural Network: Ave is for (nested) CV model performance summary, else\n",
               "                         naive summary for fit on all data \n" ) 
        } else {
          cat( "\n  Artificial Neural Network: Naive model performance summary for fit on all data\n" ) 
        }
      }
      if ( table %in% c(1,2,3) ) { 
        rownames(annr) = substring(rownames(annr),1,mlen)
        print( annr) 
        cat("\n") 
      }
    }
    
    if (dorpart == 1) {
      if (table == 3) {
        if (do_ncv == 1) {
          cat( "\n  Recursive Partitioning: Ave is for CV model performance summary, else\n", 
               "                         naive summary for fit on all data \n" ) 
        } else {
          cat( "\n  Recursive Partitioning: Naive model performance summary for fit on all data\n" ) 
        }
      }
#      if ((en2 == 1) | (en3 == 1)) { cat( "          l. feat for lasso included as feature, l. offs for lasso included as offset\n" )
#      } else if (en2 == 1) { cat( "                  l. feat for lasso included as feature\n" )
#      } else if (en3 == 1) { cat( "                  l. offs for lasso included as offset\n" ) }
      if ( table %in% c(1,2,3) ) { 
        rownames(rpartr) = substring(rownames(rpartr),1,mlen)
        print( rpartr ) 
        cat("\n") 
      } 
    } 
    
    if ((dostep == 1) | (doaic == 1)) {
      if ( table %in% c(3) ) {
        if (do_ncv == 1) {
          if        ((dostep == 1) & (doaic == 1)) { 
            cat( "\n  Stepwise tuned and AIC: Ave is for (nested) CV model performance summary, else\n")
          } else if ((dostep == 1) & (doaic == 0)) { 
            cat( "\n  Stepwise tuned: Ave is for (nested) CV model performance summary, else\n")
          } else if ((dostep == 0) & (doaic == 1)) { 
            cat( "\n  Stepwise AIC: Ave is for (nested) CV model performance summary, else\n")
          }
          cat( "                          naive summary for fit on all data \n" ) 
        } else {
          cat( "\n  Stepwise tuned and AIC: Naive model performance summary for fit on all data\n" )
        } 
      }
      if ( table %in% c(1,2,3) ) {
        rownames(stepr) = substring(rownames(stepr),1,mlen)
        print( stepr ) 
        cat("\n") 
      }
    }
    
#    cat("\n  Correlation of predictors \n")  
#    xbetas = object$xbetas  
#    xbetas0 = xbetas[,(diag(cov(xbetas)) > 0)] 
#    pearson = cor(xbetas0, method="pearson") 
#    spearman = cor(xbetas0, method="spearman") 
#    corbeta = upper.tri(pearson,diag=TRUE)*pearson + lower.tri(spearman)*spearman
#    round( corbeta, digits=4)
#    if (table != 0) { cat("\n") }
#    (names(lasso) == names(ann))*1 
    if ( table %in% c(0,2) ) { return( rbind(lasso, xgb, rf, orf, ann, rpart, step) ) }
  } #### end of if summary 
} 

####################################################################################################################################
####################################################################################################################################

