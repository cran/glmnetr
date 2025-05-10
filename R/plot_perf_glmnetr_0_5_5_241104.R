################################################################################
##### plot.perf.glmnetr_yymmdd.R ###############################################
################################################################################
#' Plot nested cross validation performance summaries
#' 
#' @description 
#' This function plots summary information from a nested.glmnetr() output object, that 
#' is from a nested cross validation performance.  Alternamvely one can output the 
#' numbers otherwise displayed to a list for extraction or customized plotting.  Performance
#' measures for plotting include "devrat" the deviance ratio, i.e. the fractional 
#' reduction in deviance relative to the null model deviance, "agree" a measure 
#' of agreement, "lincal" the slope from a linear calibration and "intcal" the 
#' intercept from a linear calibration.  Performance measure estimates 
#' from the individual (outer) cross validation fold are depicted by thin lines 
#' of different colors and styles, while the composite value from all folds is 
#' depicted by a thicker black line, and the performance measures naively 
#' calculated on the all data using the model derived from all data is 
#' depicted by a thicker red line. 
#' 
#' @param x A nested.glmnetr output object
#' @param type determines what type of nested cross validation performance measures are 
#' plotted.  Possible values are 
#' "devrat" to plot the deviance ratio, i.e. the fractional reduction in 
#' deviance relative to the null model deviance, 
#' "agree" to plot agreement in terms of concordance, correlation or R-square, 
#' "lincal" to plot the linear calibration slope coefficients, 
#' "intcal" to plot the linear calibration intercept coefficients, 
#' from the (nested) cross validation. 
#' @param pow Power to which agreement is to be raised when the "gaussian" model 
#' is fit, i.e. 2 for R-square, 1 for correlation.  Does not apply to type = "lasso". 
#' @param ylim y axis limits for model perforamnce plots, i.e. does not apply to 
#' type = "lasso".  The ridge model may calibrate very poorly obscuring plots for 
#' type of "lincal" or "intcal", so one may specify the ylim value.  If ylim is 
#' set to 1, then the program will derive a reasonable range for ylim.  If ylim is 
#' set to 0, then the entire range for all models will be displayed.  Does not 
#' apply to type = "lasso". 
#' @param fold By default 1 to display using a spaghetti the performance as 
#' calculated from the individual folds, 0 to display using dots only the composite 
#' values calculated using all folds.  
#' @param xgbsimple 1 (default) to include results for the untuned XGB model, 0 to not include.  
#' @param plot By default 1 to produce a plot, 0 to return the data used in the 
#' plot in the form of a list. 
#' 
#' @return This program returns a plot to the graphics window by default, and returns
#' a list with data used in teh plots if the plot=1 is specified.
#' 
#' @seealso
#'   \code{\link{plot.nested.glmnetr}} , \code{\link{nested.glmnetr}}
#'
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @importFrom stats interaction.plot 
#' 
#' @export
#' 
plot_perf_glmnetr_0_5_5 = function( x, type="devrat", pow=2, ylim=1, fold=1, xgbsimple=0, plot=1 ) {
  
  object = x 
  
  resample = object$resample 
  if (is.null(resample)) { resample = 1
  } else if (is.na(resample)) { resample = 1 }
  
  if (resample == 0) {
    warning("  Model performance plots are not generated when resample = 0")
  } else {
    if (is.null(type)) { type == "devrat" }
    if (is.null(fold)) { fold = 1 
    } else {
      if (fold > 0.6) { fold=round(fold)
      } else if (fold < 0) { fold = 1 
      } else { fold=round(fold,digits=1) }
    }

    family = object$sample[1] 
    
    if (is.null(pow)) { 
      if (family != "gaussian") { pow = 1
      } else {pow = 2} 
    } 
    
    if (type != "agree") { pow = 1} 
    
    if (type == "agree") { 
      if (family != "gaussian") { 
        pow = 1 
        ylab = "NCV Concordance"
      } else if (pow == 1) {
        ylab = "NCV Correlation"  
      } else if (pow==2) {
        ylab = "NCV R-square"  
      } else {
        ylab = paste0("Correlation^",pow)  
      }
    } else if (type == "intcal") { ylab = "NCV Calibration Intercept"
    } else if (type == "lincal") { ylab = "NCV Calibration Slope"
    } else if (type == "devian") { ylab = "NCV Deviance" 
    } else if (type == "devrat") { ylab = "NCV Deviance Ratio" 
    } else { ylab = "None" }
    
    ensemble = object$ensemble 
    ensemble
    
    if (is.null(object$null.m2LogLik.rep)) { object$null.m2LogLik.rep = object$null.m2LogLik.cv }
    if (is.null(object$sat.m2LogLik.rep )) { object$sat.m2LogLik.rep  = object$sat.m2LogLik.cv  }
    if (is.null(object$n.rep)) { object$n.rep = object$n.cv }
    m2.ll.null = object$null.m2LogLik.rep
    m2.ll.sat  = object$sat.m2LogLik.rep
    n__ = object$n.rep
    
    null.dev = as.numeric(object$sample[6])
    sat.m2LogLik = as.numeric(object$sample[8])
    null.m2LogLik = null.dev - sat.m2LogLik
    
    toplot1= c(NA,NA,NA)
    toplot2= c(NA,NA,NA)
    toplot3= c(NA,NA,NA,NA)
    nms = NA
    frst = 1 
    i_ = 0 
    
    ##### LASSO ################################################################
    
    if (object$fits[1] == 1) {
      
      if (is.null(object$lasso.devian.rep)) { object$lasso.devian.rep = object$lasso.devian.cv }
      if (is.null(object$lasso.intcal.rep)) { object$lasso.intcal.rep = object$lasso.intcal.cv }
      if (is.null(object$lasso.lincal.rep)) { object$lasso.lincal.rep = object$lasso.lincal.cv }
      if (is.null(object$lasso.agree.rep )) { object$lasso.agree.rep  = object$lasso.agree.cv  }
      
      if         (type == "agree") {  object1 = object$lasso.agree.rep 
      } else if (type == "lincal") {  object1 = object$lasso.lincal.rep 
      } else if (type == "intcal") {  object1 = object$lasso.intcal.rep 
      } else if (type == "devian") {  object1 = object$lasso.devian.rep 
      } else if (type == "devrat") {
        object1 = object$lasso.devian.rep
        devratl = list()
        for (j_ in c(1:7)) {
          devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
        }
      }
      nfold = dim(object1)[1]
      toplot1= cbind( rep(1,nfold) , c(1:nfold), object1[,2]) 
      toplot1= rbind( toplot1, cbind( rep(2,nfold) , c(1:nfold), object1[,4]) )
      toplot1= rbind( toplot1, cbind( rep(3,nfold) , c(1:nfold), object1[,6]) )                
      toplot1= rbind( toplot1, cbind( rep(4,nfold) , c(1:nfold), object1[,7]) )                                
      
      if        (type ==  "agree") {  object2 = object$lasso.agree.naive 
      } else if (type == "lincal") {  object2 = object$lasso.lincal.naive 
      } else if (type == "intcal") {  object2 = object$lasso.intcal.naive 
      } else if (type == "devian") {  object2 = object$lasso.devian.naive
      } else if (type == "devrat") {
        object2 = object$lasso.devian.naive
        for (j_ in c(1:7)) {
          object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
        }
      }
      toplot2 = cbind( 1 , 1, object2[2]) 
      toplot2 = rbind( toplot2, cbind( 2 , 1, object2[4]) )
      toplot2 = rbind( toplot2, cbind( 3 , 1, object2[6]) )                
      toplot2 = rbind( toplot2, cbind( 4 , 1, object2[7]) )                                
      
      object3 = colMeans(object1,na.rm=T)      
      if (type == "devrat") {
        for (j_ in c(1:7)) {
          object3[j_] = devratl[[j_]][[2]]
        }
      }
      toplot3 = cbind( 1 , 1, object3[2], 1) 
      toplot3 = rbind( toplot3, cbind( 2 , 1, object3[4], 1) ) 
      toplot3 = rbind( toplot3, cbind( 3 , 1, object3[6], 1) )
      toplot3 = rbind( toplot3, cbind( 4 , 1, object3[7], 1) ) 
      
      nms = c("Lasso", "Relaxed", "G0 Relax", "Ridge")
      i_ = 4 
    }
    
    toplot2
    toplot3
    
    ##### XGB ####################################################################
    if (object$fits[2] == 1) {
      
      if (is.null(object$xgb.devian.rep)) { object$xgb.devian.rep = object$xgb.devian.cv }
      if (is.null(object$xgb.intcal.rep)) { object$xgb.intcal.rep = object$xgb.intcal.cv }
      if (is.null(object$xgb.lincal.rep)) { object$xgb.lincal.rep = object$xgb.lincal.cv }
      if (is.null(object$xgb.agree.rep )) { object$xgb.agree.rep  = object$xgb.agree.cv  }
      
      if         (type == "agree") {  object1 = object$xgb.agree.rep 
      } else if (type == "lincal") {  object1 = object$xgb.lincal.rep 
      } else if (type == "intcal") {  object1 = object$xgb.intcal.rep 
      } else if (type == "devian") {  object1 = object$xgb.devian.rep 
      } else if (type == "devrat") {
        object1 = object$xgb.devian.rep
        devratl = list()
        for (j_ in c(1:6)) {
          devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
        }
      }
      nfold = dim(object1)[1]    
      
      if        (type ==  "agree") {  object2 = object$xgb.agree.naive 
      } else if (type == "lincal") {  object2 = object$xgb.lincal.naive 
      } else if (type == "intcal") {  object2 = object$xgb.intcal.naive 
      } else if (type == "devian") {  object2 = object$xgb.devian.naive 
      } else if (type == "devrat") {
        object2 = object$xgb.devian.naive
        for (j_ in c(1:6)) {
          object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
        }
      }
      
      object3 = colMeans(object1,na.rm=T)      
      if (type == "devrat") {
        for (j_ in c(1:6)) {
          object3[j_] = devratl[[j_]][[2]]
        }
      }

      if (xgbsimple == 1) {
        if ( sum(ensemble[c(1,5)]) >=1 ) { 
          i_= i_ + 1  
          toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,1]) )
          toplot2 = rbind( toplot2, cbind( i_ , 1, object2[1]) )
          toplot3 = rbind( toplot3, cbind( i_ , 1, object3[1], 1) )
          nms[i_] = "XGB"
        } 
        
        if ( sum(ensemble[c(2,6)]) >=1 ) { 
          i_ = i_ + 1 
          toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,2]) )
          toplot2 = rbind( toplot2, cbind( i_ , 1, object2[2]) )
          toplot3 = rbind( toplot3, cbind( i_, 1, object3[2], 2) )
          nms[i_] = "XGB feat."
        }
        
        if ( sum(ensemble[c(3,4,7,8)]) >=1 ) {   
          i_ = i_ + 1
          toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,3]) )
          toplot2 = rbind( toplot2, cbind( i_ , 1, object2[3]) )
          toplot3 = rbind( toplot3, cbind( i_, 1, object3[3], 4) )
          nms[i_] = "XGB offs."
        }
      }
      
      if ( sum(ensemble[c(1,5)]) >=1 ) { 
        i_= i_ + 1  
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,4]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[4]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[4], 1) )
        nms[i_] = "XGB Tune"
      } 
      
      if ( sum(ensemble[c(2,6)]) >=1 ) { 
        i_ = i_ + 1 
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,5]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[5]) )
        toplot3 = rbind( toplot3, cbind( i_, 1, object3[5], 2) )
        nms[i_] = "XGB T. fe"
      }
      
      if ( sum(ensemble[c(3,4,7,8)]) >=1 ) {   
        i_ = i_ + 1
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,6]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[6]) )
        toplot3 = rbind( toplot3, cbind( i_, 1, object3[6], 4) )
        nms[i_] = "XGB T. of"
      }
      frst = 0 
    }
    
    toplot2
    toplot3
    
    ##### Random Forest ##########################################################
    if (object$fits[3] ==1) {
      
      if (is.null(object$rf.devian.rep)) { object$rf.devian.rep = object$rf.devian.cv }
      if (is.null(object$rf.intcal.rep)) { object$rf.intcal.rep = object$rf.intcal.cv }
      if (is.null(object$rf.lincal.rep)) { object$rf.lincal.rep = object$rf.lincal.cv }
      if (is.null(object$rf.agree.rep )) { object$rf.agree.rep  = object$rf.agree.cv  }
      
      if         (type == "agree") {  object1 = object$rf.agree.rep 
      } else if (type == "lincal") {  object1 = object$rf.lincal.rep 
      } else if (type == "intcal") {  object1 = object$rf.intcal.rep 
      } else if (type == "devian") {  object1 = object$rf.devian.rep
      } else if (type == "devrat") {
        object1 = object$rf.devian.rep
        devratl = list()
        for (j_ in c(1:3)) {
          devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
        }
      }
      nfold = dim(object1)[1]    
      
      if        (type ==  "agree") {  object2 = object$rf.agree.naive 
      } else if (type == "lincal") {  object2 = object$rf.lincal.naive 
      } else if (type == "intcal") {  object2 = object$rf.intcal.naive 
      } else if (type == "devian") {  object2 = object$rf.devian.naive 
      } else if (type == "devrat") {
        object2 = object$rf.devian.naive
        for (j_ in c(1:3)) {
          object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
        }
      }
      
      object3 = colMeans(object1,na.rm=T)      
      if (type == "devrat") {
        for (j_ in c(1:3)) {
          object3[j_] = devratl[[j_]][[2]]
        }
      }
      
      if ( sum(ensemble[c(1,5)]) >=1 ) { 
        i_ = i_ + 1
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,1]) )
        toplot2 = rbind( toplot2, cbind(  i_ , 1, object2[1]) )
        toplot3 = rbind( toplot3, cbind(  i_ , 1, object3[1],1) )
        nms[i_] = "RF"
      }
      if ( sum(ensemble[c(2,6)]) >=1 ) { 
        i_ = i_ + 1
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,2]) )
        toplot2 = rbind( toplot2, cbind(  i_ , 1, object2[2]) )
        toplot3 = rbind( toplot3, cbind(  i_ , 1, object3[2],2) )  
        nms[i_] = "RF feat."
      }
      if ((sum(ensemble[c(3,4,7,8)]) >=1) & (family == "gaussian")) { 
        i_ = i_ + 1
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,3]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[3]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[3],4) )
        nms[i_] = "RF offs."
      }
      frst = 0 
    }
    
    toplot3
    
    ##### Oblique Random Forest ################################################
    if ( is.na( object$fits[8] ) ) { object$fits[8] = 0 }
    if (object$fits[8] ==1) {
      
      if (is.null(object$orf.devian.rep)) { object$orf.devian.rep = object$orf.devian.cv }
      if (is.null(object$orf.intcal.rep)) { object$orf.intcal.rep = object$orf.intcal.cv }
      if (is.null(object$orf.lincal.rep)) { object$orf.lincal.rep = object$orf.lincal.cv }
      if (is.null(object$orf.agree.rep )) { object$orf.agree.rep  = object$orf.agree.cv  }
      
      if         (type == "agree") {  object1 = object$orf.agree.rep 
      } else if (type == "lincal") {  object1 = object$orf.lincal.rep 
      } else if (type == "intcal") {  object1 = object$orf.intcal.rep 
      } else if (type == "devian") {  object1 = object$orf.devian.rep
      } else if (type == "devrat") {
        object1 = object$orf.devian.rep
        devratl = list()
        for (j_ in c(1:3)) {
          devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
        }
      }
      nfold = dim(object1)[1]    
      
      if        (type ==  "agree") {  object2 = object$orf.agree.naive 
      } else if (type == "lincal") {  object2 = object$orf.lincal.naive 
      } else if (type == "intcal") {  object2 = object$orf.intcal.naive 
      } else if (type == "devian") {  object2 = object$orf.devian.naive 
      } else if (type == "devrat") {
        object2 = object$orf.devian.naive
        for (j_ in c(1:3)) {
          object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
        }
      }
      
      object3 = colMeans(object1,na.rm=T)      
      if (type == "devrat") {
        for (j_ in c(1:3)) {
          object3[j_] = devratl[[j_]][[2]]
        }
      }
      
      if ( sum(ensemble[c(1,5)]) >=1 ) { 
        i_ = i_ + 1
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,1]) )
        toplot2 = rbind( toplot2, cbind(  i_ , 1, object2[1]) )
        toplot3 = rbind( toplot3, cbind(  i_ , 1, object3[1],1) )
        nms[i_] = "ORF"
      }
      if ( sum(ensemble[c(2,6)]) >=1 ) { 
        i_ = i_ + 1
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,2]) )
        toplot2 = rbind( toplot2, cbind(  i_ , 1, object2[2]) )
        toplot3 = rbind( toplot3, cbind(  i_ , 1, object3[2],2) )  
        nms[i_] = "ORF feat."
      }
      if ((sum(ensemble[c(3,4,7,8)]) >=1) & (family == "gaussian")) { 
        i_ = i_ + 1
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,3]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[3]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[3],4) )
        nms[i_] = "ORF offs."
      }
      frst = 0 
    }
    
    toplot3
    
    ##### Neural Network #########################################################
    if (object$fits[5] == 1) {
      
      if (is.null(object$ann.devian.rep)) { object$ann.devian.rep = object$ann.devian.cv }
      if (is.null(object$ann.intcal.rep)) { object$ann.intcal.rep = object$ann.intcal.cv }
      if (is.null(object$ann.lincal.rep)) { object$ann.lincal.rep = object$ann.lincal.cv }
      if (is.null(object$ann.agree.rep )) { object$ann.agree.rep  = object$ann.agree.cv  }

      if        (type == "agree" ) {  object1 = object$ann.agree.rep 
      } else if (type == "lincal") {  object1 = object$ann.lincal.rep 
      } else if (type == "intcal") {  object1 = object$ann.intcal.rep 
      } else if (type == "devian") {  object1 = object$ann.devian.rep 
      } else if (type == "devrat") {
        object1 = object$ann.devian.rep
        devratl = list()
        for (j_ in c(1:8)) {
          devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
        }
      }
      nfold = dim(object1)[1]    
      
      if        (type ==  "agree") {  object2 = object$ann.agree.naive 
      } else if (type == "lincal") {  object2 = object$ann.lincal.naive 
      } else if (type == "intcal") {  object2 = object$ann.intcal.naive 
      } else if (type == "devian") {  object2 = object$ann.devian.naive 
      } else if (type == "devrat") {
        object2 = object$ann.devian.naive
        for (j_ in c(1:8)) {
          object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
        }
      }
      
      object3 = colMeans(object1,na.rm=T)      
      if (type == "devrat") {
        for (j_ in c(1:8)) {
          object3[j_] = devratl[[j_]][[2]]
        }
      }
      
      for (j_ in c(1:8)) {
        if ( ensemble[j_]==1 ) { 
          i_ = i_ + 1
          if (j_ %in% c(1,5)) { colnum = 1
          } else if (j_ %in% c(2,6)) { colnum = 2
          } else { colnum = 4 } 
          toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,j_]) )
          toplot2 = rbind( toplot2, cbind( i_ , 1, object2[j_]) )    
          toplot3 = rbind( toplot3, cbind( i_ , 1, object3[j_], colnum) )
          nms[i_] = paste0("ANN ", j_) 
        }
      }
    }
    
    toplot3
    
    ##### RPART ##################################################################
    if (object$fits[4] == 1) {
      
      if (is.null(object$rpart.devian.rep)) { object$rpart.devian.rep = object$rpart.devian.cv }
      if (is.null(object$rpart.intcal.rep)) { object$rpart.intcal.rep = object$rpart.intcal.cv }
      if (is.null(object$rpart.lincal.rep)) { object$rpart.lincal.rep = object$rpart.lincal.cv }
      if (is.null(object$rpart.agree.rep )) { object$rpart.agree.rep  = object$rpart.agree.cv  }

      if        (type == "agree" ) {  object1 = object$rpart.agree.rep 
      } else if (type == "lincal") {  object1 = object$rpart.lincal.rep 
      } else if (type == "intcal") {  object1 = object$rpart.intcal.rep 
      } else if (type == "devian") {  object1 = object$rpart.devian.rep 
      } else if (type == "devrat") {
        object1 = object$rpart.devian.rep
        devratl = list()
        for (j_ in c(1:9)) {
          devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
        }
      }
      nfold = dim(object1)[1]    
      
      if        (type ==  "agree") {  object2 = object$rpart.agree.naive 
      } else if (type == "lincal") {  object2 = object$rpart.lincal.naive 
      } else if (type == "intcal") {  object2 = object$rpart.intcal.naive 
      } else if (type == "devian") {  object2 = object$rpart.devian.naive 
      } else if (type == "devrat") {
        object2 = object$rpart.devian.naive
        for (j_ in c(1:9)) {
          object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
        }
      }
      
      object3 = colMeans(object1,na.rm=T)      
      if (type == "devrat") {
        for (j_ in c(1:9)) {
          object3[j_] = devratl[[j_]][[2]]
        }
      }
      
      if ( sum(ensemble[c(1,5)]) >=1 ) { 
        i_ = i_ + 1 
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,1]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[1]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[1], 1) )
        nms[i_] = "RPART .00"
        i_ = i_ + 1     
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,2]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[2]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[2], 1) )
        nms[i_] = "RPART .01"
        i_ = i_ + 1     
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,3]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[3]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[3], 1) )
        nms[i_] = "RPART .02"
      }
      if ( sum(ensemble[c(2,6)]) >=1 ) { 
        i_ = i_ + 1 
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,4]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[4]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[4], 2) )
        nms[i_] = "RPa Fe .00"
        i_ = i_ + 1     
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,5]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[5]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[5], 2) )
        nms[i_] = "RPa Fe .01"
        i_ = i_ + 1     
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,6]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[6]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[6], 2) )
        nms[i_] = "RPa Fe .02"
      }
      
      if ((sum(ensemble[c(3,4,7,8)]) >= 1) & (!(family %in% c("binomial")))) { 
        i_ = i_ + 1 
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,7]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[7]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[7], 4) )
        nms[i_] = "RPa Of .00"
        i_ = i_ + 1     
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,8]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[8]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[8], 4) )
        nms[i_] = "RPa Of .01"
        i_ = i_ + 1     
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,8]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[9]) )
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[9], 4) )
        nms[i_] = "RPa Of .02"
      }
    }
    
    toplot3
    
    ##### Step Wise ##############################################################
    if ( (object$fits[6] == 1) | (object$fits[7] == 1) ) {
      if (is.null(object$step.devian.rep)) { object$step.devian.rep = object$step.devian.cv }
      if (is.null(object$step.lincal.rep)) { object$step.lincal.rep = object$step.lincal.cv }
      if (is.null(object$step.intcal.rep)) { object$step.intcal.rep = object$step.intcal.cv }
      if (is.null(object$step.agree.rep )) { object$step.agree.rep  = object$step.agree.cv  }
      
      object1 = object$step.agree.rep
      object2 = object$step.agree.naive
      object3 = object$step.agree.rep
      object3 = colMeans(object3,na.rm=T)
      
      if        (type == "agree" ) {  object1 = object$step.agree.rep 
      } else if (type == "lincal") {  object1 = object$step.lincal.rep 
      } else if (type == "intcal") {  object1 = object$step.intcal.rep 
      } else if (type == "devian") {  object1 = object$step.devian.rep
      } else if (type == "devrat") {
        object1 = object$step.devian.rep
        devratl = list()
        for (j_ in c(1:3)) {
          devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
        }
      }
      nfold = dim(object1)[1]    
      
      if        (type ==  "agree") {  object2 = object$step.agree.naive 
      } else if (type == "lincal") {  object2 = c(1,1,1)
      } else if (type == "intcal") {  object2 = c(0,0,0)
      } else if (type == "devian") {  object2 = object$step.devian.naive 
      } else if (type == "devrat") {
        object2 = object$step.devian.naive
        for (j_ in c(1:3)) {
          object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
        }
      }
      
      object3 = colMeans(object1,na.rm=T)      
      if (type == "devrat") {
        for (j_ in c(1:3)) {
          object3[j_] = devratl[[j_]][[2]]
        }
      }
    }
    
    #  if ( (object$fits[6] == 1) & (type %in% c("intcal","lincal"))) {
    if (object$fits[6] == 1) {
      i_ = i_ + 1 
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,1]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[1]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[1], 1) )
      nms[i_] = "Step df"
      
      i_ = i_ + 1 
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,2]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[2]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[2], 1) )
      nms[i_] = "Step p"
    }
    
    toplot3
    
    #  if ( (object$fits[7] == 1) & (type %in% c("intcal","lincal"))) {
    if  (object$fits[7] == 1) {
      i_ = i_ + 1 
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,3]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[3]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[3], 1) )
      nms[i_] = "AIC"
      aic_which = i_ 
    }
    
    toplot3
    
    cbind(nms,toplot2)
    
    ################################################################################
    ################################################################################
    
    if (is.null(ylim)) {  
      ylim = 1 
    } 
    if (length(ylim) == 1) {
      if (ylim == 0) { 
        limridge = 0
      } else { 
        ylim = 1 
        limridge = 1  
      }
    }
    
    if (length(ylim) != 2) {
      if (type == "agree") { 
        if (family == "gaussian") {
          ylim = c(0,1) 
        } else {
          ylim = c(0.5,1) 
        }
      } else if (type == "intcal") { 
        ylim = range(rbind(toplot1[,3],toplot2[,3],toplot3[,3]),na.rm=T)
        if (limridge == 1) {
          ylim2 = range(rbind(toplot1[(toplot1[,1]!=4),3],toplot2[ (toplot2[,1] != 4) ,3], toplot3[ (toplot3[,1]!=4) ,3]),na.rm=T) 
          ylim[2] = min( ylim[2], 2 * ylim2[2])
        }
      } else if (type == "lincal") { 
        ylim = range(rbind(toplot1[,3],toplot2[,3],toplot3[,3]), na.rm=T) 
        if (limridge == 1) {
          ylim2 = range(rbind(toplot1[(toplot1[,1]!=4),3],toplot2[ (toplot2[,1] != 4) ,3], toplot3[ (toplot3[,1]!=4) ,3]),na.rm=T) 
          ylim[2] = min( ylim[2], 2 * ylim2[2])
        }
        ylim[2] = max(ylim,1)
      } else if (type == "devian") { 
        if (family == "cox") {
          ylim = range(rbind(toplot1[,3],toplot3[,3]),na.rm=1) 
        } else {
          ylim = range(rbind(toplot1[,3],toplot2[,3],toplot3[,3]),na.rm=1) 
        }
        if ((limridge == 1) & (object$fits[7] == 1)) {
          ylim2 = range(rbind(toplot1[(toplot1[,1] != aic_which),3],toplot2[ (toplot2[,1] != aic_which) ,3], toplot3[ (toplot3[,1] != aic_which) ,3]),na.rm=1) 
          ylim[1] = max( ylim[1], 0.8 * ylim2[1])
          
          if (family == "cox") {
            ylim[2] = min( ylim[2], 1   * ylim2[2])
          } else {
            ylim[2] = min( ylim[2], 1.2 * ylim2[2])
          }
        }
        #      ylim[1] = 0 ; 
      } else if (type == "devrat") {
        ylim = range(rbind(toplot1[,3],toplot2[,3],toplot3[,3]),na.rm=1) 
        ylim[1] = 0 
      } 
    }
    
    #  print(ylim)
    
    if ((family == "cox") & (type == "intcal")) {
      cat(paste(" For the Cox model there is no intercept term and all calibration intercpet\n",
                " coefficients are by convention 0\n"))
    } else {
      if (plot != 0) {
        if (fold > 0) {
          lwd2 = 2 
          if (fold == 1 ) {
            if ( length(unique(toplot1[,2])) >=25) { lwd2 = 3 } else { lwd2 = 2 }
            interaction.plot(toplot1[,1], toplot1[,2], (toplot1[,3])^pow, col=toplot1[,2], 
                             ylab=ylab, xlab="",legend=F, ylim=ylim,
                             xaxt = "n")
            lines(toplot3[,1], toplot3[,3]^pow,col=1,lwd=lwd2) 
          } else if (fold > 1) {
            #            select fold random lines for plotting
            nunique = length( unique(toplot1[,2]) )
            if (fold < nunique) { 
              set.seed( object$seed$seedr[1] )
              smp = sample(nunique, fold, replace=0)
            } else { 
              smp = c(1:nunique) 
            }
            slct = toplot1[,2] %in% smp
            if ( fold >=25) { lwd2 = 3 } else { lwd2 = 2 }
            interaction.plot(toplot1[slct,1], toplot1[slct,2], (toplot1[slct,3])^pow, col=toplot1[,2], 
                             ylab=ylab, xlab="",legend=F, ylim=ylim,
                             xaxt = "n")
            lines(toplot3[,1], toplot3[,3]^pow,col=1,lwd=lwd2) 
          } else {
            lwd2 = 2
            plot(toplot3[,1], toplot3[,3]^pow, type='l', col=1, lwd=2, ylim=ylim, axes=0, xlab="", ylab=ylab) 
            axis(2)
            box()
          }
          axis(1, at = c(1:length(nms)), labels = nms, las=2)
          if (!((type %in% c("devian","devrat")) & (family == "cox"))) { 
            lines(toplot2[,1], toplot2[,3]^pow,col=2,lwd=lwd2) 
          }
          
          if (type == "lincal") { abline(h=1, lty=3, col=1, lwd=2) }
          if ((type == "intcal") | (type == "devrat")) { abline(h=0, lty=3, col=1, lwd=2) }
        } 
        if (fold == 0) {
          pch = toplot3[,4]
          pch = ifelse( pch == 1, 19, pch)
          pch = ifelse( pch == 2, 17, pch)
          pch = ifelse( pch == 4, 15, pch)
          #        print(pch)
          plot( toplot3[,1], (toplot3[,3])^pow, pch=pch, col=toplot3[,4], 
                ylab=ylab, xlab="",ylim=ylim,
                xaxt = "n")
          axis(1, at = c(1:length(nms)), labels = nms, las=2)
        }
      }
    }
    if (plot %in% c(0,2)) { return( list(toplot1=toplot1, toplot2=toplot2, toplot3=toplot3, nms=nms) ) }
    }
  }



