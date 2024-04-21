################################################################################

#  library(glmnetr)
#  library(splines)
#  library(survival)

################################################################################
################################################################################

#' Construct a single calibration plot for a nested.glmnetr output object 
#'
#' @param object1 A nested.glmnetr() output object for calibration
#' @param wbeta  Which Beta should be plotted, an integer.  This will depend on
#' which machine learning models were run when creating the output object.   If 
#' unsure the user can run the function without specifying wbeta and a legend 
#' will be directed to the console.  
#' @param df The degrees of freedom for the spline function, default=5. 
#' @param cv 1 to base the splines on the leave out X*Beta's ($xbetas.cv), or 0 to 
#' use the naive X*Beta's ($xbetas).
#' @param trim the percent of top and bottom of the data to be trimed way when 
#' producing plots.
#' @param fold the leave out fold to be selected for calibration or NULL to 
#' calibrate with all data.
#' @param collapse Collapse or combine pairs of hold out folds before dong the 
#' separate spline fits for calibration.  Intuitively it might make the spline 
#' fits more stable but empirically it doesn't seem to.
#' @param plot 1 by default to plot the calibration
#' @param lwd line width for plots, 1 by default
#' @param title title for inclusion on plots, default NULL 
#' @param xlim xlim for the plots.  This does not effect the curves within 
#' the plotted region.  Caution, for the "cox" framework the xlim are specified
#' in terms of the X*beta and not the HR, even when HR is described on the axes.
#' @param ylim ylim for the plots, which will usually only be specified in a 
#' second run of for the same data.  This does not effect the curves within 
#' the plotted region.  Caution, for the "cox" framework the ylim are specified
#' in terms of the X*beta and not the HR, even when HR is described on the axes.
#' @param xlab a user specified label for the x axis
#' @param ylab a user specified label for the y axis
#' @param col.term a number for the line depicting the overall calibration estimates
#' @param col.se a number for the line depicting the +/- 2 * standard error 
#' lines for the overall calibration estimates 
#' @param rug 1 to plot a rug for the model x*betas, 0 (default) to not.
#' @param plothr  a power > 1 determining the spacing of the values 
#' on the axes, e.g. 2, exp(1), sqrt(10) or 10.  The default of 0 plots the 
#' X*Beta.  This only applies for "cox" survival data models.   
#' @param ... allowance to pass terms to the invoked plot function
#'
#' @return Plots are returned by default, and optionally data for plots are 
#' output to a list.
#' 
#' @importFrom stats quantile 
#' @importFrom graphics box rug title
#' @importFrom survival coxph pspline
#' @importFrom splines ns
#' 
#' 
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
calplot0 = function(object1, wbeta=NULL, df=5, cv=1, trim=0, fold=NULL, collapse=0, plot=1, lwd=1, title=NULL, 
                    xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                    col.term=1, col.se=2, rug=1, plothr=0, ... ) {
# calplot0 = function(object1, wbeta=NULL, df=5, cv=1, trim=1, fold=NULL, collapse=0, plot=1, lwd=1, ... ) {
# wbeta = 5 ; df = 5 ; cv = 0 ; trim = 1 ; fold = 1 ;  plot = 1 ; lwd = 1 
# object1 = object ; trim=0 ; fold=i_ ; plot=plotfold_ ; title=paste("fold",i_) ;    
# fold=i_ ; plot=plotfold_ ; title=paste("fold",i_) ; lwd = 1 
  
# rug = 1   

  if ( (is.null(object1$xbetas.cv)) & (cv == 1) ) {
    (cv == 1) 
    cv = 0 
    cat(paste("  xbeta.cv is NULL so cv set to 0\n"))
  }
  
  stop_ = 0 
  if (is.null(wbeta)) {
    cat(paste(" specify num for wbeta = \n"))
    colnames(object1$xbetas.cv)[1] = "null"
    if (cv == 1) {  
      colnames(object$xbetas.cv)[1] = "null"
      varxb = diag(cov((object1$xbetas.cv)))
    } else { 
      colnames(object$xbetas)[1] = "null"
      varxb = diag(cov((object1$xbetas)))
    }
    print( data.frame(cbind( varxb, num=1:length(varxb))) )
    stop_ = 1 
  } 
  
  if (stop_ == 0) {
    
   family = object1$sample[1] 

    if (cv == 1) { origxb = object1$xbetas.cv[,wbeta ] ; xbname = colnames(object1$xbetas.cv)[wbeta] ;
    } else { origxb = object1$xbetas[,wbeta] ; xbname = colnames(object1$xbetas)[wbeta] ; }
   
    origy_ = object1$y_ ; 
    maxxb = max(origxb)
    minxb = min(origxb)
    c( minxb , maxxb ) 
    
    folds_n = as.numeric( object1$tuning[1] )
    if (!is.null(fold)) {
      foldid = object1$foldid 
      table( foldid )
      if ( (collapse == 1) & (folds_n%%2 == 0) ) {
        foldidnew = floor((foldid + 1)/2)
        table(foldidnew,foldid)  
        foldid = foldidnew   
      }
      foldy_ = origy_[ (foldid == fold) ]      
      foldxb = origxb[ (foldid == fold) ]
    } else {
      foldxb = origxb
      foldy_ = origy_ 
    }
#    foldxb
    minfoldxb = min(foldxb)
    maxfoldxb = max(foldxb)
    c( minfoldxb , maxfoldxb )
    
    ordr = order(foldxb)
    foldxb = foldxb[ordr]
    foldy_ = foldy_[ordr]
    if (family == "cox") {
      length(table(table(foldxb))) 
      lfoldxb = length(unique(foldxb))
      if ( lfoldxb < (df+1))  { df = lfoldxb - 1 }
      if (df > 1) {
        object2 = coxph(foldy_ ~ pspline(foldxb,df=df), data=as.data.frame(cbind(foldy_, foldxb)) ) 
      } else {
        object2 = coxph(foldy_ ~ foldxb, data=as.data.frame(cbind(foldy_, foldxb)) ) 
      }
    } else {
      object2 = glm(foldy_ ~ ns(foldxb,df=df) , family=family ) # data=as.data.frame(cbind(foldy_, xb))  
      summary(object2)
    }

    plotxb = minxb + rep(0:100)*(maxxb-minxb)/100

    ## mean(predict(object2)) ## 9e-16
    meancalxb =  mean(predict(object2,data.frame(foldxb=origxb)))
    #      cat( paste(title, meancalxb, "\n" ) )  ## 9e-16
    meanfoldxb = mean(foldxb)
    dfr = data.frame(foldxb=plotxb)
    est = predict(object2, newdata=dfr)   
    predicteds = predict(object2, newdata=data.frame(foldxb=plotxb), se=TRUE)  
    est = predicteds$fit
    se = predicteds$se.fit
    estci = cbind(plotxb, est, se, est - 2 * se, est + 2 * se)
    
    lty1 = 1 ; lty2 = 2 ;   
    if (is.null(ylab)) {
      if (family %in% c("cox") ) {
        if (plothr > 1) {
          ylab = "Calibrated HR"
        } else {
          ylab = "Calibrated log(HR)"
        }
      } else if (family %in% c("binomial","gaussian") ) {
        ylab = "Calibrated X*Beta"
      } 
    }
    if (is.null(xlab)) {
      if (family %in% c("cox") ) {
        if (plothr > 1) {
          if (cv==1) { 
            xlab = paste(xbname,"HR (CV)") 
          } else { 
            xlab = paste(xbname,"HR (naive)")
          }
        } else {
          if (cv==1) { 
            xlab = paste(xbname,"log(HR) (CV)") 
          } else { 
            xlab = paste(xbname,"log(HR) (naive)")
          }
        }
      } else if (family %in% c("binomial","gaussian") ) {
        if (cv==1) { xlab = paste(xbname,"X*Beta (CV)") } else { xlab = paste(xbname,"X*Beta (naive)") }
      } 
    }
    xlims = quantile(foldxb, (c(trim,(100-trim))/100))
#    print(xlims)
    summary(estci[,1])
#    estci_ = estci[((estci[,1]>=xlims[1]) & (estci[,1]<= xlims[2])),]
#    summary(estci_[,1])
    if (is.null(xlim)) { xlim = c(minxb, maxxb) }
    if (is.null(ylim)) { ylim = c(min(estci[,2]), max(estci[,2])) }
    if (plot %in% c(1,2)) {
        if ((minxb < minfoldxb) | (maxxb > maxfoldxb)) { 
          indxlo = c(1:length(estci[,1])) [ estci[,1] < minfoldxb ] 
          indxhi = c(1:length(estci[,1])) [ estci[,1] > maxfoldxb ] 
          if ( length( c(indxlo,indxhi) ) > 0) {
            indxfold = c(1:length(estci[,1])) [-c(indxlo,indxhi) ] 
          } else {
            indxfold = c(1:length(estci[,1]))
          }
          plot(  estci[indxlo,1], estci[indxlo,2] , type="l", 
                 lty=lty2, col=col.term, lwd=lwd, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=0, ... ) 
          lines( estci[indxlo,1], estci[indxlo,4], lty=lty2, lwd=lwd, col=col.se ) 
          lines( estci[indxlo,1], estci[indxlo,5], lty=lty2, lwd=lwd, col=col.se ) 
          lines( estci[indxfold,1], estci[indxfold,2], lty=lty1, lwd=lwd, col=col.term ) 
          lines( estci[indxfold,1], estci[indxfold,4], lty=lty1, lwd=lwd, col=col.se ) 
          lines( estci[indxfold,1], estci[indxfold,5], lty=lty1, lwd=lwd, col=col.se ) 
          lines( estci[indxhi,1], estci[indxhi,2], lty=lty2, lwd=lwd, col=col.term ) 
          lines( estci[indxhi,1], estci[indxhi,4], lty=lty2, lwd=lwd, col=col.se ) 
          lines( estci[indxhi,1], estci[indxhi,5], lty=lty2, lwd=lwd, col=col.se ) 
        } else {
         plot(  estci[,1], estci[,2] , type="l", col=col.term, lwd=lwd, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=0, ... ) 
         lines( estci[,1], estci[,4] , lty=lty1, lwd=lwd, col=col.se ) 
         lines( estci[,1], estci[,5] , lty=lty1, lwd=lwd, col=col.se )
        }
      
         if (family == "cox") {
           abline(a=-meanfoldxb, b=1, lty=3)
           abline(h=0, lty=1)
           myaxis(ylim, plothr, est , side=2)
           myaxis(xlim, plothr, plotxb, side=1)
         } else { 
           abline(a=0, b=1, lty=1)
           axis(side=2)
           axis(side=1)
         }
         box() 
         
         keep_  = ((foldxb >= xlim[1]) & (foldxb <= xlim[2]))
         if (sum(keep_==0) > 0) {
           cat(paste("  Due to user specified xlim ", sum(keep_==0) , "tick marks are not displayed (fold", fold, ")\n"))
         }
         foldxb_ = foldxb[ keep_ ]
         foldy__ = foldy_[ keep_ ]
         if (rug == 1) { myrug(foldxb_, foldy_, family) } 
         if (!is.null(title)) { title ( title ) }
#        termplot(object2, terms=1,se=TRUE, xlab=xlab, ylab="Calibrated log(HR)", lwd.term=lwd, lwd.se=lwd )
    }
    if (plot %in% c(0,2)) {
      return( list( estci=estci, meancalxb=meancalxb, meanfoldxb=meanfoldxb, minmax=c(minfoldxb,maxfoldxb) ) )
    }
  }
}
#
################################################################################
#
#' Construct calibration plots for a nested.glmnetr output object 
#' 
#' Using k-fold cross validation this function constructs calibration plots for 
#' a nested.glmnetr output object. Each hold out subset of the k-fold cross 
#' validation is regressed on the x*beta predicteds based upon the 
#' model fit using the non-hold out data using splines.  This yields k spline 
#' functions for evaluating model performance.  These k spline functions are 
#' averaged to provide an overall model calibration.  Standard deviations of 
#' the k spline fits are also calculated as a function of the predicted X*beta, 
#' and these are used to derive and plot approximate 95% confidence intervals 
#' (mean +/- 2 * SD/sqrt(k)).  Because regression equations can be unreliable 
#' when extrapolating beyond the data range used in model derivation, we 
#' display this overall calibration fit and CIs with solid lines only for the 
#' region which lies within the ranges of the predicted x*betas for  
#' all the k leave out sets.  The spline fits are made using the same framework 
#' as in the original machine learning model fits, i.e. one of "cox", "binomial" 
#' or "gaussian"family. For the "cox" famework the pspline() funciton is used, 
#' and for the "binomial" and "gaussian" frameworks the ns() function is 
#' used.  Predicted X*betas beyond the range of any of the hold 
#' out sets are displayed by dashed lines to reflect the lessor certainty when 
#' extrapolating even for a single hold out set.   
#' 
#' Optionally, for comparison,  
#' the program can fit a spline based upon the predicted x*betas ignoring the 
#' cross validation structure, or one can fit a spline using the x*betas 
#' calculated using the model based upon all data. 
#'
#' @param object A nested.glmnetr() output object for calibration
#' @param wbeta  Which Beta should be plotted, an integer.  This will depend on
#' which machine learning models were run when creating the output object.   If 
#' unsure the user can run the function without specifying wbeta and a legend 
#' will be directed to the console.  
#' @param df The degrees of freedom for the spline function 
#' @param cv 1 to base the splines on the leave out X*Beta's ($xbetas.cv), or 0 to 
#' use the naive X*Beta's ($xbetas).
#' @param trim the percent of top and bottom of the data to be trimmed away when 
#' producing plots.  The original data are still used used calcualting the curves 
#' for plotting. 
#' @param vref Similar to trim but instead of trimming the spline lines, plots 
#' vertical refence lines aht the top vref and bottom vref percent of the model
#' X*Betas's
#' @param fold 1 (default) to base the calibration by first fitting splines for 
#' each individual fold and then averaging, 0 to base the calibration on a single
#' spline fit using all X*Beta. 
#' @param collapse 1 to collapse or combine pairs of hold out folds before doing 
#' the separate spline fits for calibration.  Intuitively it might make the spline 
#' fits more stable but empirically it doesn't seem to.
#' @param plot 1 by default to produce plots, 0 to output data for plots only, 
#' 2 to plot and output data. 
#' @param plotfold 0 by default to not plot the individual fold calibrations, 1 
#' to produce pltos for each of the k hold out sets, and 2 to overlay the k 
#' leave out fits to the overall calibration curve. 
#' @param xlim xlim for the plots.  This does not effect the curves within 
#' the plotted region.  Caution, for the "cox" framework the xlim are specified
#' in terms of the X*beta and not the HR, even when HR is described on the axes.
#' @param ylim ylim for the plots, which will usually only be specified in a 
#' second run of for the same data.  This does not effect the curves within 
#' the plotted region.  Caution, for the "cox" framework the ylim are specified
#' in terms of the X*beta and not the HR, even when HR is described on the axes.
#' @param xlab a user specified label for the x axis
#' @param ylab a user specified label for the y axis
#' @param col.term a number for the line depicting the overall calibration estimates
#' @param col.se a number for the line depicting the +/- 2 * standard error 
#' lines for the overall calibration estimates 
#' @param rug 1 to plot a rug for the model x*betas, 0 (default) to not.
#' @param plothr  a power > 1 determining the spacing of the values 
#' on the axes, e.g. 2, exp(1), sqrt(10) or 10.  The default of 0 plots the 
#' X*Beta.  This only applies fore "cox" survival data models.  
#' @param ... allowance to pass terms to the invoked plot function  
#'
#' @return Calibration plots are returned by default, and optionally data for plots 
#' are output to a list.
#' 
#' @seealso
#'   \code{\link{plot.nested.glmnetr}} , \code{\link{summary.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
#'   
#' @author Walter Kremers (kremers.walter@mayo.edu)
#'
#' @export
#'
calplot = function(object, wbeta=NULL, df=3, cv=1, trim=0, vref=0, fold=1, collapse=0, 
                   plot=1, plotfold=0, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                   col.term=1, col.se=2, rug=1, plothr=0, ... ) {
  # wbeta = 5 ; df = 5 ; cv = 1 ; trim = 0 ; vref=0 ; fold = 1 ;  collapse=0 ; 
  # plot = 1 ; plotfold=0 ; xlim=NULL ;  ylim=NULL ; xlab=NULL ;  ylab=NULL ; 
  # col.term=1 ;  col.se=2 ;  plothr=0 ; 
  
  # object = nps_csf_glo_ratcap_fit
  # object = nested_glmnetr_fit_hp_las
  # object1 = object ;
  
  #  object1 = object ;  lwd=2 ;  trim = 0 ; fold=NULL ;  
  #  wbeta=5 ; cv=0  ; fold=0 ; plotfold=0 ; ylim=c(-0.2,0.1) ; 
  
  if ( plothr == 1 ) { plothr = exp(1) }
  
  if (is.null(object$xbetas.cv)) {
    cv = 0 
    cat(paste("  xbeta.cv is NULL so cv set to 0\n"))
  }

  stop_ = 0 
  if (is.null(wbeta)) {
    cat(paste(" specify num for wbeta = \n"))
    if (cv == 1) {  
      colnames(object$xbetas.cv)[1] = "null"
      Var = diag(cov((object$xbetas.cv)))
    } else { 
      colnames(object$xbetas)[1] = "null"
      Var = diag(cov((object$xbetas)))
    }
    print( data.frame(cbind( Var, num=1:length(Var))) )
    stop_ = 1 
  } 
  
  if (stop_ == 0) {
    
    family = object$sample[1] 
    
    if (cv == 1) { xbeta = object$xbetas.cv[,wbeta ] ; xbname = colnames(object$xbetas.cv)[wbeta] ;
    } else { xbeta = object$xbetas[,wbeta] ; xbname = colnames(object$xbetas)[wbeta] ; }
#    print(summary(xbeta))
    y_ = object$y_ ; 
    
    if (var(xbeta) == 0) {
      warning("  There may be no varaibalbity in the model X*Beta's\n",
              "    The calibration plot may be nonsensical")
    }
    
    maxxb = max(xbeta)
    minxb = min(xbeta)
    if (is.null(xlim)) {  xlim = c( minxb, maxxb )  }
    
    if (is.null(ylab)) {
      if (family == "cox") {
        if (plothr > 1) { ylab = "Calibrated HR"
        } else { ylab = "Calibrated log(HR)" } 
      } else {
        ylab = "Calibrated X*Beta"
      }
    }
    if (is.null(xlab)) {
      if (family == "cox") {
        if (plothr > 1) {
          if (cv==1) { xlab = paste(xbname,"HR (CV)") } else { xlab = paste(xbname,"HR (naive)") }
        } else {
          if (cv==1) { xlab = paste(xbname,"log(HR) (CV)") } else { xlab = paste(xbname,"log(HR) (naive)") }
        }
      } else {
        if (cv==1) { xlab = paste(xbname,"X*Beta (CV)") } else { xlab = paste(xbname,"X*Beta (naive)") }
      }
    }
    
#    print(xlab) 
#    print(ylab)

    if (is.null(fold)) {
      if (family == "cox") { fold = 0 } else { fold = 1 }
      fold = 1 
    }
    
    folds_n = as.numeric( object$tuning[1] )
    if (!(plotfold %in% c(0,1,2))) { plotfold = 0 }                             ## save estimates
    plotfold_ = plotfold 
    if (plotfold == 1) { plotfold_ = 2 }                                        ## make individual plots
    if (plotfold == 2) { plotfold_ = 0 }                                        ## make single plot with folds overlaid 
    
    if (cv == 1) { xbname = colnames(object$xbetas.cv)[wbeta] 
    } else { xbname = colnames(object$xbetas)[wbeta] }
    
    if ( (collapse == 1) & (folds_n%%2 == 0) ) { folds_n = folds_n / 2 }
    
    if (fold == 1) {
      est.cv = matrix(0,nrow=folds_n, ncol=101)
      se.cv = est.cv
      lower.cv = est.cv
      upper.cv = est.cv
      meancalxb  = rep(0,folds_n)
      meanfoldxb = rep(0,folds_n)
      minplot = min(minxb)
      maxplot = max(maxxb)
      for (i_ in c(1:folds_n)) { 
        calplotout = calplot0(object, wbeta=wbeta, df=df, cv=cv, trim=trim, fold=i_ , 
                     collapse=collapse, plot=plotfold_, col.term=col.term, col.se=col.se, 
                     title=paste("fold",i_), xlim=xlim, ylim=ylim, plothr=plothr, rug=rug, ... )
        ## get back the xbeta from the whole data set -- no predicteds used
        if (i_ == 1) { plotxb = calplotout$estci[,1] } 
        est.cv[i_,]  = calplotout$estci[,2]
        se.cv[i_,]   = calplotout$estci[,3]
        lower.cv[i_,] = calplotout$estci[,4]
        upper.cv[i_,] = calplotout$estci[,5]
        
        meancalxb[i_] = calplotout$meancalxb 
        meanfoldxb[i_] = calplotout$meanfoldxb 
        minplot = max( minplot, calplotout$minmax[1] ) 
        maxplot = min( maxplot, calplotout$minmax[2] ) 
      }
#      plotxb.cv = calplotout$estci[,1]
      c(minplot, maxplot)
      ## get average and SE from the folds_n folds 
      est = colMeans(est.cv)
      se = sqrt(apply(est.cv,2,var)/folds_n)    
      upper = est + 2 * se 
      lower = est - 2 * se 
      
      if (plot %in% c(1,2)) {
        lty1 = 1 ; lty2 = 2 ; 
        
        xlims = quantile(xbeta, (c(trim,(100-trim))/100))
        keep_ = ( ((plotxb >= xlims[1]) & (plotxb <= xlims[2])) & ((plotxb >= xlim[1]) & (plotxb <= xlim[2])) ) 
        if ( trim > 0) { cat("  Spline curves are trimed at x's of", xlims[1], "and", xlims[2], ", or", trim, "% hi and lo\n") } 
#        print(summary(xbeta))
#        print(summary(plotxb))
#        print(table(keep_))
        
        plotxb_ = plotxb[keep_]
        est_ = est[keep_]
        lower_ = lower[keep_]
        upper_ = upper[keep_]
        
        if (vref > 0) { vrefs = quantile(xbeta, (c(vref,(100-vref))/100)) }
        
        if (plotfold == 2) { lwd = 3 } else { lwd = 2 }

#          plot(plotxb_, est_, type="l", lty=lty1, col=col.term, lwd=lwd, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=0, ... )
#          lines( plotxb_, lower_, lty=lty1, col=col.se, lwd=lwd) 
#          lines( plotxb_, upper_, lty=lty1, col=col.se, lwd=lwd) 
          
          indxlo = c(1:length(plotxb_)) [ plotxb_ < minplot ] 
          indxhi = c(1:length(plotxb_)) [ plotxb_ > maxplot ] 
          if ( length( c(indxlo,indxhi) ) > 0) {
            indxfold = c(1:length(plotxb_)) [-c(indxlo,indxhi) ] 
          } else {
            indxfold = c(1:length(plotxb_))
          }
          
          xlim0 = c(min(plotxb_), max(plotxb_)) 
          if (!is.null(xlim)) { xlim0 = xlim }
          if (is.null(ylim)) {
            ylim = c(min(est_), max(est_))
            if ( min(lower_[indxfold] ) < ylim[1] ) { ylim[1] = min(lower_[indxfold] ) }
            if ( max(upper_[indxfold] ) > ylim[2] ) { ylim[2] = max(upper_[indxfold] ) }
          }
           
          plot(  plotxb_[indxlo], est_[indxlo] , type="l", 
                 lty=lty2, col=col.term, lwd=lwd, xlab=xlab, ylab=ylab, xlim=xlim0, ylim=ylim, axes=0, ... ) 
          lines( plotxb_[indxlo]  , lower_[indxlo  ], lty=lty2, lwd=lwd, col=col.se ) 
          lines( plotxb_[indxlo]  , upper_[indxlo  ], lty=lty2, lwd=lwd, col=col.se ) 
          lines( plotxb_[indxfold], est_  [indxfold], lty=lty1, lwd=lwd, col=col.term ) 
          lines( plotxb_[indxfold], lower_[indxfold], lty=lty1, lwd=lwd, col=col.se ) 
          lines( plotxb_[indxfold], upper_[indxfold], lty=lty1, lwd=lwd, col=col.se ) 
          lines( plotxb_[indxhi]  , est_  [indxhi  ], lty=lty2, lwd=lwd, col=col.term ) 
          lines( plotxb_[indxhi]  , lower_[indxhi  ], lty=lty2, lwd=lwd, col=col.se ) 
          lines( plotxb_[indxhi]  , upper_[indxhi  ], lty=lty2, lwd=lwd, col=col.se ) 
          
        if (family == "cox") { 
          abline(a=-mean(xbeta), b=1, lty=3)
          abline(h=0, lty=1)
          myaxis(ylim, plothr, est)
          myaxis(xlim, plothr, NULL, side=1)
        } else { 
          abline(a=0, b=1, lty=3)
          axis(side=1)
          axis(side=2)
        }
        box()
        keep_  = ((xbeta >= xlim[1]) & (xbeta <= xlim[2]))
#        length( xbeta ) 
#        length(keep_)
#        sum(keep_)
        if (sum(keep_==0) > 0) {
          cat(paste("  Due to user specified xlim ", sum(keep_==0) , "tick marks are not displayed in rug\n"))
        }
        xbeta_ = xbeta[ keep_ ]
#        length(xbeta_)
        y__    = y_[ keep_ ]
#        c( min(xbeta), min(xbeta_) ) 
        if (rug == 1) { myrug(xbeta_, y__, family) } 
        if (vref > 0) { abline(v=vrefs, lty=3, col=2) }
      }
      if (plotfold == 2) {
        for (i_ in c(1:folds_n)) {  
          lines( plotxb, est.cv[i_,], lty=i_, col=i_) 
        }
      }
    }
    
    if (fold == 0) {
      estci = calplot0(object, wbeta=wbeta, df=df, cv=cv, trim=trim, lwd=2, plot=plot, xlim=xlim, ylim=ylim, 
                       col.term=col.term, col.se=col.se, xlab=xlab, ylab=ylab, plothr=plothr, rug=rug, ... )$estci  
    } 
    
    if (plot %in% c(0,2)) {
      if (fold == 1) {  
        return_list = list(estimates=cbind(plotxb, est, se, lower, upper), 
                           est.cv=est.cv, se.cv=se.cv, lower.cv=lower.cv, upper.cv=upper.cv)  
      } else {
        return_list = list( estci = estci )  
      }
      return( return_list )
    }
  }
} 

#===============================================================================
# x_ = rnorm(100)^2
# y_ = rbinnom(100,1,(1/(1+exp(x_)))) 

#' A customized rug 
#'
#' @param x_ the x variable to be included in the rug
#' @param y_ an event variable for placing the x_ values top for events and bottom
#' for non-events 
#' @param family one of "cox", "binomial" or "gaussian" 
#'
#' @return A rug to plot 
#'
myrug = function(x_, y_=NULL, family="gaussian") {
  if (family == "cox") {
    if (dim(y_)[2]== 3) { e_ = y_[,3] 
    } else { e_ = y_[,2] }
#    ticksize = (0.1*y_[,1]/max(y_[,1])) 
#    ticksize[ticksize < 0.001] = 0.001
#    print(c(length(ticksize), length(y_)))
#    print(summary(ticksize))
#    ticksize0 = ticksize[e_==0]
#    ticksize1 = ticksize[e_==1]
    ticksize0 = 0.03 
    ticksize1 = 0.03 
    rug(x_[e_==0], side=1, col=1, ticksize=ticksize0)
    rug(x_[e_==1], side=3, col=1, ticksize=ticksize1)
  } else if (family == "binomial") {
    e_ = y_ 
    rug(x_[e_==0], side=1, col=1,ticksize=0.03)
    rug(x_[e_==1], side=3, col=1,ticksize=0.03)
  } else {
    rug(x_, side=1, col=1, ticksize=0.03)
  }
}

#===============================================================================
#' Un-log the log(HR)'s for plotting
#'
#' @param zlim the y limits on the log(HR) scale A nested.glmnetr() output object for calibration
#' @param plothr  a power of determining the placement of y axis 
#' values, e.g. 2, sqrt(10) or 10.  The default of 0 plots the X*Beta.
#' @param est values of X*Beta estimates to be plotted which will determine the 
#' range to plot in case zlim is NULL. 
#' @param side 1 for horizontal axis, 2 for vertical 
#' 
#' @return An axis statement 
#'  
myaxis = function(zlim,plothr,est=NULL,side=2) {
  if (is.null(plothr)) { plothr = 0 }
  if ((is.null(zlim)) & (!is.null(est))) {
    zlim = c(min(est),max(est)) 
  } else if ((is.null(zlim)) & (is.null(est))) {
    plothr = 0
  }
  if (plothr > 1) {
    xbs = zlim
    hrs = exp(zlim) 
    logXhrs = log(hrs) / log(plothr)
    logXhrs_ = c(ceiling((logXhrs[1])): floor(logXhrs[2]))
    hrs_ = plothr^(logXhrs_)
    loghrs_ = log(hrs_)
    labls = round(hrs_,digits=2) 
    myaxis = axis(side, at=loghrs_, labels = labls, las=1 )  
  } else {
    myaxis = axis(side)  
  }
  return(myaxis)
}

################################################################################
################################################################################

  
  
