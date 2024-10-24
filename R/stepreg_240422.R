################################################################################
##### stepreg_yymmdd.R #########################################################
################################################################################
#' Fit the steps of a stepwise regression.  
#'
#' @param xs_st predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
#' the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
#' and not a data frame. 
#' @param start_time_st start time, Cox model only - class numeric of length same as number of patients (n)
#' @param y_st output vector: time, or stop time for Cox model, y_st 0 or 1 for binomal (logistic), numeric for gaussian. 
#' Must be a vector of length same as number of sample size. 
#' @param event_st event_st indicator, 1 for event, 0 for census, Cox model only.
#' Must be a numeric vector of length same as sample size.
#' @param steps_n number of steps done in stepwise regression fitting 
#' @param method method for choosing model in stepwise procedure, "loglik" or "concordance".
#' Other procedures use the "loglik".     
#' @param family model family, "cox", "binomial" or "gaussian" 
#' @param track 1 to output stepwise fit program, 0 (default) to suppress 
#'
#' @return does a stepwise regression of depth maximum depth steps_n
#' 
#' @seealso 
#'    \code{\link{summary.stepreg}} , \code{\link{aicreg}} , \code{\link{cv.stepreg}} , \code{\link{nested.glmnetr}}   
#' 
#' @export
#' 
#' @importFrom stats binomial glm.fit 
#'
#' @examples
#' set.seed(18306296)
#' sim.data=glmnetr.simdata(nrows=100, ncols=100, beta=c(0,1,1))
#' # this gives a more intersting case but takes longer to run
#' xs=sim.data$xs           
#' # this will work numerically
#' xs=sim.data$xs[,c(2,3,50:55)] 
#' y_=sim.data$yt  
#' event=sim.data$event
#' # for a Cox model 
#' cox.step.fit = stepreg(xs, NULL, y_, event, family="cox", steps_n=40) 
#' # ... and for a linear model 
#' y_=sim.data$yt  
#' norm.step.fit = stepreg(xs, NULL, y_, NULL, family="gaussian", steps_n=40) 
#' 
stepreg = function(xs_st, start_time_st=NULL, y_st, event_st, steps_n=0, method="loglik", family=NULL, track=0) {
  
#  if  ( is.null(start_time_st)) { print("HERE step 0 : start_time_st=NULL") }
  
#  xx1 = start_time_st[1]
#  xx2 = start_time_st[2]
#  if  (!is.null(start_time_st)) { print("HERE step 1 : start_time_st=", start_time_st[1], start_time_st[2], start_time_st[3], start_time_st[4]) }
#  if  (!is.null(start_time_st)) { print("HERE step 1 : start_time_st=", xx1, xx2) }
  
  ## xs_st or predictors - must be as matrix 
  xsk = dim(xs_st)[2]                                                              ## number of candidate predictor variables 
  nobs = dim(xs_st)[1] 
  one = c( rep(1, nobs) )
  
  if (steps_n==0) { steps_n=xsk }
  steps_n = min(steps_n, xsk)
  
  if ( is.null(family)) {family = NA}
  if ( is.na(family) & (!is.null(start_time_st)) ) { family = "cox" }
  if ( is.na(family) & (!is.null(event_st)) & is.null(y_st)) { 
    y_st = event_st 
    event_st = NULL 
    family = "binomial" 
  } else if ( is.na(family) & (!is.null(event_st     )) ) { 
    family = "cox" 
  }
  if ( is.na(family) & ( length(table(y_st)) == 2 ) ) { family = "binomial" }
  if ( is.na(family) & ( length(table(y_st)) >  2 ) ) { family = "gaussian" }
  if ( family == "poison" ) { family = "poisson" }
  if ( family %in% c("gausian", "gause", "normal") ) { family = "gaussian" }
  if ( family %in% c("cox", "Cox", "coxph") ) { family = "cox" }

  coxcontrol = survival::coxph.control() 
  glmcontrol = glm.control()

  riskvarsp =  c( rep(0, xsk)  )                                                ## iniitialize risk factors for model with no terms - presence indicator 
  
  xsnames = colnames(xs_st)
  dframe   = as.data.frame( cbind(y_st, xs_st) ) 
#  colnames(dframe)
  
  betalast = NULL 
# i_ = 1 ; j_ = 1 ;   ## for testing only 
# i_ = 2 ; j_ = 1 ; 
  ##### begin i_ loop ##########################################################
  if (track>=1) { cat ("Calculating stepwise model step : ") }
  for (i_ in 1:steps_n)  { 
#  for (i_ in 1:10)  { 
    if (track>=1) { cat ( i_, " ") } 
    for (j_ in 1:xsk) {
      if (riskvarsp[j_] == 0) {
        riskvarsnewp = as.numeric( riskvarsp ) ;                                ## initialize risk factor for new model (present indicator) from last best fit
        riskvarsnewp[j_] = 1                                                    ## add current variable to risk factors for new model fit (present indicators)
        riskvarsnewp 

        sum( riskvarsnewp )
        dim(xs_st)
        
        riskm = as.matrix( xs_st[,(riskvarsnewp==1)] )                          ## risk factors for model fit in current step 
        riskm [ 1:20,]
        dim(riskm)

        riskvarsnew = xsnames[riskvarsnewp==1]
        riskvarsnew
      
##        if ( (i_==1) & (j_ %in% c((0:20)*10+1) ) ) cat ("HERE start_time_st = ", start_time_st[1:10], "\n" )
        if (family == "cox") {
          if (is.null(betalast)) { init = rep(0,dim(riskm)[2]) 
          } else { init = betalast[ (riskvarsnewp==1) ] ; length(init) } 
          if (is.null(start_time_st)) { 
            fit1 = survival::coxph.fit(riskm ,Surv(                y_st,  event_st), strata=NULL, init=init,
                                       control=coxcontrol, method="efron", rownames=NULL, resid=FALSE)  
            clist =  concordancefit(Surv(       y_st, event_st), -fit1$linear.predictors)   
          } else {  
            init = c(rep(0,dim(riskm)[2]))
            fit1 = survival::agreg.fit(riskm ,Surv( start_time_st, y_st, event_st), strata=NULL, init=init,  
                                      control=coxcontrol, method="efron", rownames=NULL, resid=FALSE) 
            clist = concordancefit(Surv(start_time_st, y_st, event_st), -fit1$linear.predictors)     ## get concordances  --  allows extension to shose best by C instead of loglik
#           if ( (i_==1) & (j_ %in% c(0:5) ) ) cat ("HERE step 0.2" , "\n" )
          }
          if (rankMatrix(riskm) < 1e-12) { fit1$coefficients = rep(0,dim(riskm)[2]) } 
          fit1$coefficients[is.na(fit1$coefficients)] = 0                       # assign missing coefficients value 0          
          loglik = fit1$loglik 
          names(loglik) = c("loglik.null", "loglik")
          concord = c(as.numeric( clist[1] ), sqrt( as.numeric( clist[4] ) ) )    ## check this XXX
          length( concord )
          if (length( concord ) == 2) { 
            names(concord) = c("concordance", "std") 
          } else {  
            print(" HERE step 2" )  
            print(  names(clist))  
          }
        } else if ( family == "binomial" ) {
          riskm1 = cbind(one, riskm)                                                               ## c(rep(1,dim(y_st)[1]))  ## start=c(0) 
          fit1 = glm.fit(riskm1 , y_st , family = binomial(link = "logit"), control = glmcontrol)  ## , start = beta_val , weights=weights_val_j
          summary(fit1)
#          beta = fit1$coefficients 
#          beta[is.na(beta)] = 0 
#          phat = riskm1 %*% beta # fit1$coefficients
#          table(phat)
          fit1$coefficients[is.na(fit1$coefficients)] = 0                       # assing missing coefficients value 0
#          if (rankMatrix(riskm1) < 1e-12) { fit1$coefficients = rep(0,dim(riskm1)[2]) } 
#          clist = concordancefit( y_st, fit1$linear.predictors)  
#          form1 = formula( paste(" y_st ~ " , paste(riskvarsnew, collapse = " + " ) ) )  
#          fit1 = glm( form1, data=dframe, family = binomial(link = "logit"))   ## , start = beta_val , weights=weights_val_j          
          clist =  concordance( y_st ~ fit1$linear.predictors )                 ## get concordances  
          concord = c(concordance=as.numeric( clist[1] ), se=sqrt( as.numeric( clist[4] ) ) )  
          loglik = c(loglik.null=-fit1$null.deviance/2, loglik=-fit1$deviance/2)
        } else if ( family == "gaussian" ) {
          ##          riskm1 = cbind(1, riskm)                                                    ## c(rep(1,dim(y_st)[1]))  ## start=c(0) 
          ##          fit1 = glm.fit(riskm1 , y_st , family = binomial(link = "logit"), control = glmcontrol)  ## , start = beta_val , weights=weights_val_j
          ##          fit1 = glm(event_st ~ riskm1 , family=family)                              ######## your model statement here #######
          ##          ##        (var(fit1$linear.predictors) + var(fit1$residuals))/var(y_st) 
          ##          rsquare = (var(fit1$linear.predictors)                      )/var(y_st)
          ##          fit1 = glm(xs_st , event_st , family = binomial(link = "logit"), control = glmcontrol)  ## , start = beta_val , weights=weights_val_j
          ##          clist =  concordancefit(Surv(  c(rep(1, length(y_st))) , y_st), -fit1$linear.predictors)     ## get concordances  --  allows extension to shose best by C instead of loglik
          ##          loglik = -fit1$deviance/2
          fit1 = lm(y_st ~ riskm)         ###### and in other calls like this ###### , weights=weights_val
          fit1$coefficients[is.na(fit1$coefficients)] = 0                       # assing missing coefficients value 0
          n = length(y_st) 
          kp1 = length(fit1$coefficients) 
          rsquare = var(fit1$fitted.values)/var(y_st) 
          rsquareadj = 1 - ((1-rsquare)*(n-1)/(n-kp1)) 
          rsquare = c(rsquare=rsquare, rsquareadj=rsquareadj)
##          names(fit1)
##          summary(fit1)      
          S2 = var(fit1$residuals)*(n-1)/n
          loglik      = -(n/2)*log(2*pi) - n*log(sqrt(S2)) -(1/(2*S2))*n*S2 
          S2 = var(y_st)*(n-1)/n
          loglik.null = -(n/2)*log(2*pi) - n*log(sqrt(S2)) -(1/(2*S2))*n*S2 
          loglik = c(loglik.null=loglik.null, loglik=loglik)
##          logLik(fit1)
        }  
    
        if (family == "cox") {
          modsum_ = c(ii=i_, jj=j_, best=0, loglik, concord ) 
          modsumt = data.frame(t(modsum_), t(riskvarsnewp),  t(rep(0, xsk)) )     ## first set up data row with place holders for coefficients 
          colnames(modsumt)[(xsk+8):(2*xsk+7)] = xsnames
        } else if (family == "binomial") { 
          modsum_ = c(ii=i_, jj=j_, best=0, loglik, concord ) 
          modsum_
          modsumt = data.frame(t(modsum_), t(riskvarsnewp),  t(rep(0, xsk+1)) )     ## first set up data row with place holders for coefficients 
          colnames(modsumt)[(xsk+8):(2*xsk+8)] = c("Int",xsnames)
        } else if (family == "gaussian") { 
#          modsum_ = c(ii=i_, jj=j_, best=0, loglik.null=loglik.null, loglik=loglik , rsquare=rsquare, rsquareadj=rsquareadj ) 
#          modsum_ = c(ii=i_, jj=j_, best=0, loglik, rsquare=rsquare, rsquareadj=rsquareadj ) 
          modsum_ = c(ii=i_, jj=j_, best=0, loglik, rsquare ) 
          modsum_
          modsumt = data.frame(t(modsum_), t(riskvarsnewp),  t(rep(0, xsk+1)) )     ## first set up data row with place holders for coefficients 
          colnames(modsumt)[(xsk+8):(2*xsk+8)] = c("Int",xsnames)
        }

        ## store model summary for one model fit in this temporary object

        dim(modsumt)                                                 
        modsumt 
        
#        if (rankMatrix(riskm) < 1e-12) { fit1$coefficients = rep(0,dim(riskm)[2]) } 
        fit1$coefficients 
        
        ## save model coefficients in model summary matrix 
#        riskvarsnew = xsnames [riskvarsnewp==1]
        riskvarsnew        
        coefsp = as.numeric( (xsnames %in% riskvarsnew ) )                      ## present(indicators)
        coefsp                             ## same as riskvarsnewp
        coefsi = c(1:xsk)[ (coefsp==1) ]                                        ## numerical index
        coefsi
        if (family == "cox") {
          coefsw = coefsi + xsk + 7                                             ## where (column) in the model summary matrix
          coefsw
          modsumt[ coefsw ] = fit1$coefficients                                 ## save coefficient values to MODel SUMmary 
          modsumt 
        } else {
          coefsw = c(1,coefsi+1) + xsk + 7                                      ## where (column) in the model summary matrix
          coefsw
          modsumt[ coefsw ] = fit1$coefficients                                 ## save coefficient values to MODel SUMmary 
          modsumt 
        }
        ## end i_, j_ coefficient calculation (not loop)  
      } else {
        if (family == "cox") {
          modsumt = data.frame( ii=i_, jj=j_, best=0, modsumlast[ 4:(2*xsk+7) ] )
          colnames(modsumt)[(xsk+8):(2*xsk+7)] = xsnames
        } else {
          modsumt = data.frame( ii=i_, jj=j_, best=0, modsumlast[ 4:(2*xsk+8) ] )
          colnames(modsumt)[(xsk+8):(2*xsk+8)] = c("Int",xsnames)
        }
      } 

      if ((i_==1) & (j_==1)) { 
        modsum  = modsumt                                                       ## store model summaries here in modsum
      } else {
        modsum = rbind( modsum, modsumt )                                       ## store model summaries here in modsum
      }

    }  ## end of j_ loop

    ## find best model from current step (number of model terms), and save for usage in next step 
    modsumi_ = modsum[(modsum$ii==i_),]                                         ## model summaries for step i_ 
#    modsumi_ 
    if (toupper(method)=="LOGLIK") { 
      imax = which.max( modsumi_$loglik )                                       ## find best model log likelihood,  in this iteration by
    } else {
      imax = which.max( modsumi_$concordance )                                  ## find best model concordance,  in this iteration by
    }
    imax
    modsum[ ((i_-1)*xsk + imax), 3] = 1                                         ## set best indicator for best model in entire model summary set  
    modsumlast = modsumi_[ imax , ]                                             ## model summary for best fit model in this LAST iteration  
    modsumlast 
    riskvarsp =  modsumlast[ 8:(xsk+7)]                                         ## predictors in best model (present indicators)
    riskvarsp
    if (family=="cox") { betalast = modsumlast[(xsk+8):(2*xsk+7)]
    } else { betalast = modsumlast[(xsk+8):(2*xsk+8)] }
#    riskvarsi = c(1:xsk)[riskvarsp==1]
#    riskvarsi 

  }  
  ##### end i_ loop ############################################################ 
  if (track>=1) { cat ("\n") }
  
##  colnames(modsum)[6] = "Concordance"
##  colnames(modsum)[7] = "C STD"
  rownames(modsum) = c(1:dim(modsum)[1]) 
#  print(class(modsum))
  class(modsum) <- c("stepreg")
#  print(class(modsum))
  return(modsum)
} 

##############################################################################################################################
##############################################################################################################################

##============================================================================================================================
##============================================================================================================================

#' Get predictors form a stepwise regression model. 
#'
#' @param modsumbest matrix with best predictors based upon number of model terms 
#' @param k_ Value for number of predictors in model 
#' @param risklist Riskset list 
#' @param risklistl Number of terms (length) in the riskset 
#'
#' @return input to best.preds()
#'
#' @noRd

preds_1 = function (modsumbest, k_, risklist, risklistl) {
  temp1 = modsumbest [ , 8:(risklistl+7) ]  
  temp2 = temp1[k_,]  
  temp3 = risklist[ (temp2 ==1) ]
  temp4 =  paste(temp3, collapse = " + " )   
  return(temp4) 
}

##### get lists for all best models of different sizes #########################

#' Get the best models for the steps of a stepreg() fit 
#'
#' @param modsum model summmary
#' @param risklist riskset list 
#'
#' @return best predictors at each step of a stepwise regression
#' 
#' @seealso 
#'    \code{\link{stepreg}} , \code{\link{cv.stepreg}} , \code{\link{nested.glmnetr}}   
#' 
#' @export
#'
best.preds = function(modsum, risklist) {
  risklistl = length( risklist )
  modsumbest = modsum[ (modsum$best==1), ]
#  print ( dim(modsumbest ))
  modsumbestd = dim(modsumbest) 
  k_ = 1 
  preds = data.frame(k_, preds_1(modsumbest, k_, risklist, risklistl) )                                          ## as.data.frame(k_, preds_1(k_) ) 
  for (k_ in 2:modsumbestd[1]) {
    preds_ = data.frame(k_, preds_1(modsumbest, k_, risklist, risklistl) ) 
    preds = rbind (preds, preds_)
  }
  return(preds)
}
##############################################################################################################################
##############################################################################################################################

#' Briefly summarize steps in a stepreg() output object, i.e. a stepwise regression fit
#'
#' @param object A stepreg() output object 
#' @param ... Additional arguments passed to the summary function.  
#'
#' @return Summarize a stepreg() object
#' 
#' @seealso 
#'    \code{\link{stepreg}} , \code{\link{cv.stepreg}} , \code{\link{nested.glmnetr}}
#' 
#' @export

summary.stepreg = function(object, ...) {
  if ( ( dim(object)[2] %% 2) ==0 ) { cox=0 }
  data1 = object[(object$best==1),]
  if (cox==1) { 
    k_ = dim(data1)[2]/2 - 3.5
    keepvars = c(1:(k_+1))[(colSums(abs(data1[,(k_+8):(2*k_+7)]))!=0)] ## + 7 + k_ 
  } else {
    k_ = dim(data1)[2]/2 - 4
    keepvars = c(1:(k_+1))[(colSums(abs(data1[,(k_+8):(2*k_+8)]))!=0)] ## + 7 + k_ 
  }
  keepvars = keepvars + 7 + k_ 
  data2 = data1[,c(1:7,keepvars)]
  return(data2)
}

##############################################################################################################################
##############################################################################################################################
