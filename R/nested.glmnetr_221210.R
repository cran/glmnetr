################################################################################################################################################################
################################################################################################################################################################

##  xs      - predictor input - an n by p matrix, rows for patients (records), columns for predictors
##          - must be n matrix form for complete data, NO NA, NO Inf, etc.  
##  start   - start time, Cox model only - class numeric of length same as number of patients (n)
##  y_      - time, or stop time for Cox model, Y_ 0 or 1 for binomal (logistic), numeric for gaussian 
##          - class numeric of length same as number of patients (n)
##  event   - event indicator, 1 for event, 0 for census, Cox model only
##          - class numeric of length same as number of patients (n)
##  steps_n - number of steps done in stepwise regression fitting 
##  folds_n - number of folds for each level of cross validation 
##  dolasso - fit and do cross validition for lasso model, 0 or 1
##  doaic   - fit and do cross validation for AIC fit, 0 or 1
##  dostep  - fit and do cross validation for stepwise regression fit, 0 or 1
##  method  - method for choosing model in stepwise procedure, "loglik" or "concordance"  
##  family  - model family, "cox", "binomial" or "gaussian" 
##  limit - 1 for only do lambda up to (gamma=1) lambda.min^2/lambda.1se when fitting the relaxed lasso, else (0) lambda.min^3/lambda.1se^2
##  fine  - 1 for smaller differences between lambdas, gives smoother plots, minor relative to variation due to random splits in k-folds

################################################################################################################################################################
################################################################################################################################################################

#' Using nested cross validation, describe the fit of a cross validated tuned relaxed lasso model fit.
#' 
#' @description Performs a nested cross validation for a cross validated 
#' informed relaxed lasso model.  
#' 
#' @param xs     - predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
#' the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
#' and not a data frame. 
#' @param start  - start time, Cox model only - class numeric of length same as number of patients (n)
#' @param y_     - output vector: time, or stop time for Cox model, Y_ 0 or 1 for binomal (logistic), numeric for gaussian. 
#' Must be a vector of length same as number of sample size. 
#' @param event   - event indicator, 1 for event, 0 for census, Cox model only.
#' Must be a numeric vector of length same as sample size.
#' @param steps_n - number of steps done in stepwise regression fitting 
#' @param folds_n - number of folds for each level of cross validation 
#' @param dolasso - fit and do cross validition for lasso model, 0 or 1
#' @param doaic   - fit and do cross validation for AIC fit, 0 or 1.  
#' This is provided for reference only and is not recommended.   
#' @param dostep  - fit and do cross validation for stepwise regression fit, 0 or 1, 
#' as discussed in James, Witten, Hastie and Tibshirani, 2nd edition.  
#' @param method  - method for choosing model in stepwise procedure, "loglik" or "concordance".
#' Other procedures use the "loglik". 
#' @param family  - model family, "cox", "binomial" or "gaussian" 
#' @param lambda  - lambda vector for teh lasso fit
#' @param gamma   - gamma vector for the relaxed lasso fit, default is c(0,0.25,0.5,0.75,1).
#' @param relax   - Fit the relaxed lasso model when fitting a lasso model.  
#' @param limit   - limit the small values for lambda after the initial fit.  
#' This will calcualtions that have minimal impact on the cross validation.  
#' Default is 2 for moderate limitation, 1 for less limitation, 0 for none.   
#' @param fine  - Use a finer step in determining lambda.  Of little value unless one 
#' repeats the cross valiaiton many times to more finely tune the hyper paramters.  
#' See the _glmnet_ documentation.  
#' @param time   - print out the time splits.  
#' @param seed A seed for set.seed() to assure one can get the same results twice.  If NULL 
#' the program will generate a random seed.  Whether specified or NULL, the seed is stored in the output
#' object for future reference.  
#' @param foldid A vector of integers to associate each record to a fold.  Should be integers between 1 and folds_n.
#' These will only be used in the outer folds. 
#'
#' @return- The fit of a cross validated tuned relaxed lasso model fit, obtained by nested cross validation.
#' 
#' @seealso
#'   \code{\link{glmnetr}} , \code{\link{cv.glmnetr}}  , \code{\link{glmnetr.simdata}} , \code{\link{summary.nested.glmnetr}} , , \code{\link{plot.nested.glmnetr}} 
#'   
#' @export
#' 
#' @importFrom stats runif logLik predict cov cor 
#' @importFrom survival Surv coxph coxph.control concordance concordancefit 
#' @importFrom glmnet cv.glmnet 
#' @importFrom Matrix rankMatrix 
#'
#' @examples
#' \donttest{
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$y_ 
#' # for this example we use a small number for folds_n to shorten run time 
#' nested.glmnetr.fit = nested.glmnetr( xs, NULL, y_, NULL, family="gaussian", folds_n=3)
#' plot(nested.glmnetr.fit) 
#' plot(nested.glmnetr.fit, coefs=TRUE) 
#' summary(nested.glmnetr.fit) 
#' summary(nested.glmnetr.fit, cvfit=TRUE) 
#' }
#' 
nested.glmnetr = function(xs, start=NULL, y_, event, family=NULL, steps_n=0, folds_n=10, dolasso=1, doaic=0, dostep=0, method="loglik", lambda=NULL, gamma=NULL, relax=TRUE, limit=1, fine=0, time=0, seed=NULL, foldid=NULL) {
  
  if ( is.null(xs) ) { print(" xs cannot be NULL, program will stop ") ; stop ; }
  if (sum(is.na(xs))> 0) { print(" xs cannot have missing values, program will stop ") ; stop ; }
  if ( is.null(y_) ) { print(" y_ cannot be missing or NULL, program will stop ") ; stop ; }
  if (sum(is.na(y_))> 0) { print(" y_ cannot have missing values, program will stop ") ; stop ; }
#  if ((family=="cox") & (is.null(event))) { print(" event cannot be NULL for Cox model, program will stop ") ; stop ; }
#  if (is.na(start)) { start = NULL }
  if ((family=="cox") & (min(y_) < 0)) { print(" survival times cannot be negative for Cox model, program will stop ") ; stop ; }
#  if ((doaic==1) | (dostep==1)) {
#    if (family!="cox") { doaic=0 ; dostep=0 ; }
#    cat(paste0(" AIC and Stepwise are currenly only implemented for Cox model, and results are only approximate","\n"))
#  }
  if (folds_n < 2) { folds_n = 2 }
  if (time>=1) { cat(paste0("\n", " ##############################################################################################" , "\n")) }
  
  if (time==1) {
    time_start = difftime2() 
    time_last = NULL 
  }
  
  folds_n=max(folds_n, 3)
  
  xs_ncol = dim(xs)[2]                                                          ## number of candidate predictor variables 
  nobs    = dim(xs)[1] 
  one     = c( rep(1, nobs) )
  
  if (steps_n==0) { steps_n=xs_ncol }
  steps_n = min(steps_n, xs_ncol) 
  
  if ( is.null(family)) {family = NA}
  if ( is.na(family) & (!is.null(start)) ) { family = "cox" }
  if ( is.na(family) & (!is.null(event)) & is.null(y_)) { 
    y_ = event 
    event = NULL 
    family = "binomial" 
  } else if ( is.na(family) & (!is.null(event     )) ) { 
    family = "cox" 
  }
  if ( is.na(family) & ( length(table(y_)) == 2 ) ) { family = "binomial" }
  if ( is.na(family) & ( length(table(y_)) >  2 ) ) { family = "gaussian" }
  if ( family == "poison" ) { family = "poisson" }
  if ( family %in% c("gausian", "gause", "normal") ) { family = "gaussian" }
  if ( family %in% c("COX", "cox", "coxph") ) { family = "cox" }
  
  if (is.null(seed)) { seed = round(runif(1)*1000000000) }
  set.seed(seed) 
  if (is.null(foldid)) { foldid = sample( rep( 1:folds_n , ceiling(nobs/folds_n) )[ 1:nobs ] , nobs ) }

  length( foldid )
  table (foldid) 

  ## log likelihoods & conncordances from LASSO and relaxed lasso by Cross Validation 
  lassodeviancv = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 ) 
  lassoagreecv  = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 )
  lassolincalcv = matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 )
  ## number of non-zero coefficients in LASSO models 
  lassonzerocv =  matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6)
  lassogammacv =  matrix ( data=rep(0,(2*folds_n)), nrow = folds_n, ncol = 2)

  ## log likelihoods & conncordances for STEPWISE by Cross Validation 
  ## and coxph models based upon LASSO terms by Cross Validation
  stepdeviancv = matrix ( data=rep(0,(4*folds_n)), nrow = folds_n, ncol = 4 )
  stepagreecv  = matrix ( data=rep(0,(folds_n  )), nrow = folds_n, ncol = 3 )
  step_dfcv    = matrix ( data=rep(0,(folds_n  )), nrow = folds_n, ncol = 3 )
  step_pcv     = matrix ( data=rep(0,(folds_n  )), nrow = folds_n, ncol = 1 )

  foldid_ = NULL 
  if (dolasso == 1) {
    if (time>=1) { cat(paste0(" Initial lasso fit of all data ", "\n")) }
########     (xs_tmp, start_tmp, y_tmp, event_tmp, lambda=NULL, gamma=c(0,0.25,0.50,0.75,1), object=NULL, time=0, family="cox" ) ; time = 1 ; 
    if (relax) { cv.glmnet.fit.f = cv.glmnetr( xs, start, y_, event, lambda=lambda, gamma=gamma, folds_n=folds_n, limit=limit, fine=fine, time=time, family=family) 
    } else {
      if (family=="cox") {
        if ( is.null(start) ) { cv.glmnet.fit.f = cv.glmnet( xs, Surv(y_, event)       , family="cox", lambda=lambda, relax=FALSE ) 
        } else                { cv.glmnet.fit.f = cv.glmnet( xs, Surv(start, y_, event), family="cox", lambda=lambda, relax=FALSE ) } 
      } else {
       cv.glmnet.fit.f = cv.glmnet( xs, y_, family=family, lambda=lambda, relax=FALSE) 
      }
    }

    if (time >= 1) { cat(paste0(" length(lambda) = " , length(cv.glmnet.fit.f$lambda), "\n" )) } 
    
    if (time==1) {
      time_last = difftime2(time_start, time_last)
    }
  } 
  
  if (dostep == 1) { 
    if (time >= 1){ cat(paste0("\n", " ########## Derive stepwise model on all data #################################################" , "\n")) }
    cv.stepreg.fit.all  = cv.stepreg(xs, start, y_, event, steps_n, folds_n, method=method, family=family,foldid=foldid_,time=time)
    if (time==1) { time_last = difftime2(time_start, time_last) }
  }
    
  if (doaic ==1) { 
    if (dostep==1) { 
      func.fit.aic = aicreg(xs, start, y_, event, family=family, object=cv.stepreg.fit.all, time=time) 
    } else {
      func.fit.aic = aicreg(xs, start, y_, event, steps_n=steps_n, family=family, object=NULL, time=time) 
    }
    if (time==1) { time_last = difftime2(time_start, time_last) } 
  }
  
  ##### FOLDS ##############################################################################################
  ##### FOLDS ##############################################################################################
  ## i_ = 1 
  for (i_ in 1:folds_n) {
    ##=== begin folds ==========================================================
    if (time>=1) { cat(paste0("\n", " ########## Entering Nested Cross Validation outer fold  ", i_, "  of  " , folds_n , "  ############################" , "\n")) }
  
    ##### set up train and test data sets in matrix form for glmnet & stepreg #####
    trainxs = xs[(foldid!=i_),]   
    testxs  = xs[(foldid==i_),]  
    test1   = c(rep(1,dim(testxs)[1])) 
    testxs1 = cbind(test1, testxs)                                              ## for models with Int ercept 
    
    dim(trainxs)[1] + dim(testxs)[1] 
    dim(trainxs) 
    dim(testxs) 
    
    if ( is.null(start) ) {
      trainstart = NULL
      teststart  = NULL
    } else {
      trainstart = start[(foldid!=i_)]    
      teststart  = start[(foldid==i_)]   
    }
  
    trainy_ = y_[(foldid!=i_)]    
    testy_  = y_[(foldid==i_)]   
    
    if (family == "cox") {
      trainevent = event[(foldid!=i_)] ;  
      testevent  = event[(foldid==i_)] ; 
    } else {
      trainevent = NULL
      testevent  = NULL 
    }    
    
    if (family=="cox") {
      if ( is.null(start) ) {
        train_data = data.frame(trainy_, trainevent, trainxs)
        test_data  = data.frame(testy_ , testevent , testxs )
      } else {
        train_data = data.frame(trainstart, trainy_, trainevent, trainxs)
        test_data  = data.frame(teststart , testy_ , testevent , testxs )
      }
    } else {
      train_data = data.frame(trainy_, trainxs)
      test_data  = data.frame(testy_ , testxs )
    }

    foldid_ = NULL
    ##### LASSO fit ##########################################################################################
    if (dolasso==1) { 
      if (time>=1) { cat( "Starting LASSO model fits ", "\n" ) } 

      if (relax) { 
     #  xs=trainxs; start=trainstart; y_=trainy_; event=trainevent; lambda=lambda; gamma=gamma; folds_n=folds_n; limit=limit; fine=fine; time=1; family=family ;         
        cv.glmnet.fit = cv.glmnetr( trainxs, trainstart, trainy_, trainevent, lambda=lambda, gamma=gamma, folds_n=folds_n, limit=limit, fine=fine, time=0, family=family) 
      } else {
        if (family=="cox") {   
          if ( is.null(start) ) { cv.glmnet.fit = cv.glmnet( trainxs, Surv(trainy_, trainevent)            , family="cox", lambda=lambda, relax=FALSE ) 
          } else                { cv.glmnet.fit = cv.glmnet( trainxs, Surv(trainstart, trainy_, trainevent), family="cox", lambda=lambda, relax=FALSE ) }
        } else {
          cv.glmnet.fit = cv.glmnet( trainxs, trainy_, family=family, lambda=lambda, relax=FALSE) 
        }
      }
      foldid_ = cv.glmnet.fit$foldid 

# summary(cv.glmnet.fit)      
# cv.glmnet.fit$relaxed 
# names(cv.glmnet.fit$relaxed)
      
      lassonzerocv[i_,1] = cv.glmnet.fit$nzero [ cv.glmnet.fit$index[2] ]
      lassonzerocv[i_,2] = cv.glmnet.fit$nzero [ cv.glmnet.fit$index[1] ]     
      if (relax) {
        lassonzerocv[i_,3] = cv.glmnet.fit$relaxed$nzero.1se 
        lassonzerocv[i_,4] = cv.glmnet.fit$relaxed$nzero.min 
        lassonzerocv[i_,5] = cv.glmnet.fit$nzero [ cv.glmnet.fit$relaxed$index.g0[2] ]
        lassonzerocv[i_,6] = cv.glmnet.fit$nzero [ cv.glmnet.fit$relaxed$index.g0[1] ]
      }
      lassogammacv[i_,1]= cv.glmnet.fit$relaxed$gamma.1se 
      lassogammacv[i_,2]= cv.glmnet.fit$relaxed$gamma.min 
      
##X     print( "HERE nested 2" )

      #### LASSO predictions on TEST data ######################################################################
      # object=cv.glmnet.fit ; xs_new=testxs ; lam="lambda.min" ; gam="gamma.min" ;  comment=TRUE
      pred1se   = predict(cv.glmnet.fit, testxs, lam=cv.glmnet.fit$lambda.1se, gam=1)  ##### 1se lambda model predictions ################
      predmin   = predict(cv.glmnet.fit, testxs, lam=cv.glmnet.fit$lambda.min, gam=1)  ##### minimizing lambda model predictions #########
      pred1seR  = predict(cv.glmnet.fit, testxs, lam="lambda.1se" , gam="gamma.1se" )  ###### 1se RELAXED lasso ##########################
      predminR  = predict(cv.glmnet.fit, testxs, lam="lambda.min" , gam="gamma.min" )  ###### min RELAXED lasso ##########################
      pred1seR0 = predict(cv.glmnet.fit, testxs, lam=cv.glmnet.fit$relaxed$lambda.1se.g0, gam=0)  ###### 1se gamma = 0 RELAXED lasso #####
      predminR0 = predict(cv.glmnet.fit, testxs, lam=cv.glmnet.fit$relaxed$lambda.min.g0, gam=0)  ###### min gamma = 0 RELAXED lasso #####
      
#      object = cv.glmnet.fit ; xs_new=testxs ; lam=cv.glmnet.fit$lambda.1se ; gam=1 ; comment = TRUE 
      
#      print( cor(cbind(pred1se, predmin, pred1seR, predminR, pred1seR0, predminR0)) ) 
#      print( cor(cbind(predmin, predminR, predminR0)) ) 
      
      if (family=="cox") {    
        if ( is.null(start) ) { SURV = Surv(testy_, testevent) 
        } else {                SURV = Surv(teststart, testy_, testevent) }
        fit1 = coxph( SURV ~ pred1se  , init=c(1))  #, control=coxph.control(iter.max=0) )  
        fit2 = coxph( SURV ~ predmin  , init=c(1))  #, control=coxph.control(iter.max=0) ) 
        fit3 = coxph( SURV ~ pred1seR , init=c(1))  #, control=coxph.control(iter.max=0) ) 
        fit4 = coxph( SURV ~ predminR , init=c(1))  #, control=coxph.control(iter.max=0) ) 
        fit5 = coxph( SURV ~ pred1seR0, init=c(1))  #, control=coxph.control(iter.max=0) ) 
        fit6 = coxph( SURV ~ predminR0, init=c(1))  #, control=coxph.control(iter.max=0) ) 
        lassodeviancv[i_,1] = -2*fit1$loglik[1] / fit1$nevent  
        lassodeviancv[i_,2] = -2*fit2$loglik[1] / fit2$nevent  
        lassodeviancv[i_,3] = -2*fit3$loglik[1] / fit3$nevent 
        lassodeviancv[i_,4] = -2*fit4$loglik[1] / fit4$nevent  
        lassodeviancv[i_,5] = -2*fit5$loglik[1] / fit5$nevent   
        lassodeviancv[i_,6] = -2*fit6$loglik[1] / fit6$nevent 
        lassolincalcv[i_,1] = fit1$coefficients
        lassolincalcv[i_,2] = fit2$coefficients
        lassolincalcv[i_,3] = fit3$coefficients
        lassolincalcv[i_,4] = fit4$coefficients
        lassolincalcv[i_,5] = fit5$coefficients
        lassolincalcv[i_,6] = fit6$coefficients
      } else if (family=="binomial") {
#        fit1 = glm( testevent ~ pred1se   , start=c(0,1), family="binomial")  #, control=coxph.control(iter.max=0) )  
        fit1 = glm( testy_ ~ pred1se   , start=c(0,1), family=family)  #, control=coxph.control(iter.max=0) )  
        fit1 = glm( testy_ ~ pred1se   , family=family)  #, control=coxph.control(iter.max=0) )  
        fit2 = glm( testy_ ~ predmin   , family=family)  #, control=coxph.control(iter.max=0) ) 
        fit3 = glm( testy_ ~ pred1seR  , family=family)  #, control=coxph.control(iter.max=0) ) 
        fit4 = glm( testy_ ~ predminR  , family=family)  #, control=coxph.control(iter.max=0) ) 
        fit5 = glm( testy_ ~ pred1seR0 , family=family)  #, control=coxph.control(iter.max=0) ) 
        fit6 = glm( testy_ ~ predminR0 , family=family)  #, control=coxph.control(iter.max=0) ) 
        p_ = 1/(1+exp(-pred1se  )) ;  lassodeviancv[i_,1] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_) ; 
        p_ = 1/(1+exp(-predmin  )) ;  lassodeviancv[i_,2] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_) ; 
        p_ = 1/(1+exp(-pred1seR )) ;  lassodeviancv[i_,3] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_) ; 
        p_ = 1/(1+exp(-predminR )) ;  lassodeviancv[i_,4] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_) ; 
        p_ = 1/(1+exp(-pred1seR0)) ;  lassodeviancv[i_,5] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_) ; 
        p_ = 1/(1+exp(-predminR0)) ;  lassodeviancv[i_,6] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_) ; 
        lassolincalcv[i_,1] = fit1$coefficients[2]
        lassolincalcv[i_,2] = fit2$coefficients[2]
        lassolincalcv[i_,3] = fit3$coefficients[2]
        lassolincalcv[i_,4] = fit4$coefficients[2]
        lassolincalcv[i_,5] = fit5$coefficients[2]
        lassolincalcv[i_,6] = fit6$coefficients[2]
      } else if (family=="gaussian")  {
        lassodeviancv[i_,1] = sum((testy_-pred1se  )^2) / length(testy_) 
        lassodeviancv[i_,2] = sum((testy_-predmin  )^2) / length(testy_) 
        lassodeviancv[i_,3] = sum((testy_-pred1seR )^2) / length(testy_) 
        lassodeviancv[i_,4] = sum((testy_-predminR )^2) / length(testy_) 
        lassodeviancv[i_,5] = sum((testy_-pred1seR0)^2) / length(testy_) 
        lassodeviancv[i_,6] = sum((testy_-predminR0)^2) / length(testy_) 
        lassolincalcv[i_,1] = cov(testy_,pred1se  )/var(pred1se  )
        lassolincalcv[i_,2] = cov(testy_,predmin  )/var(predmin  )
        lassolincalcv[i_,3] = cov(testy_,pred1seR  )/var(pred1seR )
        lassolincalcv[i_,4] = cov(testy_,predminR )/var(predminR )
        lassolincalcv[i_,5] = cov(testy_,pred1seR0)/var(pred1seR0)
        lassolincalcv[i_,6] = cov(testy_,predminR0)/var(predminR0)
      }
      
      if (family == "cox") {
        lassoagreecv[i_,1] = fit1$concordance[6] 
        lassoagreecv[i_,2] = fit2$concordance[6]  
        lassoagreecv[i_,3] = fit3$concordance[6] 
        lassoagreecv[i_,4] = fit4$concordance[6]  
        lassoagreecv[i_,5] = fit5$concordance[6] 
        lassoagreecv[i_,6] = fit6$concordance[6] 
      } else if (family == "binomial") { 
        lassoagreecv[i_,1] = concordance(testy_ ~ pred1se  )[[1]] 
        lassoagreecv[i_,2] = concordance(testy_ ~ predmin  )[[1]]   
        lassoagreecv[i_,3] = concordance(testy_ ~ pred1seR )[[1]]   
        lassoagreecv[i_,4] = concordance(testy_ ~ predminR )[[1]]   
        lassoagreecv[i_,5] = concordance(testy_ ~ pred1seR0)[[1]]   
        lassoagreecv[i_,6] = concordance(testy_ ~ predminR0)[[1]]           
#        lassoagreecv[i_,1] = with(summary(fit1), 1-deviance/null.deviance)         
      } else if (family == "gaussian") { 
        if (var(pred1se)>0)   { lassoagreecv[i_,1] =  cor(x=testy_, y=pred1se  )^2 }
        if (var(predmin)>0)   { lassoagreecv[i_,2] =  cor(x=testy_, y=predmin  )^2 }
        if (var(pred1seR)>0)  { lassoagreecv[i_,3] =  cor(x=testy_, y=pred1seR )^2 }
        if (var(predminR)>0)  { lassoagreecv[i_,4] =  cor(x=testy_, y=predminR )^2 }
        if (var(pred1seR0)>0) { lassoagreecv[i_,5] =  cor(x=testy_, y=pred1seR0)^2 }
        if (var(predminR0)>0) { lassoagreecv[i_,6] =  cor(x=testy_, y=predminR0)^2 }
      }
  
      ##### print some summary results ##################################################
##X      cat(paste0("\n", " Number of nonzero model terms for LASSO (1se and min): ", lassonzerocv[i_,1] , ", ", lassonzerocv[i_,2], "\n"))
##X      cat(paste0(      " Gamma term for maximizing relaxed lass model (1se and min): ", lassogammacv[i_,1] , ", ", lassogammacv[i_,2], "\n"))
##X      cat(paste0("\n", " LASSO loglik, 1se and min and cox models based upon temrs in 1se and min", "\n"))
##X      print( round( lassodeviancv[i_,] , digits=4) )
##X      cat(paste0(      " Concordance, 1se and min and cox models baed upon terms in 1se and min", "\n"))
##X      print( round( lassoagreecv[i_,] , digits=4) ) 
##X      cat("\n")
      
##    cov(cbind(pred1se, predmin, pred1seR, predminR, pred1seR0, predminR0))
##    lassogammacv
    }
    
    ###### START STEPWISE fits #################################################
    if (dostep == 1) { 
      if (time >= 1) { cat(paste0("\n", " Stepwise regression calculations ... ", "\n")) }
      ## xs_cv=trainxs ; start_cv=trainstart ; y_cv=trainy_ ; event_cv=trainevent ; steps_n=steps_n ; folds_n_cv=folds_n ; method=method ; family=family ; 
      cv.stepreg.fit = cv.stepreg(trainxs, trainstart, trainy_, trainevent, steps_n, folds_n, method=method, family=family, foldid=foldid_, time=time) 
      stepreg.fit.all.best = cv.stepreg.fit$stepreg.fit.all.best 
#      names(cv.stepreg.fit) ; cv.stepreg.fit$best.p ; cv.stepreg.fit$df.p ; cv.stepreg.fit$pvalue
#            names(cv.stepreg.fit) ; cv.stepreg.fit$best.p ; cv.stepreg.fit$func.fit.p      
#      str(cv.stepreg.fit)
#      cat(paste0( cv.stepreg.fit$stepreg.fit.all.best[1,1] , cv.stepreg.fit$stepreg.fit.all.best[2,2], cv.stepreg.fit$stepreg.fit.all.best[3,3], cv.stepreg.fit$stepreg.fit.all.best[4,4])) 

#     stepreg.fit.all.best[1:10, 1:7] 
      #### stepwise tuned by DF ################################################
#      df = cv.stepreg.fit$best.df
#      stepdeviancv [i_,4] = stepreg.fit.all.best [ df , 4 ]        
#      stepdeviancv [i_,1] = stepreg.fit.all.best [ df , 5 ]  
#      stepagreecv  [i_,1] = stepreg.fit.all.best [ df , 6 ] 
#      step_dfcv    [i_,1] = df 

      func.fit.df = cv.stepreg.fit$func.fit.df                                  ##  cox model suggested by Cross Validation 
      mod_df      = length(func.fit.df$coefficients)                            ##  Number of terms in Cross Validation model 
      testscore  = predict( func.fit.df, as.data.frame( testxs) )              ##  predicted scores for test data based upon train data model 
      if (family == "cox") {
        if (is.null(start)) { 
          fit0 = coxph( Surv(testy_, testevent) ~ rep(1,length(testscore)), init=c(1), control=coxph.control(iter.max=0))
          fit3 = coxph( Surv(testy_, testevent) ~ testscore, init=c(1), control=coxph.control(iter.max=0))
        } else { 
          fit0 = coxph( Surv(testy_, teststart, testevent) ~ rep(1,length(testscore)), init=c(1), control=coxph.control(iter.max=0))
          fit3 = coxph( Surv(testy_, teststart, testevent) ~ testscore, init=c(1), control=coxph.control(iter.max=0))  
        }
        stepdeviancv[i_,4] = -2*fit0$loglik[1] / fit0$nevent  
        stepdeviancv[i_,1] = -2*fit3$loglik[1] / fit3$nevent  
        stepagreecv [i_,1] = fit3$concordance[6]   
        step_dfcv   [i_,1] = mod_df 
      } else if (family == "binomial") {      
        p_ = sum(testy_)/length(testy_) 
        stepdeviancv[i_,4] = -2 * sum( log(p_)*testy_ + log(1-p_)*(1-testy_) ) / length(testy_)
        p_ = 1/(1+exp(-testscore)) 
        stepdeviancv[i_,1] = -2 * ( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_)
        stepagreecv [i_,1] = concordance(testy_ ~ testscore)[[1]]
        step_dfcv   [i_,1] = mod_df - 1 
      }  else if (family == "gaussian") {  
        stepdeviancv[i_,4] = (var(testy_) * (length(testy_)-1)) / length(testy_) 
        stepdeviancv[i_,1] = sum((testy_-testscore)^2)  / length(testscore) 
        stepagreecv [i_,1] = cor(x=testy_,y=testscore)^2         
        step_dfcv   [i_,1] = mod_df - 1 
      }
      
      #### stepwise tuned by p (critical value) ################################
      func.fit.p = cv.stepreg.fit$func.fit.p                                  ##  cox model suggested by Cross Validation 
      mod_df     = length(func.fit.p$coefficients)                            ##  Number of terms in Cross Validation model 
      testscore  = predict( func.fit.p, as.data.frame( testxs) )              ##  predicted scores for test data based upon train data model 
      if (family == "cox") {
        if (is.null(start)) { 
          fit3 = coxph( Surv(testy_, testevent) ~ testscore, init=c(1), control=coxph.control(iter.max=0))
        } else { 
          fit3 = coxph( Surv(testy_, teststart, testevent) ~ testscore, init=c(1), control=coxph.control(iter.max=0))  
        }
        stepdeviancv[i_,2] = -2*fit3$loglik[1] / fit3$nevent  
        stepagreecv [i_,2] = fit3$concordance[6]   
        step_dfcv   [i_,2] = mod_df 
      } else if (family == "binomial") {      
        p_ = 1/(1+exp(-testscore)) 
        stepdeviancv[i_,2] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testy_)
        stepagreecv [i_,2] = concordance(testy_ ~ testscore)[[1]]
        step_dfcv   [i_,2] = mod_df - 1 
      }  else if (family == "gaussian") {  
        stepdeviancv[i_,2] = sum((testy_-testscore)^2) / length(testscore) 
        stepagreecv [i_,2] = cor(x=testy_,y=testscore)^2         
        step_dfcv   [i_,2] = mod_df - 1 
      }
      step_pcv    [i_,1] = cv.stepreg.fit$best.p 
      
      if (time==1) { time_last = difftime2(time_start, time_last) }
    }
    
    ## BEGIN AIC FIT (fold) ####################################################
    if (doaic==1){
      if (dostep==1) {
        stepreg.fit.best = as.matrix( cv.stepreg.fit$stepreg.fit.all.best ) 
#        class(stepreg.fit) 
      } else {
        stepreg.fit = stepreg(trainxs, trainstart, trainy_, trainevent, steps_n=steps_n)  
#        class(stepreg.fit) 
        class(stepreg.fit) = "data.frame" 
        stepreg.fit.best = stepreg.fit[(stepreg.fit$best==1) , ] 
#        class( stepreg.fit.best )
      }
#      stepreg.fit.best[1:10,1:7]
      aic = 2*stepreg.fit.best[,1] - 2*stepreg.fit.best[,5]                     
      mod_df  = which.min( aic )                                                    ## may be differnt from first p > 0.15 XX
   
      if (family=="cox") {
        beta = stepreg.fit.best [ mod_df , (xs_ncol+8):(2*xs_ncol+7) ]              ## get beta, the regression estimate
#        class(beta)
#        beta = matrix(beta, nrow=xs_ncol, ncol=1) 
      } else {
        beta = stepreg.fit.best [ mod_df , (xs_ncol+8):(2*xs_ncol+8) ]              ## get beta, the regression estimate
#        class(beta)
#        beta = matrix(beta, nrow=(xs_ncol+1), ncol=1) 
      }
      
      beta = as.numeric(beta)
      if (family=="cox") { testscore  = testxs  %*% beta
      } else             { testscore  = testxs1 %*% beta }                      ## get predicteds 
      
      if (family=="cox") {
        if (is.null(start)) { SURV = Surv(           testy_, testevent)
        } else {              SURV = Surv(teststart, testy_, testevent) }
        fit0  = coxph( SURV ~ rep(1,length(testscore)), control=coxph.control(iter.max=0) )  
        fit1  = coxph( SURV ~ testscore, init=c(1), control=coxph.control(iter.max=0) )  
        stepdeviancv[i_,4] = -2*fit0$loglik[1] / fit0$nevent  
        stepdeviancv[i_,3] = -2*fit1$loglik[1] / fit1$nevent
        stepagreecv [i_,3] = fit1$concordance[6]   
        step_dfcv   [i_,3] = mod_df
      } else if (family=="binomial") {
        if (dostep!=1) {
          p_ = sum(testy_)/length(testy_) 
          stepdeviancv[i_,4] = -2*sum( log(p_)*testy_ + log(1-p_)*(1-testy_) )
        }
        p_ = 1/(1+exp(-testscore)) 
        stepdeviancv[i_,3] = -2*( t(log(p_))%*%testy_ + t(log(1-p_))%*%(1-testy_) ) / length(testscore)  
        stepagreecv [i_,3] = concordance(testy_ ~ testscore)[[1]]   
        step_dfcv   [i_,3] = mod_df - 1 
      } else if (family=="gaussian") {
        if (dostep!=1) {
          stepdeviancv[i_,4] = (var(testy_) * (length(testy_)-1)) / length(testy_)
        }
        stepdeviancv[i_,3] = sum((testy_-testscore)^2)  / length(testscore) 
        stepagreecv [i_,3] = cor(x=testy_,y=testscore)^2 
        step_dfcv   [i_,3] = mod_df - 1 
      } 

      if ((time==1) & (dostep==0)) { time_last = difftime2(time_start, time_last) }
    }
    ## end within fold AIC FIT #################################################
  }
  ##### END FOLDS ##########################################################################################
  ##### END FOLDS ##########################################################################################

  if (time >= 1) { cat(paste0("\n", " ################################################################################################" , "\n")) } 
  
  if (dolasso==1) { 
    if (time >= 1) { cat(paste0("\n", " ########## calculate lasso fit on all data ########################################################" , "\n")) }

    cv.glmnet.fit = cv.glmnet.fit.f 
    
#    lassonzerocv[i_,1]= cv.glmnet.fit$relaxed$nzero.1se 
#    lassonzerocv[i_,2]= cv.glmnet.fit$relaxed$nzero.min 
#    lassogammacv[i_,1]= cv.glmnet.fit$relaxed$gamma.1se 
#    lassogammacv[i_,2]= cv.glmnet.fit$relaxed$gamma.min 

    pred1se   = predict(cv.glmnet.fit, xs, lam=cv.glmnet.fit$lambda.1se, gam=1)
    predmin   = predict(cv.glmnet.fit, xs, lam=cv.glmnet.fit$lambda.min, gam=1)
    pred1seR  = predict(cv.glmnet.fit, xs, lam="lambda.1se" , gam="gamma.1se" )
    predminR  = predict(cv.glmnet.fit, xs, lam="lambda.min" , gam="gamma.min" )
    pred1seR0 = predict(cv.glmnet.fit, xs, lam=cv.glmnet.fit$lambda.1se, gam=0)
    predminR0 = predict(cv.glmnet.fit, xs, lam=cv.glmnet.fit$lambda.min, gam=0)
    
    if (family=="cox") {    
      if ( is.null(start) ) { 
        fit1 = coxph( Surv(y_, event) ~ pred1se   ) 
        fit2 = coxph( Surv(y_, event) ~ predmin   ) 
        fit3 = coxph( Surv(y_, event) ~ pred1seR  ) 
        fit4 = coxph( Surv(y_, event) ~ predminR  ) 
        fit5 = coxph( Surv(y_, event) ~ pred1seR0 ) 
        fit6 = coxph( Surv(y_, event) ~ predminR0 ) 
      } else {     
        fit1 = coxph( Surv(start, y_, event) ~ pred1se   ) 
        fit2 = coxph( Surv(start, y_, event) ~ predmin   ) 
        fit3 = coxph( Surv(start, y_, event) ~ pred1seR  ) 
        fit4 = coxph( Surv(start, y_, event) ~ predminR  ) 
        fit5 = coxph( Surv(start, y_, event) ~ pred1seR0 ) 
        fit6 = coxph( Surv(start, y_, event) ~ predminR0 ) 
        }
    } else  {      
      fit1 = glm( y_ ~ pred1se   , family=family)               
      fit2 = glm( y_ ~ predmin   , family=family)               
      fit3 = glm( y_ ~ pred1seR  , family=family)               
      fit4 = glm( y_ ~ predminR  , family=family)               
      fit5 = glm( y_ ~ pred1seR0 , family=family)               
      fit6 = glm( y_ ~ predminR0 , family=family)               
    }
    
    if (family %in% c("cox")) {
      lasso.naive.c.1se   = fit1$concordance[6] 
      lasso.naive.c.min   = fit2$concordance[6]  
      lasso.naive.c.1seR  = fit3$concordance[6] 
      lasso.naive.c.minR  = fit4$concordance[6]  
      lasso.naive.c.1seR0 = fit5$concordance[6] 
      lasso.naive.c.minR0 = fit6$concordance[6] 
    } else if (family %in% c("binomial")) {
      lasso.naive.c.1se   = concordance(fit1)[[1]] 
      lasso.naive.c.min   = concordance(fit2)[[1]]
      lasso.naive.c.1seR  = concordance(fit3)[[1]]
      lasso.naive.c.minR  = concordance(fit4)[[1]]
      lasso.naive.c.1seR0 = concordance(fit5)[[1]]
      lasso.naive.c.minR0 = concordance(fit6)[[1]]
    } else if (family %in% c("gaussian")) {
      lasso.naive.c.1se   = with(summary(fit1), 1-deviance/null.deviance)  
      lasso.naive.c.min   = with(summary(fit2), 1-deviance/null.deviance)
      lasso.naive.c.1seR  = with(summary(fit3), 1-deviance/null.deviance)
      lasso.naive.c.minR  = with(summary(fit4), 1-deviance/null.deviance)
      lasso.naive.c.1seR0 = with(summary(fit5), 1-deviance/null.deviance)
      lasso.naive.c.minR0 = with(summary(fit6), 1-deviance/null.deviance)
    } 
    lasso.naive.agree = c(naive.c.1se  =lasso.naive.c.1se  , naive.c.min  =lasso.naive.c.min , 
                          naive.c.1seR =lasso.naive.c.1seR , naive.c.minR =lasso.naive.c.minR, 
                          naive.c.1seR0=lasso.naive.c.1seR0, naive.c.minR0=lasso.naive.c.minR0 )
    names(lasso.naive.agree) = c("naive.c.1se", "naive.c.min", "naive.c.1seR", 
                                 "naive.c.minR", "naive.c.1seR0", "naive.c.minR0" ) 
  }

  if (time >= 1) { cat(paste0("\n", " ##############################################################################################" , "\n")) }

#  colnames(lassonzerocv) = c("LASSO.1se", "LASSO.min")
  colnames(lassodeviancv) = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0" )
  colnames(lassoagreecv ) = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0" )
  colnames(lassolincalcv) = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0" )
  colnames(lassonzerocv)  = c("lasso.1se", "lasso.min", "lasso.1seR", "lasso.minR", "lasso.1seR0", "lasso.minR0" )
  colnames(stepdeviancv) = c("df", "p", "AIC", "null")
  colnames(stepagreecv ) = c("df", "p", "AIC")
  colnames( step_dfcv  ) = c("df", "p", "AIC")
  colnames( step_pcv   ) = c("p stepwise")

#  row.names(c.naive) = NULL
  if        (family=="cox"     ) { nevents=sum(event)
  } else if (family=="binomial") { nevents=sum(y_) 
  } else                         { nevents=NA         } 
  nestedcv = list( sample=c(family=family, n=dim(xs)[1], nevents=nevents, xs.columns= dim(xs)[2], xs.df=rankMatrix(xs)[1]), 
                   tuning=c(steps_n=steps_n,folds_n=folds_n,method=method,dolasso=dolasso,doaic=doaic,dostep=dostep,seed=seed),
                   foldid=foldid ) 

  if (dolasso == 1) { nestedcv = c(nestedcv, lassodeviancv=list(lassodeviancv), lassoagreecv=list(lassoagreecv),
                                   lassolincalcv=list(lassolincalcv),
                                   lassonzerocv=list(lassonzerocv) , lasso.naive.agree=list(lasso.naive.agree),
                                   cv.glmnet.fit=list(cv.glmnet.fit) , cv.glmnet.fit.f=list(cv.glmnet.fit.f) ) } 

  if ((dostep == 1) | (doaic==1)) { nestedcv = c(nestedcv, stepdeviancv =list(stepdeviancv) , 
                                                 stepagreecv =list(stepagreecv), step_dfcv  =list(step_dfcv), step_pcv=list(step_pcv) ) }
  
  if (dostep  == 1) { nestedcv = c(nestedcv, cv.stepreg.fit=list(cv.stepreg.fit.all) ) }
  
  if (doaic   == 1) { nestedcv = c(nestedcv, func.fit.aic=list(func.fit.aic) ) } 

  class(nestedcv) <- c("nested.glmnetr")

  names (nestedcv)

  return( nestedcv )

}

################################################################################################################################################################
################################################################################################################################################################
