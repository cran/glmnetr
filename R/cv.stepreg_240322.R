################################################################################
##### cv.stepreg_yymmdd.R ######################################################
################################################################################
#' Cross validation informed stepwise regression model fit.   
#'
#' @param xs_cv   predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
#' the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
#' and not a data frame. 
#' @param start_cv  start time, Cox model only - class numeric of length same as number of patients (n)
#' @param y_cv      output vector: time, or stop time for Cox model, Y_ 0 or 1 for binomal (logistic), numeric for gaussian. 
#' #' Must be a vector of length same as number of sample size. 
#' @param event_cv  event indicator, 1 for event, 0 for census, Cox model only.
#' Must be a numeric vector of length same as sample size.
#' @param steps_n  Maximun number of steps done in stepwise regression fitting.  If 0, then takes the value rank(xs_cv).   
#' @param folds_n  number of folds for cross validation 
#' @param method   method for choosing model in stepwise procedure, "loglik" or "concordance".
#' Other procedures use the "loglik".     
#' @param family   model family, "cox", "binomial" or "gaussian" 
#' @param seed a seed for set.seed() to assure one can get the same results twice.  If NULL 
#' the program will generate a random seed.  Whether specified or NULL, the seed is stored in the output
#' object for future reference.  
#' @param foldid a vector of integers to associate each record to a fold.  The integers should be between 1 and folds_n.
#' @param stratified folds are to be constructed stratified on an indicator outcome 
#' 1 (default) for yes, 0 for no.  Pertains to event variable for "cox" and y_ for 
#' "binomial" family. 
#' @param track indicate whether or not to update progress in the console.  Default of
#' 0 suppresses these updates.  The option of 1 provides these updates.  
#' In fitting clinical data with non full rank design matrix we have found some 
#' R-packages to take a very long time.  Therefore we allow the user to track the 
#' program progress and judge whether things are 
#' moving forward or if the process should be stopped. 
#'
#' @return cross validation infomred stepwise regression model fit tuned by number of model terms or p-value for inclusion.
#' 
#' @seealso 
#'    \code{\link{predict.cv.stepreg}} , \code{\link{summary.cv.stepreg}}, \code{\link{stepreg}} , \code{\link{aicreg}} , \code{\link{nested.glmnetr}}
#' 
#' @export
#' 
#' @importFrom stats formula lm pchisq runif var glm glm.control   
#' @importFrom survival concordance coxph 
#'
#' @examples 
#' set.seed(955702213)
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=c(0,1,1))
#' # this gives a more interesting case but takes longer to run
#' xs=sim.data$xs           
#' # this will work numerically as an example 
#' xs=sim.data$xs[,c(2,3,50:55)] 
#' dim(xs)
#' y_=sim.data$yt 
#' event=sim.data$event
#' # for this example we use small numbers for steps_n and folds_n to shorten run time 
#' cv.stepreg.fit = cv.stepreg(xs, NULL, y_, event, steps_n=10, folds_n=3, track=0)
#' summary(cv.stepreg.fit)
#' 
cv.stepreg = function(xs_cv, start_cv=NULL, y_cv, event_cv,  family="cox", 
                      steps_n=0, folds_n=10, method="loglik",
                      seed=NULL, foldid=NULL, stratified=1, track=0) {

  if ( is.logical(y_cv) & (family == "binomial")) { 
    cat("\n  The y_cv variable is of class logical and is converted to numeric\n")
    y_cv = y_cv * 1 
  }
  
  xs_cv_ncol = dim(xs_cv)[2]    
  nobs_cv = dim(xs_cv)[1]   
  
  folds_n=max(folds_n, 3)
  
  if (!is.null(foldid)) {
    folds_index_cv = foldid 
  } else {
#    folds_index_cv = sample( rep( 1:folds_n , ceiling(nobs_cv/folds_n) )[ 1:nobs_cv ] , nobs_cv ) 
    if  (is.null(seed))   { seed = round(runif(1)*1e9) 
    } else if (seed == 0) { seed = round(runif(1)*1e9) 
    }
    set.seed(seed) 
    folds_index_cv = get.foldid(y_cv, event_cv, family, folds_n, stratified) 
  }

  if (steps_n==0) { steps_n=xs_cv_ncol } 
  steps_n = min(steps_n, xs_cv_ncol)
  
  ## log likelihoods & conncordances for STEPWISE by cross validation 
  ## and coxph models based upon LASSO terms by cross validation
  steploglikcv.df = matrix ( data=rep(0,(folds_n*(steps_n+1))), nrow = folds_n, ncol = (steps_n+1))
  stepagreecv.df  = matrix ( data=rep(0,(folds_n*(steps_n  ))), nrow = folds_n, ncol = (steps_n  ))
  step_ncv        = matrix ( data=rep(0,(folds_n*(steps_n  ))), nrow = folds_n, ncol = (steps_n  ))
  steploglikcv.p  = matrix( rep(0,folds_n*40), nrow=folds_n, ncol=40 ) ;        ## 40 steps for p from 0.01 to 0.40 
  stepagreecv.p   = matrix( rep(0,folds_n*40), nrow=folds_n, ncol=40 ) ; 
  step_dfcv.p     = matrix( rep(0,folds_n*40), nrow=folds_n, ncol=40 ) ; 

  ##===========================================================================================##
  
  ## i_cv = 1 
  for (i_cv in 1:folds_n) {
    if (track>=1) { cat(paste0(" ## Cross Validation fold  ", i_cv, "  of  " , folds_n , "  ##" , "\n")) }
    
    ##### set up train and test data sets in matrix form for glmnet & stepreg #####
    trainxs_cv = xs_cv[(folds_index_cv!=i_cv),]   
    testxs_cv  = xs_cv[(folds_index_cv==i_cv),]  
    dim(trainxs_cv) 
    dim(testxs_cv) 
    dim(trainxs_cv)[1] + dim(testxs_cv)[1] 
    
    if (is.null(start_cv)) {
      trainstart_cv = NULL
      teststart_cv  = NULL
    } else {
      trainstart_cv = start_cv[(folds_index_cv!=i_cv)]    
      teststart_cv  = start_cv[(folds_index_cv==i_cv)]   
    }
    
    trainy_cv = y_cv[(folds_index_cv!=i_cv)]    
    testy_cv  = y_cv[(folds_index_cv==i_cv)]   
  
    if (family == "cox") {
      trainevent_cv = event_cv[(folds_index_cv!=i_cv)] ;  
      testevent_cv  = event_cv[(folds_index_cv==i_cv)] ; 
    } else {
      trainevent_cv = NULL
      testevent_cv  = NULL 
    }
  
    ###### START STEPWISE fits #####################################################
##    print( "HERE cv 2" ) 
##  xs_st=trainxs_cv ; start_time_st=trainstart_cv ; y_st=trainy_cv ; event_st=trainevent_cv ; steps_n=steps_n ; method=method ; family=family ;
    stepreg.fit = stepreg(xs_st=trainxs_cv, start_time_st=trainstart_cv, y_st=trainy_cv, event_st=trainevent_cv, steps_n=steps_n, method=method, family=family, track=track)  
#    class( stepreg.fit )
    class( stepreg.fit ) = "data.frame" 
#    typeof(stepreg.fit)
    ## get statistics (loglik and concordance or r-square) from CV test data
    ## one step at a time 

    stepreg.fit.best = stepreg.fit[stepreg.fit$best==1,]                                        ## identify those models that were best at each stage of the stepwise procedure
    stepreg.fit.best [ , 1:7]  
    
    ##========================================================================##

    for (j_cv in 1:steps_n) { 
      if (family == "cox") {
        testbeta_cv = stepreg.fit.best [ j_cv , (xs_cv_ncol+8):(2*xs_cv_ncol+7) ]   ## get beta, the regression estimate
        testbeta_cv[ is.na( testbeta_cv ) ] = 0                                   ## set terms missing, because of overparamterization, to 0
        testscore_cv  = testxs_cv %*% t( testbeta_cv )                            ## get predicteds 
      } else { 
        testbeta_cv = stepreg.fit.best [ j_cv , (xs_cv_ncol+8):(2*xs_cv_ncol+8) ]   ## get beta, the regression estimate, including intercept 
        testbeta_cv[ is.na( testbeta_cv ) ] = 0                                   ## set terms missing, because of overparamterization, to 0
        testscore_cv  = cbind(1,testxs_cv) %*% t( testbeta_cv )                            ## get predicteds 
      }
      
      testscore_cv1 = as.numeric(testscore_cv)
#      testbeta_cv[1:10]

#      testscore_cv  = testxs_cv %*% t( stepreg.fit.best [ j_cv , (xs_cv_ncol+8):(2*xs_cv_ncol+7) ] ) 
      
      if (family == "cox") {
        if (is.null(start_cv)) { 
          testdata_cv   = as.data.frame( cbind(testy_cv, testevent_cv, testscore_cv) )
          colnames(testdata_cv)[3] = c("testscore_cv")                            ## why no name of testscore_cv by default ??
        } else {                                                                  ## why doesn;t the assign work ??
          testdata_cv   = as.data.frame( cbind(teststart_cv, testy_cv, testevent_cv, testscore_cv) )
#          print( colnames(testdata_cv) )
          colnames(testdata_cv)[4] = c("testscore_cv") 
        } 
      } else {
        testdata_cv   = as.data.frame( cbind(testy_cv, testscore_cv) )
        names(testdata_cv)[2] = "testscore_cv"
      }
#      testdata_cv[1:8,]     
      ## test data for SRTEPwise model fit 
      
#      colnames(testdata_cv) = c("testy_cv", "testevent_cv", "testscore_cv") 
      colnames(testdata_cv)
      nx = sum(is.na(testscore_cv)==1)
      
#      print( "HERE cv 3" )
      if (nx > 0) {  cat(paste0("\n", "n non missing", nx ))  }
#      print( testdata_cv[1:4,] )
      
      if (family =="cox") {
        if (is.null(start_cv)) {  fittest = coxph(Surv(testy_cv, testevent_cv) ~ testscore_cv, data=testdata_cv)  
        } else {    fittest = coxph(Surv(teststart_cv, testy_cv, testevent_cv) ~ testscore_cv, data=testdata_cv) }
      
        if (j_cv==1) {
          steploglikcv.df[i_cv, j_cv] = fittest$loglik[1] / fittest$n  
        }
        steploglikcv.df[i_cv, j_cv+1] = fittest$loglik[2] / fittest$n  
        stepagreecv.df [i_cv, j_cv  ] = fittest$concordance[6]   
        step_ncv      [i_cv , j_cv  ] = fittest$n
      }
      
      if ( family == "binomial" ) {
        fittest = glm(testy_cv ~ testscore_cv, data=testdata_cv, family="binomial")  
        loglik = c(loglik.null=-fittest$null.deviance/2, loglik=-fittest$deviance/2)
        clist  =  concordance(fittest)
        concord = c(concordance=as.numeric( clist[1] ), se=sqrt( as.numeric( clist[4] ) ) )  

        if (j_cv==1) {
          steploglikcv.df[i_cv, j_cv] = loglik[1] / dim(testdata_cv)[1]  
        }
        steploglikcv.df[i_cv, j_cv+1] = loglik[2] / dim(testdata_cv)[1]  
        stepagreecv.df [i_cv, j_cv  ] = concord[1]   
        step_ncv       [i_cv, j_cv  ] = dim(testdata_cv)[1] 
#        loglik 
#        concord 
      }  
      
      if ( family == "gaussian" ) {
        fittest = lm(testy_cv ~ testscore_cv)         ###### and in other calls like this ###### , weights=weights_val
#        fittest = glm(testy_cv ~ testscore_cv, data=testdata_cv, family="gaussian")          
#        summary(fittest) ; logLik(fittest)
        n = length(testy_cv) 
        kp1 = length(fittest$coefficients) 
        
        S2MLE  =  (t(fittest$residuals)%*%fittest$residuals)/n  ;              S2MLE 
        loglik = -(n/2)*log(2*pi) - n*log(sqrt(S2MLE)) - (1/(2*S2MLE))*n*S2MLE 
        S2 = var(y_cv)*(n-1)/n
        loglik.null = -(n/2)*log(2*pi) - n*log(sqrt(S2)) - (1/(2*S2))*n*S2 
        loglik = c(loglik.null=loglik.null, loglik=loglik)
#        loglik
#        logLik(fittest)
        rsquare = var(fittest$fitted.values)/var(y_cv) 
        rsquareadj = 1 - ((1-rsquare)*(n-1)/(n-kp1)) 
        rsquare = c(rsquare=rsquare, rsquareadj=rsquareadj)
        if (j_cv==1) {
          steploglikcv.df[i_cv, j_cv] = loglik[1] / dim(testdata_cv)[1]  
        }
        steploglikcv.df[i_cv, j_cv+1] = loglik[2] / dim(testdata_cv)[1]  
        stepagreecv.df [i_cv, j_cv  ] = rsquare[1]   
        step_ncv       [i_cv, j_cv  ] = dim(testdata_cv)[1] 
      }  
      
    }  ##### end j_cv loop #####
    
    ##========================================================================##
    
    ##===== GET FIT AS FUNCTION OF P CRITICAL VALUE FOR ENTRY  ===============##
    
    stepreg.fit.best.p = as.matrix( stepreg.fit.best )                          ## same as stepreg.fit.best but in matrix format 
    
#    stepreg.fit.best [ , 1:7]
    stepreg.fit.best.p [ , 1:7]   
    
    logliklast = c(stepreg.fit.best.p[1,4], (stepreg.fit.best.p[1:(dim(stepreg.fit.best.p)[1]-1),5]) )
    logliklast
    twologlikdif  = 2 * (stepreg.fit.best.p[,5] - logliklast)
    class( twologlikdif )
    twologlikdif 
    #twologlikdifx = as.matrix( twologlikdif )
    #twologlikdifx = as.numeric( twologlikdifx  )
    #class( twologlikdifx )
    #pvalue = 1 - pchisq(as.vector( twologlikdifx ), df=1 )
    pvalue = 1 - pchisq(as.vector( twologlikdif ), df=1 )
    pvalue
    class ( pvalue ) 
    # t(pvalue)
#    plot(1:length(pvalue), pvalue)
    
    psteps = 40 
    psteps
#    stepreg.fit.best.p.p = matrix(data=rep(1, psteps * (dim(stepreg.fit.best.p)[2]+2)), nrow=psteps, ncol = (dim(stepreg.fit.best.p)[2]+2) )
#    stepreg.fit.best.p.p = matrix(data=rep(1, psteps * 9), nrow=psteps, ncol = 9 )
#    dim(stepreg.fit.best.p.p)
#    stepreg.fit.best.p.p
    
    # steps_n =5 
    
    for (pstep_ in 1:psteps){ 
      pcrit = pstep_/100 
#      pgpcriti = ( pvalue > pcrit)
#      pgpcriti 
#      pgpcrit = sum(( pvalue > pcrit)*1)
      pgpcritsum = sum(( pvalue > pcrit)*1)
      pgpcritsum 
      if (pgpcritsum == 0) {  df = steps_n   
      } else {
        df = min( c(1:steps_n)[  ( pvalue > pcrit)  ] ) - 1
      }
      df 
      steploglikcv.p[i_cv, pstep_] = steploglikcv.df[i_cv, df+1 ]  
      if (df>0) { stepagreecv.p [i_cv, pstep_] = stepagreecv.df[i_cv, df  ]  
      } else if (family %in% c("cox","binomial")) { stepagreecv.p [i_cv, pstep_] = 0.5 }
      step_dfcv.p   [i_cv, pstep_] = df
#      cat(paste0( pstep_, "  pcrit=", pcrit, "  df=", df, "  C=", round(stepagreecv.p[i_cv, pstep_], digits=4), "\n"))
    } 
    
    steploglikcv.p
    stepagreecv.p
    step_dfcv.p
    
#    cat(paste0("\n", " loglik and concordance or r-square as function of steps in stepwise fit for fold  ", i_cv, "\n"))
#    print( round( 1000*steploglikcv.df[i_cv,])/1000 ) 
#    print( round( 1000*stepagreecv[i_cv,])/1000 ) 
#    cat(paste0("\n"))
  }  ##### END i_cv loop #####
  
  ##===========================================================================================##

  if (track>=1) { cat(paste0(" ## Stepwise fit on full data set" , "\n")) } 
  stepreg.fit.all = stepreg(xs_st=xs_cv, start_time_st=start_cv, y_st=y_cv, event_st=event_cv, steps_n=steps_n, family=family, track=track)  
  class(stepreg.fit.all) = "data.frame" 
  stepreg.fit.all.best = stepreg.fit.all[(stepreg.fit.all$best==1) , ] ; 
  ##### get model specification same as that of in stepreg() function ##### 
  #  stepreg.fit.all.best[,1:7]
  #  cat(paste0( " dim(stepreg.fit.all)     = " ,  dim(stepreg.fit.all    )[1] , "  ", dim(stepreg.fit.all    )[2] ,  "\n") )
  #  cat(paste0( " dim(stepreg.fit.all.best) = " ,  dim(stepreg.fit.all.best)[1] , "  ", dim(stepreg.fit.all.best)[2] , "\n") )
  ##  print( stepreg.fit.all[1:7]  )
  ##  print( stepreg.fit.all.best[best.df,] ) 
  #  cat(paste0( " stepreg.fit.all.best = " , "\n", stepreg.fit.all.best ) )

  ##### get model fit as function of number of terms in model #############################################
  #########################################################################################################  
  meanloglikcv.df = colMeans( steploglikcv.df )
  meanagreecv.df  = colMeans( stepagreecv.df ) 
  mean_ncv        = colMeans( step_ncv) 
  
  meanloglikcv.df
  meanagreecv.df
  
#  plot(1:length( meanloglikcv.df ), meanloglikcv.df )
#  plot(1:length( meanagreecv.df ), meanagreecv.df )

  ##### get number of model terms maximizing fit ##################
  if (toupper(method)=="LOGLIK") {
    best.df = which.max( meanloglikcv.df[2:steps_n] )  
  } else {
    best.df = which.max( meanagreecv.df ) 
  } 
#  cat(paste0("best.df = " , best.df, "   xs_cv_ncol = " , xs_cv_ncol, "\n")) 
#  cat(summary(stepreg.fit.all.best[best.df,]))
  
  ##### get model based upon "best" number of model terms ###########
#  cvfit = as.vector( stepreg.fit.all.best[best.df,] )
  cvfit.df =            stepreg.fit.all.best[best.df,] 
#  print( cvfit.df  )

  if (best.df == 1) {
    bestvars = colnames(xs_cv)[cvfit.df[1,8:(xs_cv_ncol+7)]==1]                  
  } else {
    bestvars = colnames(xs_cv)[cvfit.df[8:(xs_cv_ncol+7)]==1]
  }
  bestvars

  if (family == "cox") {
    if (is.null(start_cv)) { 
      dataf = as.data.frame( cbind( y_cv, event_cv, xs_cv  ) )
      form1 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2] , ") ~ ", paste(bestvars, collapse = " + " ) ) )  
    } else { 
      dataf = as.data.frame( cbind( start_cv, y_cv, event_cv, xs_cv  ) ) 
      form1 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2],  ", " , colnames(dataf)[3] , ") ~ ", paste(bestvars, collapse = " + " ) ) )  
    }
    func.fit.df = coxph(form1, dataf)
  } else {
    dataf = as.data.frame( cbind( y_cv, xs_cv  ) )
    form1 = formula( paste(colnames(dataf)[1] , " ~ ", paste(bestvars, collapse = " + " ) ) )  
    func.fit.df = glm(form1, data=dataf, family=family)
  }

  ##### get model fit as function of p cutoff value #######################################################
  #########################################################################################################  
  
  meanloglikcv.p = colMeans( steploglikcv.p ) 
  meanagreecv.p  = colMeans( stepagreecv.p ) 
#  meandfcv.p     = colMeans( step_dfcv.p ) 

  if (toupper(method)=="LOGLIK") {
    best.p    = which.max( meanloglikcv.p )
    maxloglik = max(meanloglikcv.p)
    for (k_ in c(1:psteps) ) {
      if ( meanloglikcv.p[k_] >= maxloglik) { best.p = k_}
    }
  } else { 
    best.p    = which.max( meanagreecv.p )
#    maxagree = max( meanagreecv.p )
#    for (k_ in c(1:psteps) ) {
#      if ( meanagreecv.p[k_] >= maxagree) { best.p = k_}
#    }
  }
    
#  cat(paste0("best p critical = " , best.p/100, "   xs_cv_ncol = " , xs_cv_ncol, "\n")) 
#  pcrit=(c(1:psteps)/100) ; AveLogLik = meanloglikcv.p ; plot(pcrit, AveLogLik) 
#  pcrit=(c(1:psteps)/100) ; AveConcordance = meanconcorcv.p ; plot(pcrit, AveConcordance) 

  logliklast = c(stepreg.fit.all.best[1,4], (stepreg.fit.all.best[1:(dim(stepreg.fit.all.best)[1]-1),5]) )
  logliklast
  twologlikdif  = 2 * (stepreg.fit.all.best[,5] - logliklast)
  class( twologlikdif )
  twologlikdif 
  pvalue = 1 - pchisq(as.vector( twologlikdif ), df=1 )
  pvalue
  class ( pvalue ) 
  
  cbind( stepreg.fit.all.best [ , 1:5], p=round(pvalue, digits=4), stepreg.fit.all.best [ , 6:7])

  psteps = 40 

  pstep_ = best.p              ## used in do loop above - resue code decimal p entry
    pcrit = pstep_/100 
    pgpcritsum = sum(( pvalue > pcrit)*1)
    pgpcritsum 
    if (pgpcritsum == 0) {  df.p = steps_n   
    } else {
      df.p = min( c(1:steps_n)[  ( pvalue > pcrit)  ] ) - 1
    }
    df.p     

  cvfit.p =  stepreg.fit.all.best[df.p,] 

  if (best.df == 1) {
    bestvars = colnames(xs_cv)[cvfit.p[1,8:(xs_cv_ncol+7)]==1]
  } else {
    bestvars = colnames(xs_cv)[cvfit.p[8:(xs_cv_ncol+7)]==1]
  }
#  print( bestvars ) 
  
  if (family == "cox") {
    if (is.null(start_cv)) { 
      dataf = as.data.frame( cbind( y_cv, event_cv, xs_cv  ) )
      form1 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2] , ") ~ ", paste(bestvars, collapse = " + " ) ) )  
    } else { 
      dataf = as.data.frame( cbind( start_cv, y_cv, event_cv, xs_cv  ) ) 
      form1 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2],  ", " , colnames(dataf)[3] , ") ~ ", paste(bestvars, collapse = " + " ) ) )  
    }
    func.fit.p = coxph(form1, dataf)
  } else {
    dataf = as.data.frame( cbind( y_cv, xs_cv  ) )
    form1 = formula( paste(colnames(dataf)[1] , " ~ ", paste(bestvars, collapse = " + " ) ) )  
    func.fit.p = glm(form1, data=dataf, family=family)
  }
  
#  summary( func.fit.p ) 
#  cat(paste0("best p critical = " , best.p/100, "  best.df = " , best.df, "  Number of candidate predictors = " , xs_cv_ncol, "\n")) 
  #########################################################################################################
  ## (xs_cv, start_cv=NULL, y_cv, event_cv, steps_n=0, folds_n=10, method="loglik", family="cox") 
  ## cvfit - modsum row for chosen row of stepreg.fit.all.best  
  
  Call <- match.call()
#  indx <- match(c("x_val","y_val","xs_val","ys_val","xs_non","ys_non","jksize","vmethod","familyr",
#                  "e_val","es_val","es_non","id_val","id_non"),
#                names(Call), nomatch=0)
#  print( Call )
  
  returnlist = list( Call = Call, 
                     tuning=c(steps_n=steps_n, folds_n=folds_n, method=method, family=family), 
                     best.df = best.df,  best.p=best.p/100, df.p=df.p, 
                     func.fit.df=func.fit.df, func.fit.p=func.fit.p, 
                     steploglikcv.df = steploglikcv.df, stepagreecv.df =stepagreecv.df, 
                     steploglikcv.p  = steploglikcv.p , stepagreecv.p  =stepagreecv.p , step_dfcv.p=step_dfcv.p, 
                     cvfit.df=cvfit.df, cvfit.p=cvfit.p,  
                     stepreg.fit.all.best=stepreg.fit.all.best , pvalue=pvalue,
                     seed=seed, foldid=folds_index_cv)            ## stepreg.fit.all=stepreg.fit.all,
  class(returnlist) <- c("cv.stepreg")
  return( returnlist ) 
}

################################################################################################################################################################
################################################################################################################################################################

#' Summarize results from a cv.stepreg() output object.  
#'
#' @param object A cv.stepreg() output object
#' @param ... Additional arguments passed to the summary function.   
#'
#' @return Summary of a stepreg() (stepwise regression) output object.   
#' 
#' @seealso 
#'    \code{\link{predict.cv.stepreg}}  , \code{\link{cv.stepreg}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#' 
#' @examples 
#' set.seed(955702213)
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=c(0,1,1))
#' # this gives a more interesting case but takes longer to run
#' xs=sim.data$xs           
#' # this will work numerically as an example 
#' xs=sim.data$xs[,c(2,3,50:55)] 
#' dim(xs)
#' y_=sim.data$yt 
#' event=sim.data$event
#' # for this example we use small numbers for steps_n and folds_n to shorten run time 
#' cv.stepreg.fit = cv.stepreg(xs, NULL, y_, event, steps_n=10, folds_n=3, track=0)
#' summary(cv.stepreg.fit)
#' 
#' 
summary.cv.stepreg = function(object, ...) {
  data1 = object$stepreg.fit.all.best ## [(object$best==1),]
  if ( ( dim(data1)[2] %% 2) ==0 ) { cox=0 } else { cox=1 }
  if (cox==1) { 
    k_ = dim(data1)[2]/2 - 3.5
    keepvars = c(1:k_)[(colSums(abs(data1[,(k_+8):(2*k_+7)]))!=0)]              ## + 7 + k_ 
  } else {
    k_ = dim(data1)[2]/2 - 4
    keepvars = c(1:(k_+1))[(colSums(abs(data1[,(k_+8):(2*k_+8)]))!=0)]          ## + 7 + k_ 
  }
  keepvars = keepvars + 7 + k_ 
  data2 = cbind( data1[,c(1:5)], pvalue=object$pvalue, data1[,c(6:7,keepvars)] )
  
  best.df = object$best.df  
  best.p  = object$best.p 
  df.p    = object$df.p 
  cat(paste0("\n", " CV best df = ", best.df, ", CV best p enter = ", best.p, " for ", df.p, " predictors \n     in the full data model, from ", k_ , " candidate predictors", "\n\n")) 
  betas=data2[c(best.df, df.p),-c(2,3)]
  colnames(betas)[1] = "df"
  rownames(betas) = NULL 
  print( betas[,(colSums(abs(betas))>0)] ) 
}

################################################################################################################################################################
################################################################################################################################################################

#' Beta's or predicteds based upon a cv.stepreg() output object. 
#' 
#' @description 
#' Give predicteds or Beta's based upon a cv.stepreg() output object. If an input data matrix is 
#' specified the X*Beta's are output.  If an input data matrix is not specified then 
#' the Beta's are output.  In the first column values are given based upon df as a tuning 
#' parameter and in the second column values based upon p as a tuning parameter.
#'
#' @param object cv.stepreg() output object 
#' @param xs dataset for predictions.  Must have the same columns as the input predictor matrix in the call to cv.stepreg(). 
#' @param ... pass through parameters 
#' 
#' @return a matrix of beta's or predicteds
#' 
#' @seealso 
#'    \code{\link{summary.cv.stepreg}}, \code{\link{cv.stepreg}} , \code{\link{nested.glmnetr}}   
#' 
#' @export
#' 
predict.cv.stepreg = function(object, xs=NULL, ...) {
  data1 = object$stepreg.fit.all.best ## [(object$best==1),]
  if ( ( dim(data1)[2] %% 2) ==0 ) { cox=0 } else { cox=1 }
  if (cox==1) { 
    k_ = dim(data1)[2]/2 - 3.5
    data2 = data1[ ,(k_+8):(2*k_+7)] 
  } else {
    k_ = dim(data1)[2]/2 - 4
    data2 = data1[ ,(k_+8):(2*k_+8)] 
  }

  best.df = object$best.df  
  best.p  = object$best.p 
  df.p    = object$df.p 
#  cat(paste0("\n", " CV best df = ", best.df, ", CV best p enter = ", best.p, " for ", df.p, " predictors \n     in the full data model, from ", k_ , " candidate predictors", "\n\n")) 
  betas = data2[c(best.df, df.p),]
  rownames(betas) = c("df","p")  
  betas = t(betas) 
  if (!is.null(xs)) { 
    if (cox!=1) {
      xs = cbind(rep(1,dim(xs)[1]),xs)
      colnames(xs)[1] = "Int" 
    }
    rturn = xs %*% betas 
  } else{ 
    rturn = betas 
  } 
  return(rturn)  
}

# object = cv.stepwise.fit ; xs=NULL ; 


