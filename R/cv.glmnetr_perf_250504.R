#' @noRd
cv.glmnetr_perf = function(cv.glmnetr.list, testxs, testy__, trainxs, trainy__, family, ensemble, lasso_nms, lasso_xb_nms ) {
  cv_glmnet_fit_ = cv.glmnetr.list$cv_glmnet_fit_
  cv_lasso_fit   = cv.glmnetr.list$cv_lasso_fit 
  cv_ridge_fit   = cv.glmnetr.list$cv_ridge_fit   
  cv_elastic_fit = cv.glmnetr.list$cv_elastic_fit
  lambda_null    = cv.glmnetr.list$lambda_null
  
  predmin    = predict.cv.glmnetr(cv_lasso_fit  , testxs, lambda=cv_lasso_fit$lambda.min, gamma=1, comment = 0) 
  predminR   = predict.cv.glmnetr(cv_lasso_fit  , testxs, lambda=cv_lasso_fit$relaxed$lambda.min , gamma=cv_lasso_fit$relaxed$gamma.min , comment = 0 ) 
  predminR0  = predict.cv.glmnetr(cv_lasso_fit  , testxs, lambda=cv_lasso_fit$relaxed$lambda.min.g0, gamma=0, comment=0)  
  predminEL  = predict.cv.glmnetr(cv_elastic_fit, testxs, lambda="lambda.min" , gamma="gamma.min", comment=0 )
  predridge  = predict(cv_ridge_fit  , testxs, s="lambda.min")                  ## default is x="lambda.1se" 
#  predridge2  = predict.cv.glmnetr(cv_ridge_fit  , testxs)
#  print(cbind(predridge, predridge2)[1:5,])
  
  predmin.tr   = predict.cv.glmnetr(cv_lasso_fit  , trainxs, lambda=cv_lasso_fit$lambda.min, gamma=1, comment=0)  ##### minimizing lambda model predictions #########
  predminR.tr  = predict.cv.glmnetr(cv_lasso_fit  , trainxs, lambda="lambda.min" , gamma="gamma.min", comment = 0 )  ###### min RELAXED lasso ##########################
  predminR0.tr = predict.cv.glmnetr(cv_lasso_fit  , trainxs, lambda=cv_lasso_fit$relaxed$lambda.min.g0, gamma=0, comment=0)  ###### min gamma = 0 RELAXED lasso #####
  predminEL.tr = predict.cv.glmnetr(cv_elastic_fit, trainxs, lambda="lambda.min" , gamma="gamma.min", comment = 0 )
  predridge.tr = predict(cv_ridge_fit  , trainxs, s="lambda.min" ) 
  
  predmin.cal   = cal_train_xbhat( predmin  , trainy__, predmin.tr  , family ) 
  predminR.cal  = cal_train_xbhat( predminR , trainy__, predminR.tr , family ) 
  predminR0.cal = cal_train_xbhat( predminR0, trainy__, predminR0.tr, family ) 
  predminEL.cal = cal_train_xbhat( predminEL, trainy__, predminEL.tr, family ) 
  predridge.cal = cal_train_xbhat( predridge, trainy__, predridge.tr, family ) 
  
  predmin.cal.tr   = cal_train_xbhat( predmin.tr  , trainy__, predmin.tr  , family ) 
  predminR.cal.tr  = cal_train_xbhat( predminR.tr , trainy__, predminR.tr , family ) 
  predminR0.cal.tr = cal_train_xbhat( predminR0.tr, trainy__, predminR0.tr, family ) 
  predminEL.cal.tr = cal_train_xbhat( predminEL.tr, trainy__, predminEL.tr, family ) 
  predridge.cal.tr = cal_train_xbhat( predridge.tr, trainy__, predridge.tr, family ) 

#  print( cor(cbind(predmin, predminR, predminR0, predminEL, predridge)) ) 
  
  perfm1 = perf_gen( testy__ , predmin   , family )
  perfm2 = perf_gen( testy__ , predminR  , family )
  perfm3 = perf_gen( testy__ , predminR0 , family )
  perfm4 = perf_gen( testy__ , predminEL , family )
  perfm5 = perf_gen( testy__ , predridge , family )
  
  perfm1.cal = perf_gen( testy__ , predmin.cal   , family ) 
  perfm2.cal = perf_gen( testy__ , predminR.cal  , family ) 
  perfm3.cal = perf_gen( testy__ , predminR0.cal , family ) 
  perfm4.cal = perf_gen( testy__ , predminEL.cal , family ) 
  perfm5.cal = perf_gen( testy__ , predridge.cal , family ) 
  
  lasso.devian     = c( perfm1[1] , perfm2[1] , perfm3[1] , perfm4[1] , perfm5[1] )
  lasso.agree      = c( perfm1[3] , perfm2[3] , perfm3[3] , perfm4[3] , perfm5[3] )
  lasso.intcal     = c( perfm1[4] , perfm2[4] , perfm3[4] , perfm4[4] , perfm5[4] )
  lasso.lincal     = c( perfm1[5] , perfm2[5] , perfm3[5] , perfm4[5] , perfm5[5] )
  
# lasso.cal.devian = c( perfm1[2] , perfm2[2] , perfm3[2] , perfm4[2] , perfm5[2] )   
  lasso.cal.devian = c( perfm1.cal[1] , perfm2.cal[1] , perfm3.cal[1] , perfm4.cal[1] , perfm5.cal[1] ) 
  lasso.cal.agree  = c( perfm1.cal[3] , perfm2.cal[3] , perfm3.cal[3] , perfm4.cal[3] , perfm5.cal[3] ) 
  lasso.cal.intcal = c( perfm1.cal[4] , perfm2.cal[4] , perfm3.cal[4] , perfm4.cal[4] , perfm5.cal[4] ) 
  lasso.cal.lincal = c( perfm1.cal[5] , perfm2.cal[5] , perfm3.cal[5] , perfm4.cal[5] , perfm5.cal[5] ) 
  
  names(lasso.devian) = lasso_nms 
  names(lasso.agree)  = lasso_nms 
  names(lasso.intcal) = lasso_nms 
  names(lasso.lincal) = lasso_nms 
  
  names(lasso.cal.devian) = lasso_nms 
  names(lasso.cal.agree)  = lasso_nms 
  names(lasso.cal.intcal) = lasso_nms 
  names(lasso.cal.lincal) = lasso_nms 

  lasso.nzero = rep(0,5) 
  lasso.nzero[1] = cv_lasso_fit$nzero [ cv_lasso_fit$index[1] ]     
  lasso.nzero[2] = cv_lasso_fit$relaxed$nzero.min 
  lasso.nzero[3] = cv_lasso_fit$relaxed$nzero.min.g0
  lasso.nzero[4] = cv_elastic_fit$relaxed$nzero.min
  lasso.nzero[5] = cv_ridge_fit$nzero[ cv_ridge_fit$index[1] ]
  names(lasso.nzero) = lasso_nms 
  
  ## calibrate lasso ##
  if (sum(ensemble[c(2:4,6:8)]) >= 1) {
    ofst = lasso.intcal[4] + lasso.lincal[4] * predminR 
    if (family=="cox") { ofst = ofst - mean(ofst) }
  } else { ofst = NULL }
  
  xbetas.lasso = cbind( predmin    , predminR    , predminR0    , predminEL    , predridge ,
                        predmin.cal, predminR.cal, predminR0.cal, predminEL.cal, predridge.cal )
  
  xbetas.lasso.tr = cbind( predmin.tr    , predminR.tr    , predminR0.tr    , predminEL.tr    , predridge.tr , 
                           predmin.cal.tr, predminR.cal.tr, predminR0.cal.tr, predminEL.cal.tr, predridge.cal.tr )
  
  colnames( xbetas.lasso  ) = lasso_xb_nms 
  
  returnlist = list( lasso.devian=lasso.devian, lasso.agree=lasso.agree, lasso.intcal=lasso.intcal, lasso.lincal=lasso.lincal, 
                     lasso.cal.devian=lasso.cal.devian, lasso.cal.agree=lasso.cal.agree, lasso.cal.intcal=lasso.cal.intcal, lasso.cal.lincal=lasso.cal.lincal, 
                     lasso.nzero=lasso.nzero, xbetas.lasso=xbetas.lasso, xbetas.lasso.tr=xbetas.lasso.tr, ofst=ofst ) 
  
  return( list( lasso.devian=lasso.devian, lasso.agree=lasso.agree, lasso.intcal=lasso.intcal, lasso.lincal=lasso.lincal, 
                lasso.cal.devian=lasso.cal.devian, lasso.cal.agree=lasso.cal.agree, lasso.cal.intcal=lasso.cal.intcal, lasso.cal.lincal=lasso.cal.lincal, 
                lasso.nzero=lasso.nzero, xbetas.lasso=xbetas.lasso, xbetas.lasso.tr=xbetas.lasso.tr, ofst=ofst ) ) 
}