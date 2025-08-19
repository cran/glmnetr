#' @noRd
cv.glmnetr_perf = function(cv.glmnetr.list, testxs, testy__, trainxs, trainy__, family, ensemble, lasso_nms, lasso_xb_nms, tol=1e-14 ) {
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

  cv_elastic_g0 = pred_elastic_gam_fix(cv.glmnetr.list, 0, xs=testxs)  ## $xb   x=cv.glmnetr.list ; gam = 0 ; xs=testxs ; 
  cv_elastic_g1 = pred_elastic_gam_fix(cv.glmnetr.list, 1, xs=testxs)  ## $xb
  
  cbind( cv_elastic_g0$beta , cv_elastic_g1$beta ) 
  
  names(cv_glmnet_fit_)
  cv_elastic_g0_fit = cv_glmnet_fit_[[ cv_elastic_g0$alpha_i ]]
  cv_elastic_g1_fit = cv_glmnet_fit_[[ cv_elastic_g1$alpha_i ]]
  
  predminELg0  = cv_elastic_g0$xb
  predminELg1  = cv_elastic_g1$xb
#  predminELg1  = predict.cv.glmnetr(cv_elastic_g1_fit, testxs, lambda=cv_elastic_g1$lambda , gamma=cv_elastic_g0_fit$lambda, comment=0 )
  
  predmin.tr   = predict.cv.glmnetr(cv_lasso_fit  , trainxs, lambda=cv_lasso_fit$lambda.min, gamma=1, comment=0)  ##### minimizing lambda model predictions #########
  predminR.tr  = predict.cv.glmnetr(cv_lasso_fit  , trainxs, lambda="lambda.min" , gamma="gamma.min", comment = 0 )  ###### min RELAXED lasso ##########################
  predminR0.tr = predict.cv.glmnetr(cv_lasso_fit  , trainxs, lambda=cv_lasso_fit$relaxed$lambda.min.g0, gamma=0, comment=0)  ###### min gamma = 0 RELAXED lasso #####
  predminEL.tr = predict.cv.glmnetr(cv_elastic_fit, trainxs, lambda="lambda.min" , gamma="gamma.min", comment = 0 )
  predridge.tr = predict(cv_ridge_fit  , trainxs, s="lambda.min" ) 
  predminELg0.tr = pred_elastic_gam_fix(cv.glmnetr.list, 0, xs=trainxs)$xb
  predminELg1.tr = pred_elastic_gam_fix(cv.glmnetr.list, 1, xs=trainxs)$xb
  
  predmin.cal   = cal_train_xbhat( predmin  , trainy__, predmin.tr  , family ) 
  predminR.cal  = cal_train_xbhat( predminR , trainy__, predminR.tr , family ) 
  predminR0.cal = cal_train_xbhat( predminR0, trainy__, predminR0.tr, family ) 
  predminEL.cal = cal_train_xbhat( predminEL, trainy__, predminEL.tr, family ) 
  predridge.cal = cal_train_xbhat( predridge, trainy__, predridge.tr, family ) 
  predminELg0.cal = cal_train_xbhat( predminELg0, trainy__, predminELg0.tr, family ) 
  predminELg1.cal = cal_train_xbhat( predminELg1, trainy__, predminELg1.tr, family ) 
  
  predmin.cal.tr   = cal_train_xbhat( predmin.tr  , trainy__, predmin.tr  , family ) 
  predminR.cal.tr  = cal_train_xbhat( predminR.tr , trainy__, predminR.tr , family ) 
  predminR0.cal.tr = cal_train_xbhat( predminR0.tr, trainy__, predminR0.tr, family ) 
  predminEL.cal.tr = cal_train_xbhat( predminEL.tr, trainy__, predminEL.tr, family ) 
  predridge.cal.tr = cal_train_xbhat( predridge.tr, trainy__, predridge.tr, family ) 
  predminELg0.cal.tr = cal_train_xbhat( predminELg0.tr, trainy__, predminELg0.tr, family ) 
  predminELg1.cal.tr = cal_train_xbhat( predminELg1.tr, trainy__, predminELg1.tr, family ) 

#  print( cor(cbind(predmin, predminR, predminR0, predminEL, predridge)) ) 
  
  perfm1 = perf_gen( testy__ , predmin   , family )
  perfm2 = perf_gen( testy__ , predminR  , family )
  perfm3 = perf_gen( testy__ , predminR0 , family )
  perfm4 = perf_gen( testy__ , predminEL , family )
  perfm5 = perf_gen( testy__ , predridge , family )
  perfm6 = perf_gen( testy__ , predminELg0 , family )
  perfm7 = perf_gen( testy__ , predminELg1 , family )
  
  perfm1.cal = perf_gen( testy__ , predmin.cal   , family ) 
  perfm2.cal = perf_gen( testy__ , predminR.cal  , family ) 
  perfm3.cal = perf_gen( testy__ , predminR0.cal , family ) 
  perfm4.cal = perf_gen( testy__ , predminEL.cal , family ) 
  perfm5.cal = perf_gen( testy__ , predridge.cal , family ) 
  perfm6.cal = perf_gen( testy__ , predminELg0.cal , family ) 
  perfm7.cal = perf_gen( testy__ , predminELg1.cal , family ) 
  
  lasso.devian     = c( perfm1[1] , perfm2[1] , perfm3[1] , perfm4[1] , perfm5[1] , perfm6[1] , perfm7[1] )
  lasso.agree      = c( perfm1[3] , perfm2[3] , perfm3[3] , perfm4[3] , perfm5[3] , perfm6[3] , perfm7[3] )
  lasso.intcal     = c( perfm1[4] , perfm2[4] , perfm3[4] , perfm4[4] , perfm5[4] , perfm6[4] , perfm7[4] )
  lasso.lincal     = c( perfm1[5] , perfm2[5] , perfm3[5] , perfm4[5] , perfm5[5] , perfm6[5] , perfm7[5] )
  
# lasso.cal.devian = c( perfm1[2] , perfm2[2] , perfm3[2] , perfm4[2] , perfm5[2] )   
  lasso.cal.devian = c( perfm1.cal[1] , perfm2.cal[1] , perfm3.cal[1] , perfm4.cal[1] , perfm5.cal[1] , perfm6.cal[1] , perfm7.cal[1] ) 
  lasso.cal.agree  = c( perfm1.cal[3] , perfm2.cal[3] , perfm3.cal[3] , perfm4.cal[3] , perfm5.cal[3] , perfm6.cal[3] , perfm7.cal[3] ) 
  lasso.cal.intcal = c( perfm1.cal[4] , perfm2.cal[4] , perfm3.cal[4] , perfm4.cal[4] , perfm5.cal[4] , perfm6.cal[4] , perfm7.cal[4] ) 
  lasso.cal.lincal = c( perfm1.cal[5] , perfm2.cal[5] , perfm3.cal[5] , perfm4.cal[5] , perfm5.cal[5] , perfm6.cal[5] , perfm7.cal[5] ) 
  
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
  lasso.nzero[6] = sum( (abs(cv_elastic_g0$beta) >= tol) )
  lasso.nzero[7] = sum( (abs(cv_elastic_g1$beta) >= tol) )
  names(lasso.nzero) = lasso_nms 
  
  ## calibrate lasso ##
  if (sum(ensemble[c(2:4,6:8)]) >= 1) {
    ofst = lasso.intcal[4] + lasso.lincal[4] * predminR 
    if (family=="cox") { ofst = ofst - mean(ofst) }
  } else { ofst = NULL }
  
  xbetas.lasso = cbind( predmin    , predminR    , predminR0    , predminEL    , predridge    , predminELg0    , predminELg1    ,
                        predmin.cal, predminR.cal, predminR0.cal, predminEL.cal, predridge.cal, predminELg0.cal, predminELg1.cal )
  
  xbetas.lasso.tr = cbind( predmin.tr    , predminR.tr    , predminR0.tr    , predminEL.tr    , predridge.tr    , predminELg0.tr    , predminELg1.tr    , 
                           predmin.cal.tr, predminR.cal.tr, predminR0.cal.tr, predminEL.cal.tr, predridge.cal.tr, predminELg0.cal.tr, predminELg1.cal.tr )
  
  colnames( xbetas.lasso  ) = lasso_xb_nms 
  
  returnlist = list( lasso.devian=lasso.devian, lasso.agree=lasso.agree, lasso.intcal=lasso.intcal, lasso.lincal=lasso.lincal, 
                     lasso.cal.devian=lasso.cal.devian, lasso.cal.agree=lasso.cal.agree, lasso.cal.intcal=lasso.cal.intcal, lasso.cal.lincal=lasso.cal.lincal, 
                     lasso.nzero=lasso.nzero, xbetas.lasso=xbetas.lasso, xbetas.lasso.tr=xbetas.lasso.tr, ofst=ofst ) 
  
  return( returnlist )
  
#  return( list( lasso.devian=lasso.devian, lasso.agree=lasso.agree, lasso.intcal=lasso.intcal, lasso.lincal=lasso.lincal, 
#                lasso.cal.devian=lasso.cal.devian, lasso.cal.agree=lasso.cal.agree, lasso.cal.intcal=lasso.cal.intcal, lasso.cal.lincal=lasso.cal.lincal, 
#                lasso.nzero=lasso.nzero, xbetas.lasso=xbetas.lasso, xbetas.lasso.tr=xbetas.lasso.tr, ofst=ofst ) ) 
}