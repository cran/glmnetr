################################################################################
##### predict_ann_tab_yymmdd.R #################################################
################################################################################
#' Get predicteds for an Artificial Neural Network model fit in nested.glmnetr()
#' 
#' All but one of the Artificial Neural Network (ANNs) fit by nested.glmnetr() are based 
#' upon a neural network model and input from a lasso model.  Thus a simple model(xs)
#' statement will not give the proper predicted values.  This function process information 
#' form the lasso and ANN model fits to give the correct predicteds.  Whereas the ann_tab_cv()
#' function ca be used to fit a model based upon an input data set it does not fit a 
#' lasso model to allow an informed starting point for the ANN fit.  The pieces fo this
#' are in nested.glmnetr().  To fit a cross validation (CV) informed ANN model fit
#' one can run nested.glmnetr() with folds_n = 0 to derive the full data models without 
#' doing a cross validation.   
#' 
#' @param object a output object from the nested.glmnetr() function  
#' @param xs new data of the same form used as input to nested.glmnetr() 
#' @param modl ANN model entry an integer from 1 to 5 indicating which "lasso informed" 
#' ANN is to be used for calculations.  The number corresponds to the position of the 
#' ensemble input from the nested.glmnetr() call.  The model must already be fit 
#' to calculate predicteds:   
#' 1 for ensemble[1] = 1, for model based upon raw data ; 
#' 2 for ensemble[2] = 1, raw data plus lasso predicteds as a predictor variable (features) ;
#' 4 for ensemble[3] = 1, raw data plus lasso predicteds and initial weights corresponding to offset and allowed to update ; 
#' 5 for ensemble[4] = 1, raw data plus lasso predicteds and initial weights corresponding to offset and not allowed to updated  ; 
#' 6 for ensemble[5] = 1, nonzero relaxed lasso terms ; 
#' 7 for ensemble[6] = 1, nonzero relaxed lasso terms plus lasso predicteds as a predictor variable (features) ;
#' 8 for ensemble[7] = 1, nonzero relaxed lasso terms plus lasso predicteds with initial weights corresponding to offset and allowed to update ; 
#' 9 for ensemble[8] = 1, nonzero relaxed lasso terms plus lasso predicteds with initial weights corresponding to offset and not allowed to update. 
#' 
#' @return a vector of predicteds
#' 
#' @seealso
#     \code{\link{predict_nested_xgb}} , \code{\link{predict_nested_rf}} , 
#'     \code{\link{ann_tab_cv}} , \code{\link{nested.glmnetr}}  
#' 
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @export
#'
predict_ann_tab = function(object, xs, modl=NULL) {
  
  ensemble = object$ensemble[1:8]
  
#  if (length(ensemble) < 8) { ensemble = c(ensemble, rep(0,8-length(ensemble))) }
  
  if (is.null(modl)) { 
    modl = max(ensemble*c(1:8))
    cat(paste("\n modl not spedified, modl = ", modl, " used \n"))
  }
  
#  if (modl %in% c(1,2,3,4,5,6,7,8)) { ensemble = object$ensemble } 
  
  if (ensemble[modl] == 0) {
      cat(paste(" Ensemble model", modl, " not run in nested CV.  Check ensemble specification. \n   "))
      cat(paste(" Ensemble: "))
      cat(paste(ensemble)) 
  } else { 
  
  family = object$sample[1]
  calbeta = object$ann_zb$calbeta 
  if ( modl > 1 ) {
    lassopred0 = predict(object$cv_lasso_fit, xs, lam="lambda.min" , gam="gamma.min", comment=0) ## V 0-6.1 cv_glmnet_fit to cv_lasso_fit
    lassobeta  = predict(object$cv_lasso_fit,     lam="lambda.min" , gam="gamma.min", comment=0) ## V 0-6.1 cv_glmnet_fit to cv_lasso_fit
    if (family == "cox") {
      lassopred = calbeta[1] + calbeta[2] * lassopred0                          ## this transforms to coxph() centered predicteds 
    } else if (family %in% c("binomial", "gaussian")) {
      lassopred = calbeta[1] + calbeta[2] * lassopred0 
    }
  }
  
  xs_means = object$ann_zb$xs_means 
  xs_sds_ = object$ann_zb$xs_sds_ 
  xs_c     = sweep(xs, 2, xs_means, FUN = "-") 
  xs_z     = sweep(xs_c, 2, xs_sds_, FUN = "/") 
  xs_z0    = xs_z[,(xs_sds_ > 1e-8)]
  # table(diag(cor(xs_z0))) 
  
  if (modl > 1) {
    if (family %in% c("binomial", "gaussian")) { 
      xs_z1 = xs_z[,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])] ]   
    } else if (family == "cox") {
      xs_z1 = xs_z[,(lassobeta[[1]] != 0)] 
    }
    xs_z0L = cbind(lasso=lassopred,xs_z0)    
    xs_z1L = cbind(lasso=lassopred,xs_z1)    
  }
  
  if        (modl == 1) {  
    xs_t = torch_tensor(xs_z0, dtype=torch_float()) ; 
    ann_fit = object$ann_fit_1
  } else if (modl == 2) {  
    xs_t = torch_tensor(xs_z0L, dtype=torch_float()) ; 
    ann_fit = object$ann_fit_2
  } else if (modl == 3) {  
    xs_t = torch_tensor(xs_z0L, dtype=torch_float()) ; 
    ann_fit = object$ann_fit_3
  } else if (modl == 4) {  
    xs_t = torch_tensor(xs_z0L, dtype=torch_float()) ; 
    ann_fit = object$ann_fit_4
  } else if (modl == 5) {  
    xs_t = torch_tensor(xs_z1 , dtype=torch_float()) ; 
    ann_fit = object$ann_fit_5
  } else if (modl == 6) {  
    xs_t = torch_tensor(xs_z1L, dtype=torch_float()) ; 
    ann_fit = object$ann_fit_6
  } else if (modl == 7) {  
    xs_t = torch_tensor(xs_z1L, dtype=torch_float()) ; 
    ann_fit = object$ann_fit_7
  } else if (modl == 8) {  
    xs_t = torch_tensor(xs_z1L, dtype=torch_float()) ; 
    ann_fit = object$ann_fit_8
  } else if (modl == 0) { pred = lassopred }  
  
  family    = ann_fit$modelsumc[1]
  actv      = ann_fit$modelsum[5]
  drpot     = ann_fit$modelsum[6] 
  parmfinal = ann_fit$parmfinal 
  
  if        (actv == 1) { act_fn = nn_relu()      ; actvc = "relu" 
  } else if (actv == 2) { act_fn = nn_gelu()      ; actvc = "gelu" 
  } else if (actv == 3) { act_fn = nn_sigmoid()   ; actvc = "sigmoid" 
  } else if (actv == 4) { act_fn = nn_sigmoid()$mul(2)$add(-1) ; actvc = "arctan" 
  } else if (actv ==-1) { act_fn = nn_identity()  ; actvc = "identity" 
  } else                { act_fn = nn_relu()      ; actvc = "relu" 
  }
  
  if        (family == "binomial"  ) { lastact_fn = nn_sigmoid() 
  } else if (family == "multiclass") { lastact_fn = nn_softmax()
  } else if (family == "gaussian"  ) { lastact_fn = nn_identity()
  } else if (family == "cox"       ) { lastact_fn = nn_identity() 
  } else { lastact_fn = nn_identity() } 
  
  weight0 = torch_tensor(parmfinal$weight0_r, dtype=torch_float(), requires_grad=TRUE)
  weight3 = torch_tensor(parmfinal$weight3_r, dtype=torch_float(), requires_grad=TRUE)
  weight6 = torch_tensor(parmfinal$weight6_r, dtype=torch_float(), requires_grad=TRUE)
  bias0   = torch_tensor(parmfinal$bias0_r  , dtype=torch_float(), requires_grad=TRUE)
  bias3   = torch_tensor(parmfinal$bias3_r  , dtype=torch_float(), requires_grad=TRUE)
  bias6   = torch_tensor(parmfinal$bias6_r  , dtype=torch_float(), requires_grad=TRUE) 
  
  my_nn_linear0 <- nn_module(
#    classname = "my_nn_linear0",
    initialize = function(weight, bias) {
      self$weight <- nn_parameter(weight)
      self$bias <- nn_parameter(bias)
    },
    forward = function(x) {
      nnf_linear(x, self$weight, self$bias)
    }
  )
  
  model_ = nn_sequential(
    my_nn_linear0(weight0, bias0) ,
    act_fn ,  
    nn_dropout(p = drpot) ,
    my_nn_linear0(weight3, bias3) ,
    act_fn ,  
    nn_dropout(p = drpot) ,
    my_nn_linear0(weight6, bias6) ,
    lastact_fn 
  )    
  
  pred = as.numeric( model_(xs_t) ) 

  return(pred)
  
  } 
}

################################################################################
