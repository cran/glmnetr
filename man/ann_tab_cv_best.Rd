% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ann_tab_cv_240422.R
\name{ann_tab_cv_best}
\alias{ann_tab_cv_best}
\title{Fit multiple Artificial Neural Network models on "tabular" provided as a matrix, and 
keep the best one.}
\usage{
ann_tab_cv_best(
  myxs,
  mystart = NULL,
  myy,
  myevent = NULL,
  myoffset = NULL,
  family = "binomial",
  fold_n = 5,
  epochs = 200,
  eppr = 40,
  lenz1 = 32,
  lenz2 = 8,
  actv = 1,
  drpot = 0,
  mylr = 0.005,
  wd = 0,
  l1 = 0,
  lasso = 0,
  lscale = 5,
  scale = 1,
  resetlw = 1,
  minloss = 1,
  gotoend = 0,
  bestof = 10,
  seed = NULL,
  foldid = NULL
)
}
\arguments{
\item{myxs}{predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
and not a data frame.}

\item{mystart}{an optional vector of start times in case of a Cox model. Class numeric of length same as number of patients (n)}

\item{myy}{dependent variable as a vector: time, or stop time for Cox model, Y_ 0 or 1 for binomial (logistic), numeric for gaussian. 
Must be a vector of length same as number of sample size.}

\item{myevent}{event indicator, 1 for event, 0 for census, Cox model only.
Must be a numeric vector of length same as sample size.}

\item{myoffset}{an offset term to be ues when fitting the ANN.  Not yet implemented.}

\item{family}{model family, "cox", "binomial" or "gaussian" (default)}

\item{fold_n}{number of folds for each level of cross validation}

\item{epochs}{number of epochs to run when tuning on number of epochs for fitting 
final model number of epochs informed by cross validation}

\item{eppr}{for EPoch PRint.  print summry info every eppr epochs. 0 will 
print first and last epochs, -1 nothing.}

\item{lenz1}{length of the first hidden layer in the neural network, default 16}

\item{lenz2}{length of the second hidden layer in the neural network, default 16}

\item{actv}{for ACTiVation function.  Activation function between layers, 
1 for relu, 2 for gelu, 3 for sigmoid.}

\item{drpot}{fraction of weights to randomly zero out.  NOT YET implemented.}

\item{mylr}{learning rate for the optimization step in teh neural network model fit}

\item{wd}{weight decay for the model fit.}

\item{l1}{a possible L1 penalty weight for the model fit, default 0 for not considered}

\item{lasso}{1 to indicate the first column of the input matrix is an offset 
term, often derived from a lasso model}

\item{lscale}{Scale used to allow ReLU to extend +/- lscale before capping the 
inputted linear estimated}

\item{scale}{Scale used to transform the initial random parameter assingments by 
dividing by scale}

\item{resetlw}{1 as default to re-adjust weights to account for the offset every 
epoch.  This is only used in case lasso is set to 1}

\item{minloss}{default of 1 for minimizing loss, else maximizing agreement (concordance 
for Cox and Binomial, R-square for Gaussian), as function of epochs by cross validation}

\item{gotoend}{fit to the end of epochs.  Good for plotting and exploration}

\item{bestof}{how many models to run, from which the best fitting model will be selected.}

\item{seed}{an optional a numerical/integer vector of length 2, for R and torch 
random generators, default NULL to generate these.  Integers should be positive 
and not more than 2147483647.}

\item{foldid}{a vector of integers to associate each record to a fold.  Should 
be integers from 1 and fold_n.}
}
\value{
an artificial neural network model fit
}
\description{
Fit an multiple Artificial Neural Network models for analysis of "tabular" 
data using ann_tab_cv() and select the best fitting model according to cross
validaiton.
}
\seealso{
\code{\link{ann_tab_cv}} , \code{\link{predict_ann_tab}}, \code{\link{nested.glmnetr}}
}
\author{
Walter Kremers (kremers.walter@mayo.edu)
}
