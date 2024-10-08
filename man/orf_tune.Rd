% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ObliqueRandomForest_240822.R
\name{orf_tune}
\alias{orf_tune}
\title{Fit a Random Forest model on data provided in matrix and vector formats.}
\usage{
orf_tune(
  xs,
  start = NULL,
  y_,
  event = NULL,
  family = NULL,
  mtryc = NULL,
  ntreec = NULL,
  nsplitc = 8,
  seed = NULL,
  tol = 1e-05,
  track = 0
)
}
\arguments{
\item{xs}{predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
and not a data frame.}

\item{start}{an optional vector of start times in case of a Cox model. Class numeric of length same as number of patients (n)}

\item{y_}{dependent variable as a vector: time, or stop time for Cox model, Y_ 0 or 1 for binomial (logistic), numeric for gaussian. 
Must be a vector of length same as number of sample size.}

\item{event}{event indicator, 1 for event, 0 for census, Cox model only.
Must be a numeric vector of length same as sample size.}

\item{family}{model family, "cox", "binomial" or "gaussian" (default)}

\item{mtryc}{a vector (numeric) of values to search over for optimization of the 
Random Forest fit.  This if for the mtry input variable of the orsf() program 
specifying the number of terms to consider in each step of teh Random Forest fit.}

\item{ntreec}{a vector (numeric) of 2 values, the first for the number of forests
(ntree from orsf()) to use when searhcing for a better bit and the second to use
when fitting the final model.  More trees should give a better fit but 
require more computations and storage for the final. 
model.}

\item{nsplitc}{This nsplit of orsf(), a non-negative integer for the  number of 
random splits for a predictor.}

\item{seed}{a seed for set.seed() so one can reproduce the model fit.  If 
NULL the program will generate a random seed.  Whether specified or NULL, the 
seed is stored in the output object for future reference.  Note,
for the default this randomly generated seed depends on the seed in memory at that 
time so will depend on any calls of set.seed prior to the call of this function.}

\item{tol}{a small number, a lower bound to avoid division by 0}

\item{track}{1 to output a brief summary of the final selected model, 2 to 
output a brief summary on each model fit in search of a better model or 0 
(default) to not output this information.}
}
\value{
a Random Forest model fit
}
\description{
Fit an Random Forest model using the orsf() function of the 
randomForestSRC package.
}
\seealso{
\code{\link{summary.orf_tune}} , \code{\link{rederive_orf}} , \code{\link{nested.glmnetr}}
}
\author{
Walter Kremers (kremers.walter@mayo.edu)
}
