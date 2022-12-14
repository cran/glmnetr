% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.glmnetr_221210.R
\name{cv.glmnetr}
\alias{cv.glmnetr}
\title{Get a cross validated tuned relaxed lasso model fit.}
\usage{
cv.glmnetr(
  xs,
  start,
  y_,
  event,
  family = "cox",
  lambda = NULL,
  gamma = c(0, 0.25, 0.5, 0.75, 1),
  folds_n = 10,
  limit = 2,
  fine = 0,
  time = 0,
  seed = NULL,
  foldid = NULL
)
}
\arguments{
\item{xs}{Predictor matrix}

\item{start}{Vector of start times or the Cox model.  May be NULL.}

\item{y_}{outcome vector}

\item{event}{event vector in case of the Cox model.  May be NULL for other models.}

\item{family}{Model family, one of "cox", "gaussian" or "binomial".}

\item{lambda}{The lambda vector.  May be NULL.}

\item{gamma}{The gamma vector.  Default is c(0,0.25,0.50,0.75,1).}

\item{folds_n}{Number of folds for cross validation.  Default and recommended is 10.}

\item{limit}{limit the small values for lambda after the initial fit.
This will calcualtions that have minimal impact on the cross validation.  
Default is 2 for moderate limitation, 1 for less limitation, 0 for none.}

\item{fine}{Use a finer step in determining lambda.  Of little value unless one 
repeats the cross valiaiton many times to more finely tune the hyper parameters.  
See the _glmnet_ documentation.}

\item{time}{Indicate whether or not to update progress in the console.  Default of
0 suppresses these updates.  The option of 1 provides these updates.  In fitting 
clinical data with non full rank design matrix we have found some R-packages to
take a vary long time or seemingly be caught in infinite loops.  Therefore we allow
the user to track the package and judge whether things are moving forward or 
if the process should be stopped.}

\item{seed}{A seed for set.seed to assure one can get the same results twice.  If NULL 
the program will generate a random seed.  Whether specified or NULL, the seed is stored in the output
object for future reference.}

\item{foldid}{A vector of integers to associate each record to a fold.  Should be integers between 1 and folds_n.}
}
\value{
A cross validated relaxed lasso model fit.
}
\description{
This function derives a relaxed lasso model and derives hyperparmaters 
using cross validaiton.  It is analogous to the glmnet() function of the 
_glmnet_ package, but handles caes where glmnet() may run slowly when using the
relaxed-=TRUE option.
}
\details{
#' This is main program for model derivation.  As currently implemented the 
package requires the data to be input as #' vectors and matrices with no missing values 
(NA).  All data vectors and matrices must be numerical.  For categorical variables one
should first construct corresponding numerical variables to represent these 
categories.  To take advantage of the lasso model, one can use one hot coding
assigning an indicator for each level of each categorical varaible, or creating 
as well other contrasts variables suggested by the subject matter.
}
\examples{
# set seed for random numbers, optionally, to get reproducible results
set.seed(82545037)
sim.data=glmnetr.simdata(nrows=100, ncols=100, beta=NULL)
xs=sim.data$xs 
y_=sim.data$y_ 
event=sim.data$event
# for this example we use a small number for folds_n to shorten run time 
cv.glmnetr.fit = cv.glmnetr(xs, NULL, y_, NULL, family="gaussian", folds_n=3, limit=2) 
plot(cv.glmnetr.fit)
plot(cv.glmnetr.fit, coefs=1)
summary(cv.glmnetr.fit)

}
\seealso{
\code{\link{glmnetr}} , \code{\link{nested.glmnetr}}  , \code{\link{glmnetr.simdata}}
}
\author{
Walter Kremers (kremers.walter@mayo.edu)
}
