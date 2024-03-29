% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.nested.glmnetr_240322.R
\name{glmnetr.compcv0}
\alias{glmnetr.compcv0}
\title{Calculate agreement differences with CI and p}
\usage{
glmnetr.compcv0(a, b, digits = 4, txt = 0, pow = 1)
}
\arguments{
\item{a}{One term}

\item{b}{A second term}

\item{digits}{digits for printing of z-scores, p-values, etc. with default of 4}

\item{txt}{1 (default) to include inline text for estimated, 95 percent CI and p}

\item{pow}{Power to which the average of correlations is to be raised.  Only 
applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
show correlations.  Pow is ignored for the family of "cox" and "binomial".}
}
\value{
An estimate, 95% CI and p for agreement comparison
}
\description{
Perform a paired t-test as called from glmnetr.compcv().
}
