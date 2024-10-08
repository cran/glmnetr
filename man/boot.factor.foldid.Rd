% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools_240822.R
\name{boot.factor.foldid}
\alias{boot.factor.foldid}
\title{Generate foldid's by 0/1 factor for bootstrap like samples where unique option between 0 and 1}
\usage{
boot.factor.foldid(event, fraction)
}
\arguments{
\item{event}{the outcome variable in a vector identifying the different potential 
levels of the outcome}

\item{fraction}{the fraction of the whole sample included in the bootstratp sample}
}
\value{
foldid's in a vector the same length as event
}
\description{
Generate foldid's by 0/1 factor for bootstrap like samples where unique option between 0 and 1
}
\seealso{
\code{\link{get.foldid}} , \code{\link{nested.glmnetr}}
}
