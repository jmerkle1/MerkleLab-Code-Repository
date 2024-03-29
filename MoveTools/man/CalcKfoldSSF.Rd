% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalcKfoldSSF.R
\name{CalcKfoldSSF}
\alias{CalcKfoldSSF}
\title{Conduct a k-fold cross validation analysis for SSF}
\usage{
CalcKfoldSSF(
  data,
  formula,
  modType = "clogit",
  binVar,
  k = 5,
  resp = "case",
  strat = "strata",
  setSeed = TRUE,
  TMBtheta = log(1000),
  TMBmaparg = list(theta = factor(c(NA, 1:1))),
  update = TRUE,
  numb.top.ranks = 1
)
}
\arguments{
\item{data}{The dataset for which models should be run (should be a regular dataframe)}

\item{formula}{The formula (as a character string) of the model to be tested. Include response, explanatory variables, random effects, cluster term etc. Do NOT include a "data" argument,and do not include function name}

\item{modType}{The type of model you want to fit. One of: "clogit" (default), or "glmmTMB".}

\item{binVar}{The variable that you want to bin observations by. Typically animal ID. Must >= than K.}

\item{k}{Number of training and testing bins to use. Default = 5.}

\item{resp}{The name of the response variable in your data. Default = "case"}

\item{strat}{The name of the column representing strata in your data. Default = "strata"}

\item{setSeed}{Whether or not to set the randomization to the same each run (default = True)}

\item{TMBtheta}{The value to pass along to "TMBStruc$parameters$theta\link{1}" before fitting a glmmTMB model. This argument is used for the Muff et al (2020) method to fit conditional mixed effects models. Defaults to log(1e3).}

\item{TMBmaparg}{The value to pass along to "TMBStruc$mapArg" before fitting a glmmTMB model. This argument is used for the Muff et al (2020) method to fit conditional mixed effects models. Defaults to list(theta=factor(c(NA,1:1))), which is 1 random effect.}

\item{update}{Whether function should update on it's progress. Default is TRUE}

\item{numb.top.ranks}{How many top ranks to include in the calculate of proportion of used steps in top ranks Defaults to 1, which is the top bin.}
}
\value{
Returns a dataframe of each fold, cross-validation scores (i.e., spearman rank correlations) for each fold, and the proportion of used steps in the specified number of the top ranks.
}
\description{
Performs cross-validation and generates the cross-validation scores based on the number of used steps in each rank within each strata. Also calculates the proportion of used steps in a specified number of the top ranks. Based on the script from Roberts et al 2017 Ecography, updated by M. Laforge and J. Merkle to follow the methods of Fortin et al. 2009 (Ecology)
}
\examples{
#To come

}
