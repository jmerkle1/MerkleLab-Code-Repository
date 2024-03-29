% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalcKfoldCV.R
\name{CalcKfoldCV}
\alias{CalcKfoldCV}
\title{Conduct a k-fold cross validation analysis}
\usage{
CalcKfoldCV(
  data,
  formula,
  modType = glmer,
  binVar,
  k = 5,
  resp = "Pres",
  modFam = "binomial",
  setSeed = T,
  TMBtheta = NULL,
  TMBmaparg = NULL,
  update = T
)
}
\arguments{
\item{data}{The dataset for which models should be run}

\item{formula}{The formula of the model to be tested. Include response and explanatory variables, do NOT include a "data" argument,and do not include function name (e.g.: "Pres~x1+x2")}

\item{modType}{The type of model you want to fit. One of: "glmer" (default), "glm","lmer", "mclogit", and "glmmTMB".}

\item{binVar}{The variable that you want to bin observations by. Typically animal ID. Must >= than K.}

\item{k}{Number of training and testing bins to use. Default = 5.}

\item{resp}{The name of the response variable in your data. Default = Pres}

\item{modFam}{The family of the model (default = binomial)}

\item{setSeed}{The value to pass along to "TMBStruc$parameters$theta\link{1}" before running the model. This argument is used for the Muff et al (2020) method to fit conditional mixed effects models.}

\item{TMBtheta}{Number of training and testing bins to use. Default = 5.}

\item{TMBmaparg}{The value to pass along to "TMBStruc$mapArg" both when setting of the funtion and before fitting. This argument is used for the Muff et al (2020) method to fit conditional mixed effects models. Seems to be necessary to get "predict" to work when setting the variance structure manually for this method.}

\item{update}{Whether function should update on it's progress. Default is T}
}
\value{
Returns a vector of cross-validation scores (i.e., spearman rank correlations) for each fold. Uses the area-adjusted frequency of categories (bins). See Boyce et al. 2002 for details.
}
\description{
Performs cross-validation and generates the cross-validation score. Based on the script from Roberts et al 2017 Ecography, updated by M. Laforge and J. Merkle
}
\examples{
#To come

}
