% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topicModel_utils.R
\name{topicNW1}
\alias{topicNW1}
\title{network showing top terms of topics}
\usage{
topicNW1(betaDF, topn)
}
\arguments{
\item{betaDF}{topic~term probability data.frame}

\item{topn}{number of top terms of each topic}
}
\description{
network showing top terms of topics
}
\examples{
\dontrun{
ldaOut <- SeuratObj@misc$ldaOut
betaDF <- tidytext::tidy(ldaOut, matrix = "beta")
topicNW1(betaDF, topn=10)
}

}
