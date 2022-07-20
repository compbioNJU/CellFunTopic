% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QualityControl.R
\name{QCfun}
\alias{QCfun}
\title{quality control according to 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'CellsPerGene'}
\usage{
QCfun(
  SeuratObj,
  minCountsPerCell = NULL,
  maxCountsPerCell = NULL,
  minFeaturesPerCell = NULL,
  maxFeaturesPerCell = NULL,
  maxPercent.mt = NULL,
  plot = TRUE
)
}
\arguments{
\item{SeuratObj}{Seurat object}

\item{minCountsPerCell}{minCountsPerCell}

\item{maxCountsPerCell}{maxCountsPerCell}

\item{minFeaturesPerCell}{minFeaturesPerCell}

\item{maxFeaturesPerCell}{maxFeaturesPerCell}

\item{maxPercent.mt}{maxPercent.mt}
}
\value{

}
\description{
quality control according to 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'CellsPerGene'
}
\examples{
\dontrun{
SeuratObj <- QCfun(SeuratObj, plot = F)
}

}