% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plots.R
\name{gseaHeatmap}
\alias{gseaHeatmap}
\title{Heatmap of GSEA result}
\usage{
gseaHeatmap(
  SeuratObj,
  by = "GO",
  pathwayIDs = NULL,
  toshow = "-logFDR",
  topPath = 10,
  colour = "Greens",
  scale = "none",
  fontsize_row = 10,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)
}
\arguments{
\item{SeuratObj}{Seurat object}

\item{by}{which GSEA result to show, one of "GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"}

\item{pathwayIDs}{pathway IDs to show. Default:NULL. If not NULL, parameter \code{topPath} becomes invalid.}

\item{toshow}{which GSEA score to show, "-logFDR", "enrichmentScore", "NES", "pvalue", "p.adjust"}

\item{topPath}{number of top pathways of each cluster to show}

\item{colour}{color of heatmap, see \code{RColorBrewer::brewer.pal.info}}

\item{scale}{if the values should be centered and scaled in either the row direction or the column direction, or none.
Corresponding values are "row", "column" and "none"}

\item{fontsize_row}{fontsize for rownames}

\item{cluster_rows}{boolean values determining if rows should be clustered.}

\item{cluster_cols}{boolean values determining if columns should be clustered.}
}
\value{
A Heatmap-class object.
}
\description{
Heatmap of GSEA result
}
\examples{
\dontrun{
gseaHeatmap(SeuratObj, by = "GO", toshow = "-logFDR", topPath = 10, colour = "Greens")
}


}
