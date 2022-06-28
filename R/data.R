#' Mitochondrial genes of all species.
#'
#' A dataset containing mitochondrial genes of all species in NCBI.
#'
#' @format A data frame with 310200 rows and 4 variables:
#' \describe{
#'   \item{txid}{taxonomy id of species}
#'   \item{ENTREZID}{ENTREZID of genes}
#'   \item{SYMBOL}{symbol of genes}
#'   \item{chromosome}{MT}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz}
"txid_gene_MT"


#' GO terms and their corresponding level or ontology.
#'
#' A dataset containing all GO terms and their corresponding level or ontology.
#'
#' @format A data frame with 137699 rows and 3 variables:
#' \describe{
#'   \item{GO}{GO IDs}
#'   \item{ont}{ontology}
#'   \item{level}{corresponding level}
#' }
#' @source level information were obtained from \code{getGOLevel()} of \code{clusterProfiler} package.
"GO2level"

