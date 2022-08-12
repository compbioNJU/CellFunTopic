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


#' Mitochondrial genes of all species.
#'
#' A dataset containing mitochondrial genes of all species in NCBI.
#'
#' @format A data frame with 1082809 rows and 5 variables:
#' \describe{
#'   \item{txid}{taxonomy id of species}
#'   \item{species}{species name}
#'   \item{ENTREZID}{ENTREZID of genes}
#'   \item{SYMBOL}{symbol of genes}
#'   \item{chromosome}{MT}
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz}
"txid_species_gene_MT"


#' Mitochondrial genes of 3 species.
#'
#' A list containing mitochondrial genes of 3 species, Homo sapiens, Mus musculus and Arabidopsis thaliana.
#'
#' @format A list containing 3 elements:
#' \describe{
#'   \item{Homo sapiens}{mitochondrial genes of several gene ID types}
#'   \item{Mus musculus}{mitochondrial genes of several gene ID types}
#'   \item{Arabidopsis thaliana}{mitochondrial genes of several gene ID types}
#' }
#' @source mitochondrial genes information derived from bioMart and NCBI.
"MTGENE"


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


#' built-in reference topic model derived from HCL.
#'
#' built-in reference topic model derived from human cell landscape (HCL).
#' citation: Han, X., Zhou, Z., Fei, L. et al. Construction of a human cell landscape at single-cell level. Nature 581, 303–309 (2020).
#' (\url{https://doi.org/10.1038/s41586-020-2157-4})
#'
#' @format A LDA_VEM topic model with 63 topics.
#'
#' @source single-cell data comes from \url{https://figshare.com/articles/dataset/HCL_DGE_Data/7235471}
"HCL_ldaOut"


#' built-in reference topic model derived from MCA.
#'
#' built-in reference topic model derived from Mouse Cell Atlas (MCA).
#' citation: \url{https://doi.org/10.1016/j.cell.2018.02.001}.
#'
#' @format A LDA_VEM topic model with 52 topics.
#'
#' @source single-cell data comes from \url{https://figshare.com/articles/dataset/HCL_DGE_Data/7235471}
"MCA_ldaOut"







