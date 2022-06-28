#' Load data from different sources
#'
#' @param data Can be a directory like "/outs/filtered_gene_bc_matrices/" provided by 10X.
#' Or an expression matrix of single-cell RNA-seq (matrix or dgCMatrix).
#' Or an object from R/python packages such as scater, scran, Seurat, monocle.
#' Or H5AD files that AnnData uses.
#' @param type 10X, 10X_h5, expMatrix, SingleCellExperiment, Seurat, CellDataSet, AnnData, loom
#' @param meta.data Additional cell-level metadata to add to the Seurat object. Should be a data frame where the rows are cell names and the columns are additional metadata fields.
#' @param species Species name, such as Homo sapiens or Mus musculus.
#' @importFrom Seurat Read10X CreateSeuratObject as.Seurat ReadH5AD
#'
#' @return Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratObj <- CellFunMap::readData(data = "~/outs/filtered_gene_bc_matrices/GRCh38/", type = "10X", species = "Homo sapiens")
#' SeuratObj <- CellFunMap::readData(data = counts, type = 'expMatrix', species = "Homo sapiens")
#' }
#'
#'
readData <- function(data,
                     type = c("10X", "10X_h5", "expMatrix", "SingleCellExperiment", "Seurat", "CellDataSet", "AnnData", "loom"),
                     meta.data = NULL,
                     species) {
  type <- match.arg(type)

  if (type == "10X") {
    if (!dir.exists(paths = data)) {
      stop("Directory provided does not exist")
    }
    counts <- Read10X(data.dir = data)
    SeuratObj <- CreateSeuratObject(counts = counts, meta.data = meta.data)
    } else if (type == "10X_h5") {
      counts <- Read10X_h5(file = data)
      SeuratObj <- CreateSeuratObject(counts = counts, meta.data = meta.data)
    }  else if (type == "expMatrix") {
      if (!inherits(x = data, what = c('matrix', 'Matrix', 'dgCMatrix', 'data.frame'))) {
        stop("data format must be 'matrix' or 'Matrix' or 'dgCMatrix' or 'data.frame' when type == 'expMatrix'")
      }
      SeuratObj <- CreateSeuratObject(counts = data, meta.data = meta.data)
    } else if (type == "Seurat") {
      if (!inherits(x = data, what = 'Seurat')) {
        stop("Please check the 'type' argument, class(data) is not 'Seurat'")
      }
      SeuratObj <- data
    } else if (type == "AnnData") {
      # message("We use the Seurat function 'ReadH5AD' to read the data from the H5AD files that AnnData uses.")
      SeuratObj <- ReadH5AD(file = data)
    } else {
      if (!inherits(x = data, what = type)) {
        stop("Wrong 'type' argument")
      }
      SeuratObj <- as.Seurat(data)
    }
  slot(object = SeuratObj, name = 'misc')[["species"]]  <- species
  return(SeuratObj)
}





