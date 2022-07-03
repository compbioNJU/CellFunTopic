
#' basic single cell analysis with Seurat
#'
#' @param SeuratObj Seurat object
#' @param nPCs Dimensions of reduction to use as input
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param plot draw plots or not. If TRUE, plots will be produced and saved in working directory.
#'
#' @importFrom dplyr group_by top_n
#' @return Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratObj <- RunSeurat(SeuratObj, nPCs = 10, resolution = 0.5, plot = FALSE)
#' }
#'
#'
#'
RunSeurat <- function(SeuratObj, nPCs = 10, resolution = 0.5, plot = FALSE) {

  suppressMessages(require("Seurat", quietly=TRUE))
  SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 2000)
  # Visualize  the 10 most highly variable genes
  top10 <- head(VariableFeatures(SeuratObj), 10)
  if (plot & !dir.exists(paths = "./CellFunMap_output/plots/RunSeurat")) {
    dir.create("./CellFunMap_output/plots/RunSeurat", recursive = TRUE)
  }
  if (plot) {
    pdf(file = "./CellFunMap_output/plots/RunSeurat/VariableFeaturePlot.pdf")
    plot1 <- VariableFeaturePlot(SeuratObj)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    print(plot1 + plot2)
    dev.off()
  }

  # Scaling the data
  all.genes <- rownames(SeuratObj)
  SeuratObj <- ScaleData(SeuratObj, features = all.genes)
  # Perform linear dimensional reduction
  SeuratObj <- RunPCA(SeuratObj, features = VariableFeatures(object = SeuratObj))
  if (plot) {
    pdf(file = "./CellFunMap_output/plots/RunSeurat/PCA.pdf")
    print(VizDimLoadings(SeuratObj, dims = 1:2, reduction = "pca"))
    print(DimPlot(SeuratObj, reduction = "pca"))
    print(DimHeatmap(SeuratObj, dims = 1, cells = 500, balanced = TRUE))
    print(DimHeatmap(SeuratObj, dims = 1:15, cells = 500, balanced = TRUE))
    dev.off()
  }

  # Determine the ‘dimensionality’ of the dataset
  # NOTE: The JackStraw procedure can take a long time for big datasets.
  # SeuratObj <- JackStraw(SeuratObj, num.replicate = 100)
  # SeuratObj <- ScoreJackStraw(SeuratObj, dims = 1:20)
  # JackStrawPlot(SeuratObj, dims = 1:15)
  if (plot) {
    pdf(file = "./CellFunMap_output/plots/RunSeurat/ElbowPlot.pdf")
    print(ElbowPlot(SeuratObj))
    dev.off()
  }

  # Cluster the cells
  SeuratObj <- FindNeighbors(SeuratObj, dims = 1:nPCs)
  SeuratObj <- FindClusters(SeuratObj, resolution = resolution)

  # Run non-linear dimensional reduction (UMAP/tSNE)
  SeuratObj <- RunUMAP(SeuratObj, dims = 1:nPCs)
  SeuratObj <- RunTSNE(SeuratObj, dims = 1:nPCs)

  if (plot) {
    pdf(file = "./CellFunMap_output/plots/RunSeurat/UMAP.pdf")
    print(DimPlot(SeuratObj, reduction = "umap"))
    dev.off()
    pdf(file = "./CellFunMap_output/plots/RunSeurat/TSNE.pdf")
    print(DimPlot(SeuratObj, reduction = "tsne"))
    dev.off()
  }

  # Finding differentially expressed features (cluster biomarkers)
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  SeuratObj.markers <- FindAllMarkers(SeuratObj, only.pos = TRUE, min.pct = 0.0001, logfc.threshold = 0.0001, return.thresh=0.9)
  top10 <- SeuratObj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  if (plot) {
    pdf(file = "./CellFunMap_output/plots/RunSeurat/top10Markers.pdf")
    print(DoHeatmap(SeuratObj, features = top10$gene) + NoLegend())
    dev.off()
  }

  slot(object = SeuratObj, name = 'misc')[["Allmarkers"]] <- SeuratObj.markers
  slot(object = SeuratObj, name = 'misc')[["top10markers"]] <- top10

  return(SeuratObj)
}
