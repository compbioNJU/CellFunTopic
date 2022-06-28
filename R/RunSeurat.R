
#' basic single cell analysis with Seurat
#'
#' @param SeuratObj Seurat object
#' @param nPCs Dimensions of reduction to use as input
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#'
#' @importFrom dplyr group_by top_n
#' @return Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratObj <- RunSeurat(SeuratObj, nPCs = 10, resolution = 0.5)
#' }
#'
#'
#'
RunSeurat <- function(SeuratObj, nPCs = 10, resolution = 0.5) {

  suppressMessages(require("Seurat", quietly=TRUE))
  SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 2000)
  # Visualize  the 10 most highly variable genes
  top10 <- head(VariableFeatures(SeuratObj), 10)
  if (!dir.exists(paths = "./scFunMap_output/plots/RunSeurat")) {
    dir.create("./scFunMap_output/plots/RunSeurat", recursive = TRUE)
  }
  pdf(file = "./scFunMap_output/plots/RunSeurat/VariableFeaturePlot.pdf")
  plot1 <- VariableFeaturePlot(SeuratObj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot1 + plot2)
  dev.off()
  # Scaling the data
  all.genes <- rownames(SeuratObj)
  SeuratObj <- ScaleData(SeuratObj, features = all.genes)
  # Perform linear dimensional reduction
  SeuratObj <- RunPCA(SeuratObj, features = VariableFeatures(object = SeuratObj))
  pdf(file = "./scFunMap_output/plots/RunSeurat/RunPCA.pdf")
  print(VizDimLoadings(SeuratObj, dims = 1:2, reduction = "pca"))
  print(DimPlot(SeuratObj, reduction = "pca"))
  print(DimHeatmap(SeuratObj, dims = 1, cells = 500, balanced = TRUE))
  print(DimHeatmap(SeuratObj, dims = 1:15, cells = 500, balanced = TRUE))
  dev.off()

  # Determine the ‘dimensionality’ of the dataset
  # NOTE: The JackStraw procedure can take a long time for big datasets.
  # SeuratObj <- JackStraw(SeuratObj, num.replicate = 100)
  # SeuratObj <- ScoreJackStraw(SeuratObj, dims = 1:20)
  # JackStrawPlot(SeuratObj, dims = 1:15)
  pdf(file = "./scFunMap_output/plots/RunSeurat/ElbowPlot.pdf")
  print(ElbowPlot(SeuratObj))
  dev.off()

  # Cluster the cells
  SeuratObj <- FindNeighbors(SeuratObj, dims = 1:nPCs)
  SeuratObj <- FindClusters(SeuratObj, resolution = resolution)

  # Run non-linear dimensional reduction (UMAP/tSNE)
  # 3D-UMAP
  SeuratObj <- RunUMAP(SeuratObj, dims = 1:nPCs, n.components=3, verbose=FALSE)
  pdf(file = "./scFunMap_output/plots/RunSeurat/RunUMAP.pdf")
  print(DimPlot(SeuratObj, reduction = "umap"))
  # DimPlot(SeuratObj, reduction = "umap", label = TRUE)
  dev.off()
  # 3D-tSNE
  SeuratObj <- RunTSNE(SeuratObj, dims = 1:nPCs, dim.embed = 3)
  pdf(file = "./scFunMap_output/plots/RunSeurat/RunTSNE.pdf")
  print(DimPlot(SeuratObj, reduction = "tsne"))
  # DimPlot(SeuratObj, reduction = "tsne", label = TRUE)
  dev.off()

  # Finding differentially expressed features (cluster biomarkers)
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  SeuratObj.markers <- FindAllMarkers(SeuratObj, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.05, return.thresh=0.8)
  top2 <- SeuratObj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  pdf(file = "./scFunMap_output/plots/RunSeurat/top2Markers.pdf")
  print(VlnPlot(SeuratObj, features = top2$gene))
  print(FeaturePlot(SeuratObj, features = top2$gene))
  dev.off()
  top10 <- SeuratObj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  pdf(file = "./scFunMap_output/plots/RunSeurat/top10Markers.pdf")
  print(DoHeatmap(SeuratObj, features = top10$gene) + NoLegend())
  dev.off()
  slot(object = SeuratObj, name = 'misc')[["Allmarkers"]] <- SeuratObj.markers
  slot(object = SeuratObj, name = 'misc')[["top10markers"]] <- top10
  # save the object
  # saveRDS(SeuratObj.markers, file = "./scFunMap_output/Allmarkers.rds")
  # saveRDS(top10, file = "./scFunMap_output/top10markers.rds")

  return(SeuratObj)
}
