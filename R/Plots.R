

#' Heatmap of GSEA result
#'
#' @param SeuratObj Seurat object
#' @param by which GSEA result to show, one of "GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"
#' @param pathwayIDs pathway IDs to show. Default:NULL. If not NULL, parameter \code{topPath} becomes invalid.
#' @param toshow which GSEA score to show, "-logFDR", "enrichmentScore", "NES", "pvalue", "p.adjust"
#' @param topPath number of top pathways of each cluster to show
#' @param colour color of heatmap, see \code{RColorBrewer::brewer.pal.info}
#' @param scale if the values should be centered and scaled in either the row direction or the column direction, or none.
#' Corresponding values are "row", "column" and "none"
#' @param fontsize_row fontsize for rownames
#' @param cluster_rows boolean values determining if rows should be clustered.
#' @param cluster_cols boolean values determining if columns should be clustered.
#'
#'
#' @return A Heatmap-class object.
#' @export
#'
#' @examples
#' \dontrun{
#' gseaHeatmap(SeuratObj, by = "GO", toshow = "-logFDR", topPath = 10, colour = "Greens")
#' }
#'
#'
gseaHeatmap <- function(SeuratObj, by = "GO", pathwayIDs = NULL, toshow = "-logFDR", topPath = 10,
                        colour = "Greens", scale = "none", fontsize_row = 10,
                        cluster_rows = TRUE, cluster_cols = TRUE) {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]] %>% dplyr::mutate(`-logFDR`=-log10(p.adjust))
  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
      stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
    }
    topath <- GSEAresult %>% dplyr::filter(ID %in% pathwayIDs)
  } else {
    # topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=topPath, wt=toshow)
    # topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = toshow, n = topPath, with_ties = F) # 对input$topPath参数没反应
    topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(toshow), .by_group = TRUE) %>% dplyr::slice_head(n = topPath)
  }

  mat <- reshape2::acast(GSEAresult, Description~cluster, value.var=toshow)
  mat[is.na(mat)] <- 0
  mat <- mat[unique(topath$Description), ]
  rownames(mat) <- ifelse(nchar(rownames(mat)) > 60, paste(strtrim(rownames(mat), 60), "..."), rownames(mat))
  # ph <- ggplotify::as.ggplot(pheatmap::pheatmap(mat, fontsize_row = fontsize_row, color = colorRampPalette(c('white', brewer.pal(n=7,name=colour)))(100),
  #                                    scale = scale, angle_col = 315 , silent = T,
  #                                    cluster_rows = cluster_rows, cluster_cols = cluster_cols))
  # return(ph)

  if (scale == "row") {
    mat <- t(scale(t(mat), center = T, scale=T))
  } else if (scale == "column") {
    mat <- scale(mat, center = T, scale=T)
  }
  col_fun = circlize::colorRamp2(seq(min(mat), max(mat), length.out = 8), c("white", brewer.pal(n = 7, name = colour)))
  ht <- ComplexHeatmap::Heatmap(mat, name = ifelse(toshow == "-logFDR", "-log10(p.adjust)", toshow),
                col = col_fun, column_names_rot = 90, row_names_gp = grid::gpar(fontsize = fontsize_row),
                cluster_rows = cluster_rows, cluster_columns = cluster_cols, rect_gp = grid::gpar(col = "grey60"), border = TRUE)
  # draw(ht)
  ht

}




#' circle plot of clusters
#'
#' Link width shows number of intersection of pathways between clusters, link color is
#' in concordance with cluster which has higher enrichment score on intersection pathways.
#'
#' @param SeuratObj Seurat object
#' @param by which GSEA result to use, one of "GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"
#' @param pvaluecutoff calculate link width with pathways lower than pvaluecutoff
#' @param pathwayIDs pathway IDs to calculate link width
#' @param link_threshold only show links whose intersection number bigger than this threshold
#'
#' @importFrom circlize chordDiagram circos.clear
#' @importFrom utils combn
#' @importFrom stats setNames
#'
#' @export
#'
#' @examples
#' \dontrun{
#' circleplot(SeuratObj, by = "GO", pvaluecutoff = 0.01, pathwayIDs = NULL)
#' }
#'
#'
circleplot <- function(SeuratObj, by = "GO", pvaluecutoff = 0.01, pathwayIDs = NULL,
                       link_threshold = 10) {
  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
      stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
    }
    GSEAresult %<>% dplyr::filter(ID %in% pathwayIDs)
  } else {
    # 阈值 pvalue <= 0.01
    GSEAresult %<>% dplyr::filter(pvalue <= pvaluecutoff)
  }

  if (length(unique(GSEAresult$cluster)) < 2) {
    stop(">=2 clusters are required to generate plot.")
  }

  df <- as.character(unique(GSEAresult$cluster)) %>% combn(2) %>% t %>% as.data.frame
  # 连线宽度为pvalue <= 0.01的GO交集数目
  df$width <- apply(df, 1, FUN = function(x){
    gg1 <- GSEAresult[GSEAresult$cluster == x[1], "ID"]
    gg2 <- GSEAresult[GSEAresult$cluster == x[2], "ID"]
    length(intersect(gg1,gg2))
  })
  df <- df[df$width > link_threshold, ]
  if (nrow(df) < 1) {
    stop("No intersection of pathways between clusters was shown under such parameter settings.")
  }
  # sector颜色
  sectorc <- setNames(scPalette2(length(unique(GSEAresult$cluster))), unique(GSEAresult$cluster))
  # link颜色与交集GO的平均富集分数高的cluster颜色一致
  linkc <- apply(df, 1, FUN = function(x){
    gg1 <- GSEAresult[GSEAresult$cluster == x[1], "ID"]
    gg2 <- GSEAresult[GSEAresult$cluster == x[2], "ID"]
    gg <- intersect(gg1,gg2)
    pp1 <- GSEAresult %>% dplyr::filter(cluster == x[1] & (ID %in% gg)) %>% dplyr::pull(p.adjust)
    pp2 <- GSEAresult %>% dplyr::filter(cluster == x[2] & (ID %in% gg)) %>% dplyr::pull(p.adjust)
    ifelse(mean(pp1) <= mean(pp2), sectorc[[x[1]]], sectorc[[x[2]]])
  })

  chordDiagram(df, grid.col = sectorc, col=linkc, transparency = 0.7)
  circos.clear()
}


#' circle plot of clusters
#'
#' Helps to infer the relationship between clusters. Link width shows Pearson correlation or
#' Jaccard coefficient between clusters, calculated with GSEA result. Node size indicates cell number of each cluster.
#' If \code{by = "GO"} and \code{pathwayIDs = NULL}, only GO terms of level 5-6 are used for calculation.
#'
#' @param SeuratObj Seurat object
#' @param by which GSEA result to use for calculation
#' @param pathwayIDs IDs of pathways to use for calculation
#' @param color.use used to color nodes, can be a named vector
#' @param weight.scale scale the width or not
#' @param label.edge label edges or not
#' @param edge.curved The degree of edge bending
#' @param shape shape of nodes, 'circle' by default
#' @param layout network layout, circle by default
#' @param margin margin of plot
#' @param vertex.size.cex node size
#' @param vertex.label.cex size of node label
#' @param vertex.label.color color of label of nodes, 'black' by default
#' @param arrow.width width of arrow
#' @param arrow.size size of arrow
#' @param edge.label.color color of label of edge, 'black' by default
#' @param edge.label.cex size of edge label
#' @param edge.max.width width of edge
#' @param vertex.label.dist distance between label and nodes
#' @param link_threshold only show links whose correlation/Jaccard-index bigger than this threshold
#'
#' @import igraph
#' @importFrom stats cor
#'
#' @rdname clustercorplot
#' @export
#'
#' @examples
#' \dontrun{
#' # show Pearson correlation between clusters
#' clustercorplot(SeuratObj, by = "GO")
#' }
#'
#'
clustercorplot <- function(SeuratObj, by = "GO", pathwayIDs = NULL, color.use = NULL,
                           weight.scale = TRUE,label.edge = FALSE,edge.curved=0.2,shape='circle',
                           layout = igraph::in_circle(),margin=0.1, vertex.size.cex=1, link_threshold=0.5,
                           vertex.label.cex=1.5,vertex.label.color = "black",
                           arrow.width=1,arrow.size = 0.2,edge.label.color='black',
                           edge.label.cex=0.5,edge.max.width=8,vertex.label.dist=2) {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  mat <- -log10(reshape2::acast(data = GSEAresult, formula = ID ~ cluster, value.var = 'p.adjust'))
  mat[is.na(mat)] <- 0

  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
      stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
    }
    topPath <- rownames(mat)[rownames(mat) %in% pathwayIDs]
  } else {
    # 用至少在1个cluster中p.adjust < 0.05的通路来计算相关系数
    topPath <- rownames(mat)[rowSums(mat > -log10(0.05)) > 0]
    if ((by == "GO") & is.null(pathwayIDs)) {
      # 层级level 5，6的GO
      GO_level_5_6 <- GO2level[GO2level$level %in% c(5,6), "GO"] %>% as.character
      topPath <- topPath[topPath %in% GO_level_5_6]
    }
  }
  mat <- mat[topPath, ]

  correlation <- cor(mat, method = "pearson")
  diag(correlation) <- 0
  correlation[correlation < link_threshold] <- 0

  cellnumbers <- as.data.frame(table(SeuratObj@active.ident))
  rownames(cellnumbers) <- cellnumbers[,1]

  ###################### circle network
  g <- graph_from_adjacency_matrix(correlation, mode = "undirected", weighted = T)
  coords <- layout_(g, layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = setNames(scPalette(length(V(g))), names(V(g)))
  } else if (is.null(attr(color.use, which = "names"))) {
    color.use = setNames(color.use, names(V(g)))
  }

  vertex.size <- cellnumbers[names(V(g)), 2]
  vertex.size <- (vertex.size/max(vertex.size)*15+5)*vertex.size.cex

  V(g)$size<-vertex.size
  V(g)$color<-color.use[names(V(g))]
  V(g)$frame.color <- color.use[names(V(g))]
  V(g)$label.color <- vertex.label.color
  V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    E(g)$label<-E(g)$weight
  }
  ww <- abs(E(g)$weight)  # 负相关则E(g)$weight为负值
  if (weight.scale == TRUE) {
    E(g)$width<-0.3+edge.max.width/(max(ww)-min(ww))*(ww-min(ww))
  }else{
    E(g)$width<-0.3+edge.max.width*ww
  }

  E(g)$arrow.width<-arrow.width
  E(g)$arrow.size<-arrow.size
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
  E(g)$color <- ifelse(E(g)$weight > 0, "#FF7F0099", "#377EB899")  # 正相关和负相关的边用不同颜色表示

  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(V(g)), direction=-1, start=0)
  label.dist <- vertex.size/max(vertex.size) + vertex.label.dist
  par(mai=c(0.1,0.1,0.1,0.1))
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs)

}


#'
#' @rdname clustercorplot
#' @export
#' @examples
#' \dontrun{
#' # show Jaccard coefficient between clusters
#' clustercorplot_jaccard(SeuratObj, by = "GO")
#' }
#'
#'
clustercorplot_jaccard <- function(SeuratObj, by = "GO", pathwayIDs = NULL, color.use = NULL,
                                   weight.scale = TRUE,label.edge = FALSE,edge.curved=0.2,shape='circle',
                                   layout=in_circle(),margin=0.1, vertex.size.cex=1, link_threshold=0.4,
                                   vertex.label.cex=1.5,vertex.label.color= "black",
                                   arrow.width=1,arrow.size = 0.2,edge.label.color='black',
                                   edge.label.cex=0.5,edge.max.width=8,vertex.label.dist=2) {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  mat <- -log10(reshape2::acast(data = GSEAresult, formula = ID ~ cluster, value.var = 'p.adjust'))
  mat[is.na(mat)] <- 0
  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
      stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
    }
    topPath <- rownames(mat)[rownames(mat) %in% pathwayIDs]
  } else {
    # 用至少在1个cluster中p.adjust < 0.05的通路来计算相关系数
    topPath <- rownames(mat)[rowSums(mat > -log10(0.05)) > 0]
    if ((by == "GO") & is.null(pathwayIDs)) {
      # 层级level 5，6的GO
      GO_level_5_6 <- GO2level[GO2level$level %in% c(5,6), "GO"] %>% as.character
      topPath <- topPath[topPath %in% GO_level_5_6]
    }
  }
  mat <- mat[topPath, ]

  # 计算jaccard系数
  mat <- mat < 0.05
  mat <- mat + 0  #二值矩阵
  ddd <- proxy::dist(t(mat), method = "Jaccard")  #Jaccard距离
  ddd <- 1-ddd  #Jaccard系数
  jaccardMatrix <- as.matrix(ddd)
  jaccardMatrix[jaccardMatrix < link_threshold] <- 0

  cellnumbers <- as.data.frame(table(SeuratObj@active.ident))
  rownames(cellnumbers) <- cellnumbers[,1]

  ###################### circle network
  g <- graph_from_adjacency_matrix(jaccardMatrix, mode = "undirected", weighted = T)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = setNames(scPalette(length(V(g))), names(V(g)))
  } else if (is.null(attr(color.use, which = "names"))) {
    color.use = setNames(color.use, names(V(g)))
  }

  vertex.size <- cellnumbers[names(V(g)), 2]
  vertex.size <- (vertex.size/max(vertex.size)*15+5)*vertex.size.cex

  V(g)$size<-vertex.size
  V(g)$color<-color.use[names(V(g))]
  V(g)$frame.color <- color.use[names(V(g))]
  V(g)$label.color <- vertex.label.color
  V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    E(g)$label<-E(g)$weight
  }
  if (weight.scale == TRUE) {
    E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
  }else{
    E(g)$width<-0.3+edge.max.width*E(g)$weight
  }

  E(g)$arrow.width<-arrow.width
  E(g)$arrow.size<-arrow.size
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
  E(g)$color <- "grey"

  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(V(g)), direction=-1, start=0)
  label.dist <- vertex.size/max(vertex.size) + vertex.label.dist
  par(mai=c(0.1,0.1,0.1,0.1))
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs)
}




#' show relationship between clusters and pathways
#'
#' Hierarchical clustering of clusters and pathways according to GSEA result.
#' If \code{by = "GO"} and \code{pathwayIDs = NULL}, only GO terms of level 5-6 are used for calculation.
#'
#' @param SeuratObj Object of class "Seurat"
#' @param by which GSEA result to show
#' @param pathwayIDs IDs of pathways to show
#' @param topaths top n pathway of each cluster
#' @param cluster_cutree_k the clusters are divided into k, based on the hierarchical clustering (using cutree)
#' @param pathway_cutree_k the pathways are divided into k, based on the hierarchical clustering (using cutree)
#' @param color.use.cluster used to color the cluster nodes, a named vector or a vector
#' @param color.use.pathway used to color the pathway nodes
#' @param weight.scale scale the width or not
#' @param vertex.size.cex node size
#' @param vertex.label.cex size of node label
#' @param edge.max.width width of edge
#' @param vertex.label.color color of label of nodes, 'black' by default
#' @param alpha.edge transparency of edge color
#'
#' @importFrom graphics arrows layout par plot points rect text
#' @importFrom stats as.dendrogram cutree dist hclust
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' hierarchyplot_tree(SeuratObj, by = "GO", topaths = 10, cluster_cutree_k = 5, pathway_cutree_k = 20)
#' }
#'
#'
hierarchyplot_tree <- function(SeuratObj, by = "GO", pathwayIDs = NULL, topaths = 5,
                               cluster_cutree_k = NULL, pathway_cutree_k = NULL,
                               color.use.cluster = NULL, color.use.pathway = NULL, weight.scale = TRUE,
                               vertex.size.cex=1, vertex.label.cex=0.9, edge.max.width=3,
                               vertex.label.color= "black", alpha.edge = 0.6) {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]

  if ((by == "GO") & is.null(pathwayIDs)) {
    # GO of level 5，6
    GO_level_5_6 <- GO2level[GO2level$level %in% c(5,6), "GO"] %>% as.character
    GSEAresult %<>% dplyr::filter(ID %in% GO_level_5_6)
  }

  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
      stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
    }
    topath <- GSEAresult %>% dplyr::filter(ID %in% pathwayIDs)
  } else {
    # topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% dplyr::slice_head(n = topaths)
    # topath <- plyr::ddply(.data = GSEAresult, .variables = .(cluster), .fun = function(df){
    #   df %>% dplyr::arrange(p.adjust) %>% head(topaths)
    # })
    topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::slice_min(order_by = p.adjust, n = topaths, with_ties = F)
  }

  # 分别对cluster和pathway层次聚类，cutree后得到排序
  mat <- -log10(reshape2::acast(data = GSEAresult, formula = Description ~ cluster, value.var = 'p.adjust'))
  mat[is.na(mat)] <- 0
  mat <- mat[unique(topath$Description), ]
  # mat <- scale(mat)
  ###### 对cluster层次聚类
  out.dist=dist(t(mat), method="euclidean")
  out.hclust=hclust(out.dist,method="average")
  hc1 <- out.hclust
  # 计算最佳簇数
  if (is.null(cluster_cutree_k)) {
    ff <- factoextra::fviz_nbclust(t(mat), hcut, method = "silhouette", k.max = (nrow(t(mat))-1))
    cluster_cutree_k <- setNames(ff$data[,2],ff$data[,1]) %>% which.max %>% as.numeric
  }
  out.id <- cutree(out.hclust,k=cluster_cutree_k)
  tt1 <- table(out.id)[unique(out.id[out.hclust$order])]
  cluster_sorted <- out.hclust$labels[out.hclust$order]  # 相当于rownames(t(mat))[out.hclust$order]
  ###### 对pathway层次聚类
  out.dist=dist(mat, method="euclidean")
  out.hclust=hclust(out.dist,method="average")
  hc2 <- out.hclust
  # 计算最佳簇数
  if (nrow(mat) < 3) {
    pathway_cutree_k <- nrow(mat)
    message("pathway_cutree_k is set to ",nrow(mat)," for number of pathways < 3.")
  } else if (is.null(pathway_cutree_k)) {
    ff <- factoextra::fviz_nbclust(mat, hcut, method = "silhouette", k.max = (nrow(mat)-1))
    pathway_cutree_k <- setNames(ff$data[,2],ff$data[,1]) %>% which.max %>% as.numeric
  }
  out.id <- cutree(out.hclust,k=pathway_cutree_k)
  tt2 <- table(out.id)[unique(out.id[out.hclust$order])]
  pathway_sorted <- out.hclust$labels[out.hclust$order]

  # 边的信息
  edgess <- topath %>% dplyr::mutate(FDR = -log10(p.adjust)) %>% dplyr::select(cluster, Description, FDR)
  # 节点的信息
  cellnumbers <- as.data.frame(table(SeuratObj@active.ident))
  c_nodes <- setNames(cellnumbers, c("node", "size"))
  c_nodes$size <- as.numeric(c_nodes$size)
  c_nodes$size <- (c_nodes$size/max(c_nodes$size)*5+3)*vertex.size.cex
  p_nodes <- data.frame(node = unique(topath$Description), size = min(c_nodes$size)/2)
  c_nodes <- c_nodes[match(cluster_sorted, c_nodes$node), ]
  p_nodes <- p_nodes[match(pathway_sorted, p_nodes$node), ]

  # 节点和边的颜色
  if (is.null(color.use.cluster)) {
    color.use.cluster <- setNames(scPalette(length(cluster_sorted)), cluster_sorted)
  }else if (is.null(attr(color.use.cluster, which = "names"))) {
    color.use.cluster = setNames(color.use.cluster, cluster_sorted)
  } else {
    color.use.cluster <- color.use.cluster[as.character(cluster_sorted)]
  }
  if (is.null(color.use.pathway)) {
    color.use.pathway = scPalette2(length(pathway_sorted))
  }
  # edgess$color <- color.use.cluster[match(edgess$cluster, cluster_sorted)]
  edgess$color <- color.use.cluster[as.character(edgess$cluster)]
  edgess$color <- grDevices::adjustcolor(edgess$color,alpha.edge)
  # 边的宽度
  weight <- edgess$FDR
  if (weight.scale == TRUE) {
    edgess$width<-0.3+edge.max.width/(max(weight)-min(weight))*(weight-min(weight))
  }else{
    edgess$width<-0.3+edge.max.width*weight
  }

  l1 <- length(hc1$labels)
  l2 <- length(hc2$labels)
  layout(matrix(1:5,nrow=1), widths = c(3,2,8,4,4), heights = rep(15,5), respect = T)
  # The first dendrogram:
  op <- par(mar=c(3,3,3,0))
  on.exit(par(op))
  plot(as.dendrogram(hc1), horiz=TRUE, cex=1, leaflab="none", ylim=c(1,l1), xaxt = "n", yaxt = "n")
  # The first serie of labels:
  par(op)
  op <- par(mar=c(3,0,3,0))
  plot(NA, bty="n", axes=FALSE, xlim=c(0,1), ylim=c(1,l1), ylab="", xlab="")
  # 左边的矩形
  rr <- c(0, Reduce(f = sum, x = tt1, accumulate = TRUE))
  rectcolor = grDevices::adjustcolor(scPalette(length(tt1)),0.3)
  rect(xleft = 0, xright = 1,
       ybottom = rr[-length(rr)]+0.6, ytop = rr[-1]+0.4, col = rectcolor, border = rectcolor)
  text(x=1, y=1:l1, labels=cluster_sorted, pos=2, cex=vertex.label.cex+0.3, col = vertex.label.color)
  # nodes and arrows:
  par(op)
  op <- par(mar=c(3,0,3,0))
  plot(NA, bty="n", axes=FALSE, xlim=c(0,1), ylim=c(1,l1),ylab="", xlab="")
  points(rep(0.05, l1), 1:l1, pch=19, col = color.use.cluster, cex = c_nodes$size)
  points(rep(0.98, l2), seq(1,l1, length.out = l2), pch=19, col = color.use.pathway, cex = p_nodes$size)
  ord_arrow <- cbind(match(edgess$cluster, cluster_sorted), seq(1,l1, length.out = l2)[match(edgess$Description, pathway_sorted)])
  arrows(0.05, ord_arrow[,1], 0.96, ord_arrow[,2], code=2, length=0.04, lwd = edgess$width, col = edgess$color)
  # The second serie of labels:
  par(op)
  op <- par(mar=c(3,0,3,0))
  plot(NA, bty="n", axes=FALSE, xlim=c(0,1), ylim=c(1,l2), ylab="", xlab="")
  # 右边的矩形
  rr <- c(0, Reduce(f = sum, x = tt2, accumulate = TRUE))
  rectcolor = grDevices::adjustcolor(scPalette(length(tt2)),0.3)
  rect(xleft = 0, xright = 1.1,
       ybottom = rr[-length(rr)]+0.6, ytop = rr[-1]+0.4, col = rectcolor, border = rectcolor)
  # ptext <- pathway_sorted   # 换行
  # ptext[nchar(ptext)>30] <- paste0(substr(ptext[nchar(ptext)>30],1,30), "\n", substring(ptext[nchar(ptext)>30],31))
  # text(x=0, y=1:l2, labels=ptext, pos=4, cex=vertex.label.cex, col = vertex.label.color)
  text(x=0, y=1:l2, labels=pathway_sorted, pos=4, cex=vertex.label.cex, col = vertex.label.color)
  # And the second dendrogram (to reverse it I reversed the xlim vector:
  par(op)
  op <- par(mar=c(3,0,3,3))
  plot(as.dendrogram(hc2), horiz=TRUE, cex=1, xlim=c(0,max(hc2$height)), leaflab="none", ylim=c(1,l2), xaxt = "n", yaxt = "n")
}


# 用cowplot层层叠加的方法，适用于添加小直方图
# embeddedplot <- function(SeuratObj,
#                          by = "GO",
#                          pathwayIDs = NULL,
#                          topaths = 1,
#                          reduction = "umap",
#                          type = "hist") {
#   by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
#   GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
#   if ((by == "GO") & is.null(pathwayIDs)) {
#     # GO of level 5，6
#     GO_level_5_6 <- GO2level[GO2level$level %in% c(5,6), "GO"] %>% as.character
#     GSEAresult %<>% dplyr::filter(ID %in% GO_level_5_6)
#   }
#
#   if (!is.null(pathwayIDs)) {
#     if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
#       stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
#     }
#     topath <- GSEAresult %>% dplyr::filter(ID %in% pathwayIDs)
#   } else {
#     # topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% dplyr::slice_head(n = topaths)
#     topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::slice_min(order_by = p.adjust, n = topaths, with_ties = F)
#   }
#   # colors
#   if (is.null(pathwayIDs) & (topaths == 1)) {
#     topath %<>% dplyr::select(ID, cluster) %>% unique %>% dplyr::ungroup() %>% dplyr::mutate(color = scPalette(length(unique(cluster))))
#     histcols <- setNames(object = topath$color, topath$ID)
#     cols <- setNames(object = topath$color, topath$cluster)
#     # cols <- topath %>% dplyr::pull(color, name=cluster)  # code of same effect
#   } else {
#     histcols <- setNames(object = scPalette(length(unique(topath$ID))), unique(topath$ID))
#     cols <- setNames(object = scPalette(length(unique(SeuratObj@active.ident))), unique(SeuratObj@active.ident))
#   }
#
#   cell_eb <- data.frame(SeuratObj@reductions[[reduction]]@cell.embeddings[,1:2], cluster=as.character(SeuratObj@active.ident))
#   colnames(cell_eb) <- c('x','y','cluster')
#   # scatter plot
#   pp <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(cluster)), alpha=0.5, size=0.5) +
#     scale_color_manual(values = cols) + theme_pubr() +
#     theme(legend.position="none") + labs(x = paste0(reduction, "_1"), y = paste0(reduction, "_2"))
#   # coordinates
#   coords <- cell_eb %>% group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), .groups = "keep") %>% as.data.frame
#   # coords <- plyr::ddply(cell_eb, .(cluster), .fun = function(dff){
#   #   data.frame(x=quantile(dff$x, probs = 0.5),
#   #              y=quantile(dff$y, probs = 0.5), cluster=unique(dff$cluster))
#   # })  # code of same effect
#   coords$x <- (coords$x-min(cell_eb$x))/(max(cell_eb$x)-min(cell_eb$x))
#   coords$y <- (coords$y-min(cell_eb$y))/(max(cell_eb$y)-min(cell_eb$y))
#   coords$y[which.max(coords$y)] <- max(coords$y)-0.1  # 最上和最下的直方图越过坐标轴了
#   coords$y[which.min(coords$y)] <- min(coords$y)+0.1
#   # child plots
#   df <- GSEAresult %>% dplyr::mutate(FDR=-log10(p.adjust)) %>%
#     dplyr::filter(ID %in% unique(topath$ID)) %>% dplyr::select(ID, cluster, FDR)
#   df <- reshape2::acast(data = df, formula = ID ~ cluster, value.var = 'FDR')
#   df[is.na(df)] <- 0
#   df %<>% reshape2::melt() %>% setNames(c("ID","cluster","FDR"))
#   type <- match.arg(type, choices = c("hist", "pie"))
#   p <- ggdraw() +draw_plot(pp,0,0,1,1)
#   if (type == "hist") {
#     # histogram
#     for (cl in unique(df$cluster)) {
#       dff <- df[df$cluster == cl, , drop = F]
#       hp <- ggplot(dff, aes(x=ID, y=FDR, fill = ID)) + geom_bar(stat="identity") + theme_pubr() + scale_fill_manual(values = histcols) +
#         theme(legend.position="none", plot.background=element_rect(I(0),linetype=0), panel.background=element_rect(I(0)),
#               panel.grid.major=element_line(colour=NA), panel.grid.minor=element_line(colour=NA),
#               axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.line.y = element_blank())
#       p <- p + draw_plot(hp,x = coords[coords$cluster == cl, "x"], y = coords[coords$cluster == cl, "y"], width = 0.1,height = 0.1)
#     }
#   } else {
#     # pie chart
#     for (cl in unique(df$cluster)) {
#       dff <- df[df$cluster == cl, , drop = F]
#       hp <- ggplot(dff, aes(x="", y=FDR, fill = ID)) + geom_bar(stat = "identity", width = 1, alpha = 0.7) +
#         coord_polar(theta = "y") + theme_pubr() + scale_fill_manual(values = histcols) +
#         theme(legend.position="none", plot.background=element_rect(I(0),linetype=0), panel.background=element_rect(I(0)),
#               panel.grid.major=element_line(colour=NA), panel.grid.minor=element_line(colour=NA),
#               axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_blank())
#       p <- p + draw_plot(hp,x = coords[coords$cluster == cl, "x"], y = coords[coords$cluster == cl, "y"], width = 0.1,height = 0.1)
#     }
#   }
#   ## add 2 legends
#   # scatter plot
#   p1 <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(cluster)), alpha=0.5, size=0.5) +
#     scale_color_manual(values = cols, breaks = names(cols)) + theme_bw() + # breaks参数 图例顺序
#     guides(color = guide_legend(title = "cluster", override.aes = list(size=3))) + # override.aes参数 图例的图标大小
#     theme(legend.title = element_text(size = 12))
#   # histogram
#   p2 <- ggplot(df, aes(x=ID, y=FDR, fill = ID)) + geom_bar(stat="identity") +
#     scale_fill_manual(values = histcols, breaks = names(histcols)) + theme_bw() +
#     guides(fill = guide_legend(title = "pathway", override.aes = list(size=3))) +
#     theme(legend.title = element_text(size = 12))
#   # extract legends and convert to ggplot object
#   l1 <- as_ggplot(get_legend(p1))
#   l2 <- as_ggplot(get_legend(p2))
#   # arrange the plot and legends
#   ppp <- ggarrange(p, l1,l2, ncol=3, widths = c(4,1,1), heights = c(3,3,3))
#   return(ppp)
# }


#' show embedded histogram or pie chart on UMAP/TSNE plot
#'
#' embedded histogram or pie chart show pathway GSEA score on each clusters.
#' If \code{by = "GO"} and \code{pathwayIDs = NULL}, only GO terms of level 5-6 are used for calculation.
#'
#' @param SeuratObj Object of class "Seurat"
#' @param by which GSEA result to show
#' @param pathwayIDs IDs of pathways to show, if provided, parameter \code{topaths} become invalid.
#' @param topaths number of top pathways of each cluster to show
#' @param reduction "umap", "tsne", "pca"
#' @param type "hist", "pie", type of embedded plots
#' @param pie.size.cex size of pie chart
#'
#' @importFrom magrittr `%<>%`
#' @importFrom cowplot ggdraw draw_plot
#' @importFrom scatterpie geom_scatterpie
#' @importFrom stats median setNames
#'
#' @export
#'
#' @examples
#' \dontrun{
#' embeddedplot(SeuratObj, type = "pie")
#' embeddedplot(SeuratObj, topaths = 3, reduction = "tsne", type = "hist")
#' }
#'
#'
embeddedplot <- function(SeuratObj,
                         by = "GO",
                         pathwayIDs = NULL,
                         topaths = 1,
                         reduction = "umap",
                         type = "pie",
                         pie.size.cex = 1) {
  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  if ((by == "GO") & is.null(pathwayIDs)) {
    # GO of level 5，6
    GO_level_5_6 <- GO2level[GO2level$level %in% c(5,6), "GO"] %>% as.character
    GSEAresult %<>% dplyr::filter(ID %in% GO_level_5_6)
  }

  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
      stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
    }
    topath <- GSEAresult %>% dplyr::filter(ID %in% pathwayIDs)
  } else {
    # topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% dplyr::slice_head(n = topaths)
    topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::slice_min(order_by = p.adjust, n = topaths, with_ties = F)
  }
  # colors
  # if (is.null(pathwayIDs) & (topaths == 1)) {
  #   topath %<>% dplyr::select(ID, cluster) %>% unique %>% dplyr::ungroup() %>% dplyr::mutate(color = scPalette(length(unique(cluster))))
  #   histcols <- setNames(object = topath$color, topath$ID) #可能有重复的ID
  #   cols <- setNames(object = topath$color, topath$cluster) #可能有的cluster没跑出来gsea结果
  #   # cols <- topath %>% dplyr::pull(color, name=cluster)  # code of same effect
  # } else {
    histcols <- setNames(object = scPalette(length(unique(topath$ID))), unique(topath$ID))
    cols <- setNames(object = scPalette(length(unique(SeuratObj@active.ident))), unique(SeuratObj@active.ident))
  # }

  cell_eb <- data.frame(SeuratObj@reductions[[reduction]]@cell.embeddings[,1:2], cluster=as.character(SeuratObj@active.ident))
  colnames(cell_eb) <- c('x','y','cluster')
  # scatter plot
  pp <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(cluster)), alpha=0.5, size=0.2) +
    scale_color_manual(values = cols) + ggpubr::theme_pubr() +
    theme(legend.position="none") + labs(x = paste0(reduction, "_1"), y = paste0(reduction, "_2"))

  type <- match.arg(type, choices = c("hist", "pie"))
  if (type == "hist") {
    #### histogram
    # coordinates
    coords <- cell_eb %>% group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), .groups = "keep") %>% as.data.frame
    # coords <- plyr::ddply(cell_eb, .(cluster), .fun = function(dff){
    #   data.frame(x=quantile(dff$x, probs = 0.5),
    #              y=quantile(dff$y, probs = 0.5), cluster=unique(dff$cluster))
    # })  # code of same effect
    coords$x <- (coords$x-min(cell_eb$x))/(max(cell_eb$x)-min(cell_eb$x))
    coords$y <- (coords$y-min(cell_eb$y))/(max(cell_eb$y)-min(cell_eb$y))
    coords$y[which.max(coords$y)] <- max(coords$y)-0.1  # 最上和最下的直方图越过坐标轴了
    coords$y[which.min(coords$y)] <- min(coords$y)+0.1
    # child plots
    df <- GSEAresult %>% dplyr::mutate(FDR=-log10(p.adjust)) %>%
      dplyr::filter(ID %in% unique(topath$ID)) %>% dplyr::select(ID, cluster, FDR)
    df <- reshape2::acast(data = df, formula = ID ~ cluster, value.var = 'FDR')
    df[is.na(df)] <- 0
    df %<>% reshape2::melt() %>% setNames(c("ID","cluster","FDR"))

    p <- ggdraw() +draw_plot(pp,0,0,1,1)
    for (cl in unique(df$cluster)) {
      dff <- df[df$cluster == cl, , drop = F]
      hp <- ggplot(dff, aes(x=ID, y=FDR, fill = ID)) + geom_bar(stat="identity") + ggpubr::theme_pubr() + scale_fill_manual(values = histcols) +
        theme(legend.position="none", plot.background=element_rect(I(0),linetype=0), panel.background=element_rect(I(0)),
              panel.grid.major=element_line(colour=NA), panel.grid.minor=element_line(colour=NA),
              axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.line.y = element_blank())
      p <- p + draw_plot(hp,x = coords[coords$cluster == cl, "x"], y = coords[coords$cluster == cl, "y"], width = 0.1,height = 0.1)
    }
  } else {
    # pie chart
    # coordinates
    coords <- cell_eb %>% group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), .groups = "keep") %>% as.data.frame
    # dataframe to draw pie
    df <- GSEAresult %>% dplyr::mutate(FDR=-log10(p.adjust)) %>%
      dplyr::filter(ID %in% unique(topath$ID)) %>% dplyr::select(ID, cluster, FDR)
    piedf <- reshape2::acast(df, cluster ~ ID, value.var = "FDR")
    piedf[is.na(piedf)] <- 0
    piedf <- cbind(piedf, coords[match(rownames(piedf), coords$cluster), -1])
    # radius饼图直径表示cluster的细胞数目
    cellnumbers <- dplyr::count(cell_eb, cluster)
    sizee <- cellnumbers[match(rownames(piedf), cellnumbers[,1]), 2]
    piedf$radius <- (sizee/max(sizee)*0.8+0.1)*pie.size.cex
    p <- pp + ggnewscale::new_scale_fill() + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=piedf,
                                                             cols=colnames(piedf)[1:(ncol(piedf)-3)],color=NA) +
      coord_equal() + scale_fill_manual(values = histcols)
  }
  ## add 2 legends
  # scatter plot
  p1 <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(cluster)), alpha=0.5, size=0.5) +
    scale_color_manual(values = cols, breaks = names(cols)) + theme_bw() + # breaks参数 图例顺序
    guides(color = guide_legend(title = "cluster", override.aes = list(size=3))) + # override.aes参数 图例的图标大小
    theme(legend.title = element_text(size = 12))
  # histogram
  p2 <- ggplot(df, aes(x=ID, y=FDR, fill = ID)) + geom_bar(stat="identity") +
    scale_fill_manual(values = histcols, breaks = names(histcols)) + theme_bw() +
    guides(fill = guide_legend(title = "pathway", override.aes = list(size=3))) +
    theme(legend.title = element_text(size = 12))
  # extract legends and convert to ggplot object
  l1 <- ggpubr::as_ggplot(ggpubr::get_legend(p1))
  l2 <- ggpubr::as_ggplot(ggpubr::get_legend(p2))
  # arrange the plot and legends
  ppp <- ggpubr::ggarrange(p, l1,l2, ncol=3, widths = c(4,1,1), heights = c(3,3,3))
  return(ppp)
}




#' Scatter plot showing pathway enrichment score
#'
#' @param SeuratObj Object of class "Seurat"
#' @param by which GSEA result to show
#' @param pathwayID ID of pathway to show
#' @param reduction "umap", "tsne", "pca"
#' @param colour see \code{RColorBrewer::brewer.pal.info}
#' @param pointsize size of point
#' @param label label cluster or not
#' @param label.size size of label
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' pathwayScatterplot(SeuratObj, by = "GO", pathwayID = "GO:0002576")
#' }
#'
pathwayScatterplot <- function(SeuratObj, by = "GO", pathwayID = NULL, reduction = "umap",
                               colour = "OrRd", pointsize = 1, label = TRUE, label.size = 4) {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  gd <- GSEAresult %>% dplyr::mutate(Score=-log10(p.adjust)) %>% reshape2::acast(formula = ID ~ cluster, value.var = 'Score')
  gd[is.na(gd)] <- 0
  if (!pathwayID %in% rownames(gd)) {
    stop(pathwayID, " is not present in GSEA result of this Seurat object.")
  }
  gd <- gd[pathwayID, ][as.character(Seurat::Idents(SeuratObj))]
  # > all.equal(rownames(SeuratObj@meta.data), names(Seurat::Idents(SeuratObj)))
  # [1] TRUE
  # SeuratObj@meta.data[["pathwayScore"]] <- unname(gd)
  names(gd) <- names(Seurat::Idents(SeuratObj))
  SeuratObj@meta.data[["pathwayScore"]] <- gd[rownames(SeuratObj@meta.data)]

  # p <- Seurat::FeaturePlot(SeuratObj, features = "pathwayScore", pt.size = pointsize, reduction = reduction, label = label) +
  #   scale_color_gradientn(colours = brewer.pal(n=7,name=colour)) +
  #   ggtitle(sprintf("%s: %s",pathwayID, GSEAresult[GSEAresult$ID == pathwayID, "Description"][1]))
  suppressMessages(p <- Seurat::FeaturePlot(SeuratObj, features = "pathwayScore", pt.size = pointsize,
                                            reduction = reduction, label = label, label.size = label.size) +
                     scale_color_gradientn(name="-log10(p.adjust)", colours = brewer.pal(n=7,name=colour)) +
                     ggtitle(sprintf("%s: %s",pathwayID, GSEAresult[GSEAresult$ID == pathwayID, "Description"][1])))
  return(p)
}

# pathwayScatterplot <- function(SeuratObj, by = "GO", pathwayID = NULL, reduction = "umap", colour = "OrRd", pointsize = 1) {
#
#   by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
#   GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
#
#   df <- data.frame(SeuratObj@reductions[[reduction]]@cell.embeddings[,1:2], cluster=as.character(SeuratObj@active.ident))
#   colnames(df) <- c('x','y','cluster')
#   gd <- GSEAresult %>% dplyr::mutate(Score=-log10(p.adjust)) %>% dplyr::filter(ID == pathwayID) %>% dplyr::select(cluster, Score)
#   rownames(gd) <- gd$cluster
#   df <- df %>% mutate(Score = gd[as.character(cluster),'Score'])
#
#   p <- ggplot(data=df, aes(x,y)) + geom_point(aes(colour=Score), alpha=1, size=pointsize) +
#     scale_color_gradientn(colours = brewer.pal(n=7,name=colour)) + theme_pubr() +
#     ggtitle(sprintf("%s: %s",pathwayID, GSEAresult[GSEAresult$ID == pathwayID, "Description"][1])) +
#     theme(legend.position="right", plot.title = element_text(hjust = 0.5, size = 20)) +
#     labs(x = paste0(reduction, "_1"), y = paste0(reduction, "_2"))
#   return(p)
# }



#' boxplot showing enrichment score of child or parent GOs of specific GO
#'
#' @param SeuratObj Object of class "Seurat"
#' @param goid specific GO ID
#' @param type "child","parent", show child or parent GOs of specific GO
#' @param pointsize size of point
#' @param flip whether to flip the coordinates
#'
#' @export
#'
#' @examples
#' \dontrun{
#' GOboxplot(SeuratObj, goid = "GO:0002576", type = "child")
#' GOboxplot(SeuratObj, goid = "GO:0002576", type = "parent", flip = T)
#' }
#'
#'
GOboxplot <- function(SeuratObj, goid = NULL, type = "child", pointsize = 1, flip = FALSE) {
  type <- match.arg(type, choices = c("child", "parent"))
  if (type == "child") {
    nodes <- GOfuncR::get_child_nodes(goid)$child_go_id
  } else {
    nodes <- GOfuncR::get_parent_nodes(goid)$parent_go_id
  }
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[["GSEAresult_GO"]] %>% dplyr::mutate(Score=-log10(p.adjust))
  godata <- GSEAresult[GSEAresult$ID %in% nodes, ] %>% dplyr::filter(!is.na(p.adjust)) %>% dplyr::select(cluster, Score)
  p <- ggpubr::ggboxplot(godata, x = "cluster", y = "Score", color = "cluster", size = pointsize,
                 palette = scPalette(length(unique(godata$cluster))), add = "jitter",
                 ylab = "-log10(p.adjust)") +
    # title = sprintf("%s: %s",goid,get_names(goid)$go_name), ylab = "-log10(p.adjust)") +
    theme(legend.position="right",
          plot.title = element_text(hjust = 0.5, size = 20)
    )
  if (flip) {
    return(p + coord_flip())
  }
  return(p)
}




#' Network showing cosine similarity between clusters according to GSEA result
#'
#' @param SeuratObj Seurat object
#' @param by GO KEGG Reactome MSigDb WikiPathways DO NCG DGN.
#' @param layout layout_nicely, layout_with_fr, etc. Either a function or a numeric matrix.
#' It specifies how the vertices will be placed on the plot.
#' @param cos_sim_thresh only draw links of cosine similarity bigger than this
#' @param p.adjust_thresh threshold of p.adjust to filter pathways used to calculate cosine similarity
#' @param SEED seed
#' @param node.cex node size
#' @param width_range range of width of links
#' @param text_size size of text
#' @param vertex.label.dist The distance of the label from the center of the vertex.
#'
#' @import igraph
#' @importFrom widyr pairwise_similarity
#' @importFrom magrittr set_colnames `%<>%`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Cosine_networkByGSEA(SeuratObj, layout=layout_nicely, cos_sim_thresh=0.6,
#'                      p.adjust_thresh=0.05, SEED=123, node.cex=3, width_range=c(0.1, 0.8))
#' Cosine_networkByGSEA(SeuratObj, layout=layout_with_fr, cos_sim_thresh=0.8,
#'                      p.adjust_thresh=0.05, SEED=123, node.cex=5, width_range=c(0.8, 4))
#' }
#'
#'
Cosine_networkByGSEA <- function(SeuratObj,
                                 by = "GO",
                                 layout=layout_nicely,
                                 cos_sim_thresh=0.6,
                                 p.adjust_thresh=0.05,
                                 SEED=123,
                                 node.cex=5,
                                 width_range=c(0.8, 4),
                                 text_size = 1,
                                 vertex.label.dist=0.5) {
  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  GSEAresult %<>% dplyr::mutate(logFDR=-log2(p.adjust)) # -log2 instead of -log10 to enlarge the difference
  topath <- GSEAresult %>% dplyr::filter(p.adjust <= p.adjust_thresh) %>% dplyr::distinct(ID) %>% dplyr::pull(ID) # only use part of pathways to calculate cosine similarity
  if (length(topath) == 0) {
    stop("Please adjust parameter 'p.adjust_thresh', no pathways left under such filtering")
  }
  dd <- GSEAresult %>% dplyr::filter(ID %in% topath) %>% dplyr::select(cluster, ID, logFDR)
  dd$cluster <- as.character(dd$cluster)
  ccc <- dd %>% widyr::pairwise_similarity(item = cluster, feature = ID, value = logFDR, upper=F)
  links <- ccc %>% dplyr::filter(similarity >= cos_sim_thresh) %>% magrittr::set_colnames(c('source', 'target', 'weight')) %>%
    dplyr::mutate(width=weight)
  if (nrow(links) == 0) {
    stop("Please adjust parameter 'cos_sim_thresh', no links left under such filtering")
  }
  links$width <- scales::rescale(links$width, width_range)
  nodes <- data.frame(name=unique(c(links$source, links$target)), stringsAsFactors = F)
  nodes$color <- scPalette(nrow(nodes))
  gg <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  V(gg)$frame.color <- NA
  V(gg)$size <- node.cex
  V(gg)$label.cex <- text_size
  V(gg)$label.color <- "black"
  V(gg)$label.font <- 2
  set.seed(SEED)
  plot(gg, layout=layout, vertex.label.dist=vertex.label.dist)
}


#' Hierarchical edge bundling plots helps visualizing correlation or similarity between clusters
#'
#' @param SeuratObj Seurat object
#' @param node.by draw each node as \code{node.by}, column of \code{SeuratObj@meta.data}
#' @param group.by color each node by \code{group.by}, column of \code{SeuratObj@meta.data}
#' @param by GO KEGG Reactome MSigDb WikiPathways DO NCG DGN.
#' @param link_threshold only show links of similarity or correlation above threshold
#' @param link_width width of links
#' @param p.adjust_thresh threshold of p.adjust to filter pathways used to calculate cosine similarity/correlation
#' @param method one of \code{c('cosine similarity', 'pearson', 'spearman')}
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom plyr ddply .
#' @importFrom ggraph ggraph get_con geom_conn_bundle scale_edge_color_distiller geom_node_text geom_node_point guide_edge_colourbar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' edge_bundling_GSEA(SeuratObj, link_threshold=0.6, p.adjust_thresh=0.05,
#'              method='cosine similarity', node.by='cluster', group.by='cellType')
#' edge_bundling_GSEA(SeuratObj, link_threshold=0.6, p.adjust_thresh=0.05,
#'              method='pearson', node.by='cluster', group.by='cellType')
#' }
#'
edge_bundling_GSEA <- function(SeuratObj, by = "GO", link_threshold=0.6, link_width=0.9, p.adjust_thresh=0.05,
                               method="cosine similarity", node.by="cluster", group.by="cellType") {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  method <- match.arg(method, choices = c('cosine similarity', 'pearson', 'spearman'))

  majorCellType <- unique(SeuratObj@meta.data[, c(group.by, node.by)]) %>% dplyr::arrange(get(group.by))
  ## define data.frame with hierarchical information
  d1 <- data.frame(from="origin", to=unique(majorCellType[, 1]), stringsAsFactors = F)
  d2 <- majorCellType %>% magrittr::set_colnames(c('from', 'to'))
  hierarchy <- rbind(d1, d2)
  ## node information, cell type as group, cell number of each cluster as node size
  vertices <- data.frame(name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))), stringsAsFactors = F)
  vertices$group <- hierarchy$from[match(vertices$name, hierarchy$to)]
  nn <- dplyr::count(SeuratObj@meta.data, get(node.by)) %>% magrittr::set_colnames(c('cluster', 'n'))
  vertices$cellnumber <- nn$n[match(vertices$name, nn$cluster)] # cell number of each cluster as node size
  # calculate angle of leaves' labels
  vertices$id <- NA
  myleaves <- which(is.na(match(vertices$name, hierarchy$from))) # only focus on leaves
  nleaves <- length(myleaves)
  vertices$id[myleaves] <- seq(1:nleaves)
  vertices$angle <- 120 - 360 * vertices$id/nleaves
  vertices$hjust <- ifelse(vertices$angle < -90 | vertices$angle > 90, 0, 1)
  vertices$angle <- ifelse(vertices$angle < -90 | vertices$angle > 90, vertices$angle+180, vertices$angle) # flip angle BY to make them readable
  GSEAresult %<>% dplyr::mutate(logFDR=-log10(p.adjust))
  topath <- GSEAresult %>% dplyr::filter(p.adjust < p.adjust_thresh) %>% dplyr::distinct(ID) %>% dplyr::pull(ID)  # only use part of pathways to calculate
  if (length(topath) < 1) {
    stop("Please adjust parameter 'p.adjust_thresh', no pathways left under such filtering")
  }
  dd <- GSEAresult %>% dplyr::filter(ID %in% topath) %>% dplyr::select(cluster, ID, logFDR)
  if (method == 'cosine similarity') {
    ccc <- dd %>% widyr::pairwise_similarity(item = cluster, feature = ID, value = logFDR, upper=F)
  } else if (method == 'pearson') {
    ccc <- dd %>% widyr::pairwise_cor(item = cluster, feature = ID, value = logFDR, method ="pearson", upper=F)
  } else if (method == 'spearman') {
    ccc <- dd %>% widyr::pairwise_cor(item = cluster, feature = ID, value = logFDR, method ="spearman", upper=F)
  }
  colnames(ccc)[3] <- 'links'
  connect <- ccc %>% dplyr::filter(links >= link_threshold)
  if (nrow(connect) < 1) {
    stop("Please adjust parameter 'link_threshold', no links left under such filtering")
  }

  from <- match(connect$item1, vertices$name)
  to <- match(connect$item2, vertices$name)
  mygraph <- graph_from_data_frame(hierarchy, vertices=vertices)
  pp <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE)
  df4 <- get_con(from = from, to = to)(pp$data)
  bb <- plyr::ddply(df4, .(con.id), .fun = function(df){
    df$links <- connect$links[connect$item1 == df$name[1] & connect$item2 == df$name[nrow(df)]]
    df
  })
  pp <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
    geom_conn_bundle(data = bb, alpha=1, width=link_width, aes(colour=links)) +
    scale_edge_color_distiller(palette = "RdPu", direction = 1) + # BuPu
    geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=name, angle = angle, hjust=hjust, colour=group), size=3, alpha=1) +
    geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=cellnumber), alpha=1) +
    scale_colour_manual(values= scPalette2(length(unique(majorCellType[, 1]))), guide=guide_legend(override.aes = list(size=5))) +
    scale_size_continuous(range = c(2,15)) +
    coord_fixed() +
    theme_void() +
    theme(plot.margin=unit(c(0,0,0,0),"cm")) +
    expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6))
  return(pp)
}


#' show pathways and genes in chord diagram
#'
#' @param SeuratObj Seurat object
#' @param genes genes to draw
#'
#' @importFrom circlize circos.par chordDiagram circos.trackPlotRegion get.cell.meta.data circlize circos.text circos.clear
#'
#' @export
#'
#' @examples
#' \dontrun{
#' GOcircleplot(SeuratObj, genes=NULL)
#' }
#'
GOcircleplot <- function(SeuratObj, genes=NULL) {
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[["GSEAresult_GO"]]
  # only show top 1 pathway of each cluster
  toppws <- GSEAresult %>% dplyr::group_by(cluster) %>%
    dplyr::slice_min(order_by = p.adjust, n = 1, with_ties = F) %>% dplyr::pull(ID) %>% unique
  if (is.null(genes)) {
    genes <- sample(Seurat::VariableFeatures(SeuratObj), 50)
  }
  species <- slot(object = SeuratObj, name = 'misc')[["species"]]
  OrgDb <- getOrgDb(species)
  df <- AnnotationDbi::select(OrgDb, keys=toppws, columns = "SYMBOL", keytype="GOALL")
  dd <- unique(df[, c('GOALL', 'SYMBOL')])
  dd <- dd %>% dplyr::filter(SYMBOL %in% genes) %>% magrittr::set_colnames(c("from", "to"))
  pwscols <- setNames(scPalette2(length(unique(dd$from))), unique(dd$from))
  dd$color <- pwscols[dd$from]

  sectorc <- c(pwscols, setNames(rep("grey", length(unique(dd$to))), unique(dd$to)))
  par(mar = rep(4, 4))
  circos.par(track.margin = c(-0.08, 0.1), points.overflow.warning = FALSE)
  chordDiagram(dd, grid.col = sectorc, col=dd$color, transparency = 0.4, directional = 1,
               direction.type = c("arrows", "diffHeight"), diffHeight  = -0.03,
               link.arr.type = "big.arrow", annotationTrack = "grid")
  # Add text and axis
  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {

      xlim = get.cell.meta.data("xlim")
      sector.index = get.cell.meta.data("sector.index")

      #text direction (dd) and adjusmtents (aa)
      theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
      dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
      aa = c(1, 0.5)
      if(theta < 90 || theta > 270)  aa = c(0, 0.5)
      circos.text(x=mean(xlim), y=1.7, labels=sector.index, facing = dd, cex=0.6,  adj = aa)
    }
  )
  circos.clear()
}


#' cluster functional terms into groups by clustering the similarity matrix of the terms
#'
#' We utilize package simplifyEnrichment to cluster GO terms into groups from the semantic similarity matrix. Summaries of GO terms in each cluster are visualized by word clouds.
#' https://bioconductor.org/packages/release/bioc/html/simplifyEnrichment.html
#' https://bioconductor.org/packages/release/bioc/vignettes/simplifyEnrichment/inst/doc/simplifyEnrichment.html
#' package "simplifyEnrichment" is developed based on R version 4.0
#'
#'
#' @param SeuratObj Seurat object
#' @param by GO, KEGG, Reactome, MSigDb
#' @param pathwayIDs a vector of pathway IDs, such as a vector of GO IDs
#' @param showCategory  number of top pathways of each cluster
#' @param GO_ont Gene Ontology
#'
#' @export
#'
#' @examples
#' \dontrun{
#' simplifyEnrichmentplot(SeuratObj, by = "GO", showCategory = 10)
#' }
#'
#'
simplifyEnrichmentplot <- function(SeuratObj, by = "GO", pathwayIDs = NULL, showCategory = 10, GO_ont = "BP") {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]

  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
      stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
    }
    pathwayIDs <- GSEAresult$ID[GSEAresult$ID %in% pathwayIDs]
  } else {
    # topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% dplyr::slice_head(n = showCategory)
    pathwayIDs <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::slice_min(order_by = p.adjust, n = showCategory, with_ties = F) %>% dplyr::pull(ID)
  }

  if (by == "GO") {
    mat <- simplifyEnrichment::GO_similarity(pathwayIDs, ont = GO_ont)
    simplifyEnrichment::simplifyGO(mat, verbose = FALSE)
  } else if (by == "KEGG") {
    mat <- simplifyEnrichment::term_similarity_from_KEGG(pathwayIDs)
    simplifyEnrichment::simplifyEnrichment(mat, verbose = FALSE)
  } else if (by == "Reactome") {
    mat <- simplifyEnrichment::term_similarity_from_Reactome(pathwayIDs)
    simplifyEnrichment::simplifyEnrichment(mat, verbose = FALSE)
  } else if (by == "MSigDb") {
    mat <- simplifyEnrichment::term_similarity_from_MSigDB(pathwayIDs)
    simplifyEnrichment::simplifyEnrichment(mat, verbose = FALSE)
  }
}
