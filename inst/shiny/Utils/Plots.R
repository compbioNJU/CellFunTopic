

gseaHeatmap <- function(SeuratObj, by = "GO", pathwayIDs = NULL, toshow = "-logFDR", topPath = 10,
                        colour = "Greens", scale = "none", fontsize_row = 10) {
  
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
    # topath <- plyr::ddply(.data = GSEAresult, .variables = .(cluster), .fun = function(df){
    #   df %>% dplyr::arrange(desc(toshow)) %>% head(topPath)
    # })
    topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::arrange(desc(toshow), .by_group = TRUE) %>% dplyr::slice_head(n = topPath)
  }
  
  mat <- reshape2::acast(GSEAresult, Description~cluster, value.var=toshow)
  mat[is.na(mat)] <- 0
  mat <- mat[unique(topath$Description), ]
  rownames(mat) <- ifelse(nchar(rownames(mat)) > 60, paste(strtrim(rownames(mat), 60), "..."), rownames(mat))
  # ph <- pheatmap::pheatmap(mat, fontsize_row = fontsize_row, color = colorRampPalette(c('white', brewer.pal(n=7,name=colour)))(100),
  #                          scale = scale, angle_col = 315)
  # ph
  ph <- as.ggplot(pheatmap::pheatmap(mat, fontsize_row = fontsize_row, color = colorRampPalette(c('white', brewer.pal(n=7,name=colour)))(100),
                                     scale = scale, angle_col = 315))
  return(ph)
}

circleplot <- function(SeuratObj, by = "GO", pvaluecutoff = 0.01, pathwayIDs = NULL) {
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
  df <- df[df$width>0, ]
  if (nrow(df) < 1) {
    stop("No intersection of pathways between clusters.")
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


clustercorplot <- function(SeuratObj, by = "GO", pathwayIDs = NULL, color.use = NULL, 
                           weight.scale = FALSE,label.edge = FALSE,edge.curved=0.2,shape='circle',
                           layout=in_circle(),margin=0.1, vertex.size.cex=1,
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
      GO2level <- readRDS("Data/GO2level.rds")
      GO_level_5_6 <- GO2level[GO2level$level %in% c(5,6), "GO"] %>% as.character
      topPath <- topPath[topPath %in% GO_level_5_6]
    }
  }
  mat <- mat[topPath, ]
  
  correlation <- cor(mat, method = "pearson")
  diag(correlation) <- 0
  
  cellnumbers <- as.data.frame(table(SeuratObj@active.ident))
  rownames(cellnumbers) <- cellnumbers[,1]
  
  ###################### circle network
  g <- graph_from_adjacency_matrix(correlation, mode = "undirected", weighted = T)
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


clustercorplot_jaccard <- function(SeuratObj, by = "GO", pathwayIDs = NULL, color.use = NULL, 
                                   weight.scale = FALSE,label.edge = FALSE,edge.curved=0.2,shape='circle',
                                   layout=in_circle(),margin=0.1, vertex.size.cex=1,
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
      GO2level <- readRDS("Data/GO2level.rds")
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




hierarchyplot_tree <- function(SeuratObj, by = "GO", pathwayIDs = NULL, topaths = 5, 
                               cluster_cutree_k = NULL, pathway_cutree_k = NULL,
                               color.use.cluster = NULL, color.use.pathway = NULL, weight.scale = TRUE, 
                               vertex.size.cex=1, vertex.label.cex=0.9, edge.max.width=3,
                               vertex.label.color= "black", alpha.edge = 0.6) {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  
  if ((by == "GO") & is.null(pathwayIDs)) {
    # GO of level 5，6
    GO2level <- readRDS("Data/GO2level.rds")
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
  text(x=1, y=1:l1, labels=cluster_sorted, pos=2, cex=vertex.label.cex+0.3, col = vertex.label.color)
  # 左边的矩形
  rr <- c(0, Reduce(f = sum, x = tt1, accumulate = TRUE))
  rectcolor = grDevices::adjustcolor(scPalette(length(tt1)),0.3)
  rect(xleft = 0, xright = 1, 
       ybottom = rr[-length(rr)]+0.6, ytop = rr[-1]+0.4, col = rectcolor, border = rectcolor)
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
  # ptext <- pathway_sorted   # 换行
  # ptext[nchar(ptext)>30] <- paste0(substr(ptext[nchar(ptext)>30],1,30), "\n", substring(ptext[nchar(ptext)>30],31))
  # text(x=0, y=1:l2, labels=ptext, pos=4, cex=vertex.label.cex, col = vertex.label.color)
  text(x=0, y=1:l2, labels=pathway_sorted, pos=4, cex=vertex.label.cex, col = vertex.label.color)
  # 右边的矩形
  rr <- c(0, Reduce(f = sum, x = tt2, accumulate = TRUE))
  rectcolor = grDevices::adjustcolor(scPalette(length(tt2)),0.3)
  rect(xleft = 0, xright = 1.1, 
       ybottom = rr[-length(rr)]+0.6, ytop = rr[-1]+0.4, col = rectcolor, border = rectcolor)
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
#     GO2level <- readRDS("Data/GO2level.rds")
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

# 添加小直方图用cowplot层层叠加的方法，添加小饼图用scatterpie包
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
    GO2level <- readRDS("Data/GO2level.rds")
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
    scale_color_manual(values = cols) + theme_pubr() +
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
      hp <- ggplot(dff, aes(x=ID, y=FDR, fill = ID)) + geom_bar(stat="identity") + theme_pubr() + scale_fill_manual(values = histcols) + 
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
  l1 <- as_ggplot(get_legend(p1))
  l2 <- as_ggplot(get_legend(p2))
  # arrange the plot and legends
  ppp <- ggarrange(p, l1,l2, ncol=3, widths = c(4,1,1), heights = c(3,3,3))
  return(ppp)
}



################ emapplotPie 

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x, y)))
}

graph_build <- function(y_union, geneSets, color, cex_line) {
  # Get graph.data.frame() result using y_union
  if (is.null(dim(y_union)) | nrow(y_union) == 1) { # when just one node
    g <- graph.empty(0, directed = FALSE)
    g <- add_vertices(g, nv = 1)
    V(g)$name <- y_union$Description
    V(g)$color <- "red"
  } else {
    id <- y_union$ID
    geneSets <- geneSets[id]
    n <- nrow(y_union)
    w <- matrix(NA, nrow = n, ncol = n)
    colnames(w) <- rownames(w) <- y_union$Description
    for (i in 1:n) {
      for (j in i:n) {
        w[i, j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
      }
    }
    wd <- reshape2::melt(w)
    wd <- wd[wd[, 1] != wd[, 2], ]
    wd <- wd[!is.na(wd[, 3]), ]
    g <- graph.data.frame(wd[, -3], directed = FALSE)
    E(g)$width = sqrt(wd[, 3] * 5) * cex_line
    g <- delete.edges(g, E(g)[wd[, 3] < 0.2])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == y_union$Description)))
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    colVar <- y_union[idx, color]
    V(g)$color <- colVar
  }
  return(g)
}

get_p <- function(y, g, y_union, cex_category, pie, layout){
  ## when y just have one line
  if(is.null(dim(y)) | nrow(y) == 1) {
    title <- y$cluster
    p <- ggraph(g) + geom_node_point(color="red", size=5 * cex_category) +
      geom_node_text(aes_(label=~name)) + theme_void() +
      labs(title=title)
    return(p)
  }
  # when there is just one unique ID(Description), just plot one pie
  if(is.null(dim(y_union)) | nrow(y_union) == 1) {
    p <- ggraph(g)
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)
    
    ID_Cluster_mat <- cbind(ID_Cluster_mat,1,1,0.1*cex_category)
    colnames(ID_Cluster_mat) <- c(colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
                                  "x", "y", "radius")
    
    
    p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                             cols=names(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
                             color=NA)+
      xlim(-3,3) + ylim(-3,3) + coord_equal()+
      geom_node_text(aes_(label=~name), repel=TRUE) +
      theme_void()+labs(fill = "cluster")
    return(p)
    
  }
  ggraph(g, layout=layout)
}

prepare_pie_category <- function(y, pie = "count") {
  pie <- match.arg(pie, c("count", "-log10FDR"))
  if (pie == "count") {
    y$count <- stringr::str_count(y$core_enrichment, pattern = "/") + 1
    pie_data <- y[,c("cluster", "Description", "count")]
    ID_Cluster_mat <- as.data.frame(reshape2::acast(pie_data, Description~cluster, value.var='count'))
    ID_Cluster_mat[is.na(ID_Cluster_mat)] <- 0
  } else if (pie == "-log10FDR") {
    pie_data <- y[,c("cluster", "Description", "p.adjust")]
    pie_data$p.adjust <- -log10(pie_data$p.adjust)
    ID_Cluster_mat <- as.data.frame(reshape2::acast(pie_data, Description~cluster, value.var='p.adjust'))
    ID_Cluster_mat[is.na(ID_Cluster_mat)] <- 0
  }
  return(ID_Cluster_mat)
}

emapplotPie <- function(SeuratObj, by = "GO", pathwayIDs = NULL, showCategory = 5, color = "p.adjust", layout = "kk",
                        node_label_cex = 1, node_size_cex = 1, pie = "count", cex_line = 1) {
  
  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  if ((by == "GO") & is.null(pathwayIDs)) {
    # GO of level 5，6
    GO2level <- readRDS("Data/GO2level.rds")
    GO_level_5_6 <- GO2level[GO2level$level %in% c(5,6), "GO"] %>% as.character
    GSEAresult %<>% dplyr::filter(ID %in% GO_level_5_6)
  }
  
  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
      stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
    }
    # y <- GSEAresult %>% dplyr::filter(ID %in% pathwayIDs)
    y <- subset(GSEAresult, ID %in% pathwayIDs)
  } else {
    # get top `showCategory` pathways of every cluster
    topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::slice_min(order_by = p.adjust, n = showCategory, with_ties = F)
    y <- subset(GSEAresult, ID %in% topath$ID)
  }
  # combine genes of the same ID (Description)
  rrr <- y[, c("ID", "Description", "p.adjust", "core_enrichment")]
  y_union <- plyr::ddply(.data = rrr, .variables = .(ID), .fun = function(x){
    ids <- unique(unlist(strsplit(x$core_enrichment, "/", fixed = TRUE)))
    cnt <- length(ids)
    ids <- paste0(ids, collapse="/")
    x <- x %>% dplyr::mutate(geneIDs = ids, count = cnt) %>% dplyr::select(-core_enrichment)
    x[!duplicated(x$ID), ,drop = FALSE]
  })
  
  y <- y[y$ID %in% y_union$ID, ]
  geneSets <- setNames(strsplit(as.character(y_union$geneIDs), "/",
                                fixed = TRUE), y_union$ID)
  
  g <- graph_build(y_union, geneSets=geneSets,color=color, cex_line=cex_line)
  
  p <- get_p(y = y, g = g, y_union = y_union, cex_category = node_size_cex,
             pie = pie, layout = layout)
  if (is.null(dim(y)) | nrow(y) == 1 | is.null(dim(y_union)) | nrow(y_union) == 1)
    return(p)
  
  if (length(E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                            colour='darkgrey')
  }
  
  ## then add the pie plot
  ## Get the matrix data for the pie plot
  ID_Cluster_mat <- prepare_pie_category(y,pie=pie)
  
  
  # plot the edge
  # get the X-coordinate and y-coordinate of pies
  aa <- p$data
  
  desc <- y_union$Description[match(rownames(ID_Cluster_mat),
                                    y_union$Description)]
  i <- match(desc, aa$name)
  
  ID_Cluster_mat$x <- aa$x[i]
  ID_Cluster_mat$y <- aa$y[i]
  
  #Change the radius value to fit the pie plot
  radius <- NULL
  ID_Cluster_mat$radius <- sqrt(aa$size[i] / sum(aa$size) * node_size_cex)
  #ID_Cluster_mat$radius <- sqrt(aa$size / pi)
  
  x_loc1 <- min(ID_Cluster_mat$x)
  y_loc1 <- min(ID_Cluster_mat$y)
  
  # if more than one cluster, plot piecharts
  if(ncol(ID_Cluster_mat) > 4) {
    p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                             cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +
      coord_equal()
    if (utils::packageVersion("ggrepel") >= "0.9.0") {
      p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
                              size = 3 * node_label_cex, bg.color = "white")
    } else {
      p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
                              size = 3 * node_label_cex)
    }
    p <- p + theme_void() +
      geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1,
                             n = 5,
                             labeller=function(x) round(sum(aa$size) * x^2 / node_size_cex)) +
      labs(fill = "cluster")
    return(p)
  }
  # else if just one cluster, plot nodes
  title <- colnames(ID_Cluster_mat)[1]
  p + geom_node_point(aes_(color=~color, size=~size))
  if (utils::packageVersion("ggrepel") >= "0.9.0") {
    p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
                            size = 3 * node_label_cex, bg.color = "white")
  } else {
    p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
                            size = 3 * node_label_cex)
  }
  p + theme_void() +
    scale_color_continuous(low="red", high="blue", name = color,
                           guide=guide_colorbar(reverse=TRUE)) +
    scale_size(range=c(3, 8) * node_size_cex)  +labs(title= title)
  
}



################ emapplot


geneInCategory <- function(x) {
  setNames(strsplit(as.character(x$core_enrichment), "/", fixed = TRUE), x$ID)
}

emapplot2 <- function(SeuratObj, cluster = NULL, by = "GO", pathwayIDs = NULL, showCategory = 20, color = "p.adjust", layout = "kk",
                      node_label_cex = 1, node_size_cex = 1, cex_line = 1) {
  
  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  if (is.null(cluster)) {
    cluster <- unique(GSEAresult$cluster)[1]
  }
  if (!(cluster %in% unique(GSEAresult$cluster))) {
    stop("no GSEA result of cluster ", cluster, " was founded")
  }
  if ((by == "GO") & is.null(pathwayIDs)) {
    # GO of level 5，6
    GO2level <- readRDS("Data/GO2level.rds")
    GO_level_5_6 <- GO2level[GO2level$level %in% c(5,6), "GO"] %>% as.character
    GSEAresult %<>% dplyr::filter(ID %in% GO_level_5_6)
  }
  y <- GSEAresult[GSEAresult$cluster == cluster, ]
  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% y$ID) < 1) {
      stop("pathwayIDs provided are not found in GSEA result of cluster ", cluster, ".")
    }
    # y <- GSEAresult %>% dplyr::filter(ID %in% pathwayIDs)
    y <- subset(y, ID %in% pathwayIDs)
  } else {
    y <- dplyr::slice_min(.data = y, order_by = p.adjust, n = showCategory, with_ties = F)
  }
  
  n <- nrow(y)
  geneSets <- geneInCategory(y)
  
  if (n == 0) {
    stop("no enriched term found...")
  } else if (n == 1) {
    g <- graph.empty(0, directed = FALSE)
    g <- add_vertices(g, nv = 1)
    V(g)$name <- y$Description
    V(g)$color <- "red"
    p <- ggraph(g) + geom_node_point(color = "red", size = 5) + geom_node_text(aes_(label = ~name))
    return(p)
  } else {
    id <- y$ID
    geneSets <- geneSets[id]
    w <- matrix(NA, nrow = n, ncol = n)
    colnames(w) <- rownames(w) <- y$Description
    for (i in 1:n) {
      for (j in i:n) {
        w[i, j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
      }
    }
    wd <- reshape2::melt(w)
    wd <- wd[wd[, 1] != wd[, 2], ]
    wd <- wd[!is.na(wd[, 3]), ]
    g <- graph.data.frame(wd[, -3], directed = FALSE)
    E(g)$width = sqrt(wd[, 3] * 5) * cex_line
    g <- delete.edges(g, E(g)[wd[, 3] < 0.2])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    colVar <- y[idx, color]
    V(g)$color <- colVar
  }
  p <- ggraph(g, layout = layout)
  if (length(E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha = 0.8, aes_(width = ~I(width)), 
                            colour = "darkgrey")
  }
  p + geom_node_point(aes_(color = ~color, size = ~size)) + 
    geom_node_text(aes_(label = ~name), repel = TRUE, size = 5 * node_label_cex) + theme_void() + 
    scale_color_continuous(low = "red", high = "blue", 
                           name = color, guide = guide_colorbar(reverse = TRUE)) + 
    scale_size(range = c(3, 8) * node_size_cex)
  
}

################ goplot
getAncestors <- function(ont) {
  Ancestors <- switch(ont,
                      MF = "GOMFANCESTOR",
                      BP = "GOBPANCESTOR",
                      CC = "GOCCANCESTOR",
  )
  return(eval(parse(text=Ancestors)))
}

goplot2 <- function(SeuratObj, cluster = NULL, pathwayIDs = NULL, ont = "BP", showCategory = 10, color = "p.adjust",
                    layout = "sugiyama", geom = "text", label_size = 3) {
  
  if (!("GSEAresult_GO" %in% names(SeuratObj@misc))) {
    stop("No GSEA result of GO were found.")
  }
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[["GSEAresult_GO"]] %>% dplyr::filter(ONTOLOGY == ont)
  if (is.null(cluster)) {
    cluster <- unique(GSEAresult$cluster)[1]
  }
  if (!(cluster %in% unique(GSEAresult$cluster))) {
    stop("no GSEA result of '", ont, "'and cluster ", cluster, " was founded")
  }
  clusterr <- cluster
  # GSEAresult %<>% dplyr::filter(cluster == clusterr)
  GSEAresult <- subset(GSEAresult, cluster == clusterr)

  if (!is.null(pathwayIDs)) {
    if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
      stop("pathwayIDs provided are not found in filtered GSEA result according to the provided parameters.")
    }
    id <- unique(GSEAresult$ID[GSEAresult$ID %in% pathwayIDs])
  } else {
    id <- GSEAresult %>% dplyr::slice_min(order_by = p.adjust, n = showCategory, with_ties = F) %>% dplyr::pull(ID) %>% unique()
  }

  GOSemSim_initial <- getFromNamespace(".initial", "GOSemSim")
  if (!exists(".GOSemSimEnv")) GOSemSim_initial()
  .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
  gotbl <- get("gotbl", envir=.GOSemSimEnv)
  
  GOANCESTOR <- getAncestors(ont)
  anc <- AnnotationDbi::mget(id, GOANCESTOR)
  ca <- Reduce(intersect, anc)
  
  uanc <- unique(unlist(anc))
  uanc <- uanc[!uanc %in% ca]
  dag <- gotbl[gotbl$go_id %in% unique(c(id, uanc)),]
  
  edge <- dag[, c(5, 1, 4)]
  edge <- edge[edge$parent != "all", ]
  node <- unique(gotbl[gotbl$go_id %in% unique(c(edge[,1], edge[,2])), 1:3])
  node$color <- GSEAresult[node$go_id, color]
  
  g <- graph.data.frame(edge, directed=TRUE, vertices=node)
  E(g)$relationship <- edge[,3]
  
  p <- ggraph(g, layout=layout) +
    ## geom_edge_link(aes_(color = ~relationship), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm')) +
    geom_edge_link(aes_(linetype = ~relationship),
                   arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm'),
                   colour="darkgrey") +
    ## geom_node_point(size = 5, aes_(fill=~color), shape=21) +
    geom_node_point(size = 5, aes_(color=~color)) +
    theme_void() +
    scale_color_continuous(low="red", high="blue", name = color,
                           guide=guide_colorbar(reverse=TRUE))
  ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE))
  
  if (geom == "label") {
    p <- p + geom_node_label(aes_(label=~Term, fill=~color), repel=TRUE) +
      scale_fill_continuous(low="red", high="blue", name = color,
                            guide=guide_colorbar(reverse=TRUE), na.value="white")
    ## scale_fill_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE), na.value='white')
  } else {
    p <- p + geom_node_text(aes_(label=~Term), size = label_size, repel=TRUE)
  }
  return(p)
}


# simplifyEnrichmentplot <- function(SeuratObj, by = "GO", pathwayIDs = NULL, showCategory = 10, GO_ont = "BP") {
#   
#   by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb"))
#   GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
# 
#   if (!is.null(pathwayIDs)) {
#     if (sum(pathwayIDs %in% GSEAresult$ID) < 1) {
#       stop("pathwayIDs provided are not found in GSEA result of Seurat object.")
#     }
#     pathwayIDs <- GSEAresult$ID[GSEAresult$ID %in% pathwayIDs]
#   } else {
#     # topath <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% dplyr::slice_head(n = showCategory)
#     pathwayIDs <- GSEAresult %>% dplyr::group_by(cluster) %>% dplyr::slice_min(order_by = p.adjust, n = showCategory, with_ties = F) %>% dplyr::pull(ID)
#   }
# 
#   if (by == "GO") {
#     mat <- GO_similarity(pathwayIDs, ont = GO_ont)
#     simplifyGO(mat, verbose = FALSE)
#   } else if (by == "KEGG") {
#     mat <- term_similarity_from_KEGG(pathwayIDs)
#     simplifyEnrichment(mat, verbose = FALSE)
#   } else if (by == "Reactome") {
#     mat <- term_similarity_from_Reactome(pathwayIDs)
#     simplifyEnrichment(mat, verbose = FALSE)
#   } else if (by == "MSigDb") {
#     mat <- term_similarity_from_MSigDB(pathwayIDs)
#     simplifyEnrichment(mat, verbose = FALSE)
#   }
# }

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

pathwayScatterplot <- function(SeuratObj, by = "GO", pathwayID = NULL, reduction = "umap", 
                               colour = "OrRd", pointsize = 1, label = TRUE, label.size = 4) {
  
  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  gd <- GSEAresult %>% dplyr::mutate(Score=-log10(p.adjust)) %>% reshape2::acast(formula = ID ~ cluster, value.var = 'Score')
  gd[is.na(gd)] <- 0
  gd <- gd[pathwayID, ][as.character(Idents(SeuratObj))]
  # > all.equal(rownames(SeuratObj@meta.data), names(Idents(SeuratObj)))
  # [1] TRUE
  # SeuratObj@meta.data[["pathwayScore"]] <- unname(gd)
  names(gd) <- names(Idents(SeuratObj))
  SeuratObj@meta.data[["pathwayScore"]] <- gd[rownames(SeuratObj@meta.data)]
 
  # p <- Seurat::FeaturePlot(SeuratObj, features = "pathwayScore", pt.size = pointsize, reduction = reduction, label = label) + 
  #   scale_color_gradientn(colours = brewer.pal(n=7,name=colour)) +
  #   ggtitle(sprintf("%s: %s",pathwayID, GSEAresult[GSEAresult$ID == pathwayID, "Description"][1]))
  suppressMessages(p <- Seurat::FeaturePlot(SeuratObj, features = "pathwayScore", pt.size = pointsize, 
                                            reduction = reduction, label = label, label.size = label.size) + 
                     scale_color_gradientn(colours = brewer.pal(n=7,name=colour)) +
                     ggtitle(sprintf("%s: %s",pathwayID, GSEAresult[GSEAresult$ID == pathwayID, "Description"][1])))
  return(p)
}


GOboxplot <- function(SeuratObj, by = "GO", goid = NULL, type = "child", pointsize = 1, flip = FALSE) {
  type <- match.arg(type, choices = c("child", "parent"))
  if (type == "child") {
    nodes <- get_child_nodes(goid)$child_go_id
  } else {
    nodes <- get_parent_nodes(goid)$parent_go_id
  }
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[["GSEAresult_GO"]] %>% dplyr::mutate(Score=-log10(p.adjust))
  godata <- GSEAresult[GSEAresult$ID %in% nodes, ] %>% filter(!is.na(p.adjust)) %>% dplyr::select(cluster, Score)
  p <- ggboxplot(godata, x = "cluster", y = "Score", color = "cluster", size = pointsize,
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






