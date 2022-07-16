
#' Transfer LDA model from reference data to query data
#'
#' @param query_SeuratObj query data
#' @param ref_ldaOUt reference LDA model
#' @param by using which GSEA result of query data to transfer, one of "GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"
#'
#' @return a Seurat object. Transferred LDA model of query data is stored in \code{query_SeuratObj@misc$dist}, which
#' is a list containing clusters-by-topics (θ) and topics-by-pathways (ϕ) matrices.
#' @export
#'
#' @examples
#' \dontrun{
#' query_SeuratObj <- transfer_LDA(query_SeuratObj, ref_ldaOUt, by = "GO")
#' }
#'
transfer_LDA <- function(query_SeuratObj, ref_ldaOUt, by = "GO") {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = query_SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  #create a DocumentTermMatrix
  dd <- GSEAresult %>% dplyr::mutate(score = as.integer(-log2(p.adjust)))
  dtm <- tidytext::cast_dtm(data = dd, document = cluster, term = ID, value = score)
  rowTotals <- apply(dtm, 1, sum)
  colTotals <- apply(dtm, 2, sum)
  dtm <- dtm[rowTotals > 0, colTotals > 0] # remove pathways and clusters that shows no enrichment (word frequency as 0)
  # transfer LDA model from reference to query data
  dist <- modeltools::posterior(ref_ldaOUt, newdata = dtm)
  slot(object = query_SeuratObj, name = 'misc')[["dist"]] <- dist
  return(query_SeuratObj)
}



#' Heatmap showing cosine similarity or pearson correlation between reference data and query data
#'
#' @param ldaOut reference LDA model
#' @param dist transferred model of query data
#' @param a_prefix prefix of label of reference data
#' @param b_prefix prefix of label of query data
#' @param method to calculate cosine similarity or pearson correlation between clusters
#' @param fontsize font size
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' dist <- query_SeuratObj@misc$dist
#' cosine_heatmap2(ldaOut=ref_ldaOUt, dist=dist, a_prefix="reference_", b_prefix="query_", method="cosine")
#' }
#'
#'
cosine_heatmap2 <- function(ldaOut,
                            dist,
                            a_prefix="reference_",
                            b_prefix="query_",
                            method=c("cosine", "pearson"),
                            fontsize=7)
{
  a_mm <- modeltools::posterior(ldaOut)$topics
  rownames(a_mm) <- paste0(a_prefix, rownames(a_mm))
  b_mm <- dist$topics
  rownames(b_mm) <- paste0(b_prefix, rownames(b_mm))
  mm <- rbind(a_mm, b_mm)
  if (method == "cosine") {
    # Calculate cosine similarity
    cc <- philentropy::distance(mm, method = "cosine")
    rownames(cc) <- colnames(cc) <- rownames(mm)
    ccc <- cc[rownames(a_mm), rownames(b_mm)]
    pheatmap::pheatmap(ccc, color = colorRampPalette(rev(c("#000000", "#9932CC", "#EF3B2C","#FFFFCC")))(100),
                       fontsize = fontsize, angle_col = "45")
  } else if (method == "pearson") {
    # Calculate pearson correlation
    pp <- cor(t(mm), method = "pearson")
    ppp <- pp[rownames(a_mm), rownames(b_mm)]
    pheatmap::pheatmap(ppp, fontsize = fontsize, angle_col = "45")
  }
}




#' Hive diagram showing cosine similarity between reference data and query data
#'
#' Hive diagram showing cosine similarity between reference data and query data,
#' calculated with reference LDA model and transferred model.
#'
#' @param ldaOut reference LDA model
#' @param dist transferred model of query data
#' @param a_prefix prefix of label of reference data
#' @param b_prefix prefix of label of query data
#' @param node_size node size
#' @param node_text_size node text size
#'
#' @importFrom dplyr group_by slice_max ungroup mutate
#' @import ggraph
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' dist <- query_SeuratObj@misc$dist
#' cosine_hive(ldaOut=ref_ldaOUt, dist=dist, a_prefix="reference_", b_prefix="query_")
#' }
#'
cosine_hive <- function(ldaOut,
                        dist,
                        a_prefix="reference_",
                        b_prefix="query_",
                        node_size=3,
                        node_text_size=2)
{
  a_gammaDF <- tidytext::tidy(ldaOut, matrix = "gamma")
  a_gammaDF <- a_gammaDF %>% group_by(document) %>% slice_max(order_by = gamma, n=1) %>% ungroup()
  a_gammaDF$document <- paste0(a_prefix, a_gammaDF$document)
  b_gammaDF <- reshape2::melt(dist$topics)
  colnames(b_gammaDF) <- c("document", "topic", "gamma")
  b_gammaDF$document <- paste0(b_prefix, b_gammaDF$document)
  b_gammaDF <- b_gammaDF %>% group_by(document) %>% slice_max(order_by = gamma, n=1) %>% ungroup()
  links <- rbind(a_gammaDF, b_gammaDF) %>% dplyr::mutate(color=topic) %>% .[, c(2,1,3,4)]
  colnames(links) <- c("from", "to", "weight", "color")
  # Calculate cosine similarity
  a_mm <- modeltools::posterior(ldaOut)$topics
  rownames(a_mm) <- paste0(a_prefix, rownames(a_mm))
  b_mm <- dist$topics
  rownames(b_mm) <- paste0(b_prefix, rownames(b_mm))
  mm <- rbind(a_mm, b_mm)
  cc <- philentropy::distance(mm, method = "cosine")
  rownames(cc) <- colnames(cc) <- rownames(mm)
  ccc <- cc[rownames(a_mm), rownames(b_mm)]
  # 只保留cluster之间的余弦相似度最高的2个连线
  cccc <- reshape2::melt(ccc, factorsAsStrings = TRUE) %>% dplyr::mutate(color="cosine", Var1=as.character(Var1), Var2=as.character(Var2))
  a <- cccc %>% group_by(Var1) %>% slice_max(order_by = value, n=1) %>% ungroup()
  b <- cccc %>% group_by(Var2) %>% slice_max(order_by = value, n=1) %>% ungroup()
  c <- unique(rbind(a,b))
  links <- data.table::rbindlist(list(links, c), use.names = F)
  # color
  ee <- unique(c(links$from, links$to))
  colorss <- setNames(scPalette2(length(ee)), ee)
  colorss["cosine"] <- "#696969"
  # colorss <- c(colorss, cosine="#696969") # 相同效果的代码
  # nodedf <- as.data.frame(rbind(cbind(unique(a_gammaDF$topic), "a"), cbind(unique(a_gammaDF$document), "b"), cbind(unique(b_gammaDF$document), "c")))
  # 节点的顺序要改一下,按照2个数据余弦相似度热图中层次聚类的顺序，比较有逻辑（或者按照topic~cluster热图顺序）
  p <- pheatmap::pheatmap(ccc, color = colorRampPalette(rev(c("#000000", "#9932CC", "#EF3B2C","#FFFFCC")))(100), silent = T)
  nodedf <- as.data.frame(rbind(cbind(unique(c(a_gammaDF$topic, b_gammaDF$topic)), "a"), cbind(colnames(ccc)[p$tree_col$order], "b"), cbind(rownames(ccc)[p$tree_row$order], "c")))
  colnames(nodedf) <- c("name", "aa")
  graph <- tidygraph::tbl_graph(nodes = nodedf, edges = links, directed = F, node_key = "name")
  p <- ggraph(graph, 'hive', axis = aa) +
    geom_edge_hive(aes(colour = factor(color), width = weight)) + scale_edge_width(range = c(0.2,1.5)) +
    scale_edge_colour_manual(values = colorss) +
    # geom_axis_hive(aes(colour = aa), size = 2, label = FALSE) +
    ggnewscale::new_scale_colour() +
    geom_node_point(aes(colour=name, shape=aa), size=node_size) + scale_color_manual(values = colorss) +
    geom_node_text(aes(label=name), size=node_text_size, nudge_x = 0.05, nudge_y = 0, repel=F) +
    theme_void() + theme(legend.position = "none") +
    coord_fixed()
  return(p)
}




#' network showing cosine similarity between reference data and query data
#'
#' Network showing cosine similarity between reference data and query data,
#' calculated with reference LDA model and transferred model.
#'
#' @param ldaOut reference LDA model
#' @param dist transferred model of query data
#' @param a_prefix prefix of label of reference data
#' @param b_prefix prefix of label of query data
#' @param layout layout of network
#' @param cos_sim_thresh only links with cosine similarity over threshold will be drawn
#' @param SEED seed
#' @param radius pie size
#' @param width_range edge width
#' @param text_size text size
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph ggraph geom_edge_link scale_edge_width_continuous
#' @importFrom scatterpie geom_scatterpie
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' dist <- query_SeuratObj@misc$dist
#' cosine_network_cluster2(ldaOut=ref_ldaOUt, dist=dist, a_prefix="reference_", b_prefix="query_",
#                         layout="fr", cos_sim_thresh=0.5, radius=0.1, text_size = 3)
#' cosine_network_cluster2(ldaOut=ref_ldaOUt, dist=dist, a_prefix="reference_", b_prefix="query_",
#                         layout="circle", cos_sim_thresh=0.5, radius=0.05, text_size = 3)
#' }
#'
cosine_network_cluster2 <- function(ldaOut,
                                    dist,
                                    a_prefix="reference_",
                                    b_prefix="query_",
                                    layout="circle",
                                    cos_sim_thresh=0.2,
                                    SEED=123,
                                    radius=0.12,
                                    width_range=c(0.1, 0.8),
                                    text_size = 3)
{
  # Calculate cosine similarity
  a_mm <- modeltools::posterior(ldaOut)$topics
  rownames(a_mm) <- paste0(a_prefix, rownames(a_mm))
  b_mm <- dist$topics
  rownames(b_mm) <- paste0(b_prefix, rownames(b_mm))
  mm <- rbind(a_mm, b_mm)
  ttt <- t(combn(rownames(mm), 2))
  # ccc <- plyr::mdply(ttt, .fun = function(x,y){
  #   x=mm[x, ]
  #   y=mm[y, ]
  #   sum(x*y)/(sqrt(sum(x^2))*sqrt(sum(y^2)))
  # })
  cc <- philentropy::distance(mm, method = "cosine")
  rownames(cc) <- colnames(cc) <- rownames(mm)
  ccc <- cbind(as.data.frame(ttt, stringsAsFactors = F), weight=cc[ttt])
  links <- ccc %>% dplyr::filter(weight >= cos_sim_thresh) %>% setNames(c('source', 'target', 'weight')) %>% dplyr::mutate(width=weight)
  nodes <- data.frame(name=unique(c(links$source, links$target)), stringsAsFactors = F)
  clusColor <- setNames(scPalette2(ncol(mm)), colnames(mm))
  links$color <- apply(links, 1, function(x){names(which.max(c(mm[x[1], ], mm[x[2], ])))})
  links$color <- grDevices::adjustcolor(clusColor[links$color], 0.5)
  gg <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  set.seed(SEED)
  p <- ggraph(gg, layout = layout)
  p <- p + geom_edge_link(alpha=0.5, aes_(width=~width, edge_colour=~I(color))) + scale_edge_width_continuous(range = width_range)
  ## Get the matrix data for the pie plot
  ID_Cluster_mat <- cbind(mm[p$data$name, ], p$data[, c("x", "y")])
  ID_Cluster_mat$radius <- radius
  p <- p + ggnewscale::new_scale_fill() + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                                                          cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA, alpha=1, sorted_by_radius=T) +
    scale_fill_manual(name="topic", values=clusColor) +
    coord_fixed() + theme_void()# + theme(legend.position = "none")

  p <- p + ggrepel::geom_text_repel(aes_(x =~ x, y =~ y, label =~ name), size = text_size, show.legend = FALSE,
                                    segment.size = 0.3, force=3, force_pull=2, max.overlaps=15)
  return(p)
}



#' Sankey diagram comparing the reference LDA model and query data
#'
#' After we transfer the reference LDA model to query data (\code{transfer_LDA()}), we can use this function to visualize
#' the comparison of reference LDA model and query data.
#'
#' @param ldaOut reference LDA model
#' @param dist transferred model of query data
#' @param a_prefix prefix of label of reference data
#' @param b_prefix prefix of label of query data
#' @param height output height of sankey diagram
#' @param width output width of sankey diagram
#' @param fontSize font size
#'
#' @importFrom networkD3 sankeyNetwork
#' @importFrom dplyr group_by slice_max ungroup mutate
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' dist <- query_SeuratObj@misc$dist
#' sankey_comparison(ldaOut=ref_ldaOUt, dist=dist, a_prefix="reference_", b_prefix="query_") %>% saveNetwork(file = "sankey_comparison.html")
#' }
#'
#'
sankey_comparison <- function(ldaOut,
                              dist,
                              a_prefix="reference_",
                              b_prefix="query_",
                              height=1600,
                              width = 800,
                              fontSize = 12)
{
  a_gammaDF <- tidytext::tidy(ldaOut, matrix = "gamma")
  a_gammaDF <- a_gammaDF %>% dplyr::mutate(document=paste0(a_prefix, document), topic=paste0("Topic_", topic))
  b_gammaDF <- reshape2::melt(dist$topics) %>% dplyr::mutate(Var1=paste0(b_prefix, Var1), Var2=paste0("Topic_", Var2))
  colnames(b_gammaDF) <- c("document", "topic", "gamma")
  links1 <- a_gammaDF %>% group_by(document) %>% slice_max(order_by = gamma, n=1) %>% ungroup() %>%
    mutate(linkcolor=topic)
  links2 <- b_gammaDF %>% group_by(document) %>% slice_max(order_by = gamma, n=1) %>% ungroup() %>%
    mutate(linkcolor=topic)
  links2 <- links2[,c(2,1,3,4)]
  links <- data.table::rbindlist(list(links1, links2), use.names = F)
  nodes <- data.frame(name=unique(c(links$document, links$topic)), stringsAsFactors = F) %>% mutate(index=0:(n()-1))
  links <- merge(links, nodes, by.x="document", by.y="name")
  links <- merge(links, nodes, by.x="topic", by.y="name")
  links <- links[,c(5,6,3,4)]
  colnames(links) <- c("Source","Target","Value", "linkcolor")
  sankeyNetwork(Links = links, Nodes = nodes, Source = "Source",
                Target = "Target", NodeID = "name", Value = "Value", NodeGroup="name", LinkGroup = "linkcolor",
                fontSize = fontSize, nodeWidth = 30, height = height, width = width, sinksRight=F)
}



#' Oversample the reference data
#'
#' In order to avoid class imbalance while training a classifier, we oversampled the training set.
#' In other words, we randomly selected cells from clusters of a particular cell type to form new clusters,
#' reaching the situation that each cell type owns the same number of clusters
#'
#' @param reference_SeuratObj reference data
#' @param number_clusters goal number of clusters that you need to oversample to
#' @param group_by the column of \code{reference_SeuratObj@meta.data} that indicates cell type
#' @param cluster_by the column of \code{reference_SeuratObj@meta.data} that indicates cluster
#'
#' @return a Seurat object with oversampled expression matrix
#'
oversample <- function(reference_SeuratObj, number_clusters = NULL,
                       group_by = 'cellType', cluster_by = 'seurat_clusters') {

  # make sure that each cell type owns the same number of clusters
  metadata <- reference_SeuratObj@meta.data
  if (is.null(number_clusters)) {
    number_clusters <- max(table(unique(metadata[, c(group_by, cluster_by)])[[group_by]]))
  }
  lll <- lapply(unique(metadata[[group_by]]), function(ct){
    cells <- rownames(metadata)[metadata[[group_by]] == ct]
    cl <- unique(metadata[[cluster_by]][metadata[[group_by]] == ct])
    if (length(cl) < number_clusters) {
      nn <- floor((length(cells)/length(cl))*(9/10)) # cell number of new cluster is 9/10 of average cell number of piticular cell type
      ll <- lapply(seq_len(number_clusters-length(cl)), function(x){
        cc <- sample(cells, size = nn, replace = F)
        data.frame(origCell=cc, newcell=paste0(cc, "pseudo", x), cluster=paste0("pseudo_", x, "_", ct),
                   cellType=ct, stringsAsFactors = F)
      })
      do.call("rbind", ll)
    }
  })
  celltypeNew <- do.call("rbind", lll)
  # new expression matrix
  exprMatrix <- reference_SeuratObj@assays$RNA@counts
  if (!inherits(x = exprMatrix, what = 'dgCMatrix')) {
    exprMatrix <- as(exprMatrix, "dgCMatrix")
  }
  newexp <- exprMatrix[, celltypeNew$origCell]
  colnames(newexp) <- celltypeNew$newcell
  exprMatrix <- cbind(exprMatrix, newexp)
  # new meta data
  metadata$origCell <- rownames(metadata)
  newmetadata <- metadata[celltypeNew$origCell, ]
  rownames(newmetadata) <- celltypeNew$newcell
  newmetadata[[group_by]] <- celltypeNew$cellType
  newmetadata[[cluster_by]] <- celltypeNew$cluster
  metadata <- rbind(metadata, newmetadata)
  # oversampled seurat object
  reference_SeuratObj <- Seurat::CreateSeuratObject(counts = exprMatrix, meta.data = metadata)
  return(reference_SeuratObj)
}



#' Oversample the reference data
#'
#' Oversample the reference data and perform GSEA, topic modelling.
#' In order to avoid class imbalance while training a classifier, we oversampled the training set.
#' In other words, we randomly selected cells from clusters of a particular cell type to form new clusters,
#' reaching the situation that each cell type owns the same number of clusters.
#'
#' @param reference_SeuratObj reference data
#' @param number_clusters goal number of clusters that you need to oversample to
#' @param group_by the column of \code{reference_SeuratObj@meta.data} that indicates cell type
#' @param cluster_by the column of \code{reference_SeuratObj@meta.data} that indicates cluster
#' @param species species of the reference data
#' @param by database used to perform GSEA. GO KEGG Reactome MSigDb WikiPathways DO NCG DGN.
#' @param k number of topics.
#' @param method method used for fitting a LDA model; currently "VEM" or "Gibbs" are supported.
#'
#' @return a Seurat object with oversampled expression matrix and topic-model result.
#' @export
#'
#' @examples
#' \dontrun{
#' reference_SeuratObj <- oversample_ref(reference_SeuratObj, number_clusters = 10,
#' group_by = 'cellType', cluster_by = 'seurat_clusters',
#' species = "Homo sapiens", by = 'GO', k = NULL, method = "VEM")
#' }
#'
#'
oversample_ref <- function(reference_SeuratObj, number_clusters = NULL,
                           group_by = 'cellType', cluster_by = 'seurat_clusters',
                           species = "Homo sapiens", by = 'GO', k = NULL, method = "VEM") {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  reference_SeuratObj <- oversample(reference_SeuratObj, number_clusters = number_clusters,
                                    group_by = group_by, cluster_by = cluster_by)
  reference_SeuratObj <- readData(data = reference_SeuratObj, type = 'Seurat', species = species)
  reference_SeuratObj <- CalMTpercent(reference_SeuratObj, by = "use_internal_data")
  reference_SeuratObj <- ScFunQC(reference_SeuratObj, plot = F)
  # We don't need to re-clustering the reference_SeuratObj
  Seurat::Idents(reference_SeuratObj) <- reference_SeuratObj@meta.data[[cluster_by]]
  reference_SeuratObj <- Seurat::NormalizeData(reference_SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  reference_SeuratObj <- Seurat::ScaleData(reference_SeuratObj, features = rownames(reference_SeuratObj))
  message("Calculating differentially expressed genes for reference data......")
  SeuratObj.markers <- Seurat::FindAllMarkers(reference_SeuratObj, only.pos = TRUE, min.pct = 0.0001, logfc.threshold = 0.0001, return.thresh=0.9)
  slot(object = reference_SeuratObj, name = 'misc')[["Allmarkers"]] <- SeuratObj.markers
  message("Performing GSEA on reference data......")
  reference_SeuratObj <- RunGSEA(reference_SeuratObj, by = by)
  message("Performing topic modelling on reference data......")
  if (is.null(k)) {
    k <- length(unique(reference_SeuratObj@meta.data[[group_by]]))
  }
  reference_SeuratObj <- runLDA(reference_SeuratObj, by = by, k = k, method = method, SEED = 1234, plot = F)
  return(reference_SeuratObj)
}




#' Predict cell types of query data based on reference data
#'
#' Predict cell type of query data based on reference data. We use caret R package to train a svm classifier.
#'
#' @param query_SeuratObj  query data
#' @param reference_SeuratObj  reference data
#' @param species species of the reference data
#' @param by  database used to perform GSEA on reference data. This parameter works when \code{reference_SeuratObj} does not contain topic model.
#' @param k  number of topics. This parameter works when \code{reference_SeuratObj} does not contain topic model.
#' @param LDAmethod  method used for fitting a LDA model on reference data; currently "VEM" or "Gibbs" are supported. This parameter works when \code{reference_SeuratObj} does not contain topic model.
#'
#' @import caret
#'
#' @return data frame with 2 columns including clusters of query data, and the predicted cell type based on reference data.
#' @export
#'
#' @examples
#' \dontrun{
#' predictFun(query_SeuratObj, reference_SeuratObj, species = "Homo sapiens", by = 'GO', k = NULL, LDAmethod = "VEM")
#' }
#'
#'
predictFun <- function(query_SeuratObj, reference_SeuratObj, species = "Homo sapiens",
                       by = 'GO', k = NULL, LDAmethod = "VEM") {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  if (!"ldaOut" %in% names(slot(object = reference_SeuratObj, name = 'misc'))) {
    if (!paste0("GSEAresult_", by) %in% names(slot(object = reference_SeuratObj, name = 'misc'))) {
      if (!"species" %in% names(slot(object = reference_SeuratObj, name = 'misc'))) {
        reference_SeuratObj <- readData(data = reference_SeuratObj, type = 'Seurat', species = species)
      }
      if (!"featureData" %in% names(slot(object = reference_SeuratObj, name = 'misc'))) {
        reference_SeuratObj <- DetectGeneIDtype(reference_SeuratObj)
      }
      if (!"Allmarkers" %in% names(slot(object = reference_SeuratObj, name = 'misc'))) {
        # We don't need to re-clustering the reference_SeuratObj
        message("Calculating differentially expressed genes for reference data......")
        SeuratObj.markers <- Seurat::FindAllMarkers(reference_SeuratObj, only.pos = TRUE, min.pct = 0.0001, logfc.threshold = 0.0001, return.thresh=0.9)
        slot(object = reference_SeuratObj, name = 'misc')[["Allmarkers"]] <- SeuratObj.markers
      }
      message("Performing GSEA on reference data......")
      reference_SeuratObj <- RunGSEA(reference_SeuratObj, by = by)
    }
    message("Performing topic modelling on reference data......")
    if (is.null(k)) {
      k <- length(unique(Seurat::Idents(reference_SeuratObj)))
    }
    reference_SeuratObj <- runLDA(reference_SeuratObj, by = by, k = k, method = LDAmethod, SEED = 1234, plot = F)
  }
  ref_ldaOUt <- reference_SeuratObj@misc$ldaOut
  query_SeuratObj <- transfer_LDA(query_SeuratObj = query_SeuratObj, ref_ldaOUt, by = "GO")
  dist <- query_SeuratObj@misc$dist
  # training set and testing set
  message("training a svm classifier......")
  training <- modeltools::posterior(ref_ldaOUt)$topics
  testing <- dist$topics
  colnames(training) <- colnames(testing) <- paste0("Topic_", colnames(training))
  training <- as.data.frame(training, stringsAsFactors = F)
  testing <- as.data.frame(testing, stringsAsFactors = F)
  training$Class<- factor(rownames(training))
  # train a svm classifier
  ctrl <- trainControl(method = "repeatedcv", number=10, repeats = 3, search = "random")
  set.seed(123)
  svm_Radial <- train(
    Class ~ .,
    data = training,
    method = "svmRadial",
    preProc = c("center", "scale"),
    tuneLength = 30,
    trControl = ctrl
  )
  message("predicting cell types of the query data......")
  prediction <- predict(svm_Radial, newdata = testing)
  df <- data.frame(query = rownames(testing), prediction = as.character(prediction), stringsAsFactors = F)
  return(df)
}









