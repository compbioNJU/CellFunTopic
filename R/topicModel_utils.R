#' calculate metrics to find the best number of topics
#'
#' @param dtm Object of class "DocumentTermMatrix"
#' @param max.k max number of k
#' @param method The method to be used for fitting; currently method = "VEM" or method= "Gibbs" are supported.
#' @param SEED seed
#'
#' @importFrom modeltools posterior
#'
#' @return list of 4 metrics: perplexity, loglikelihood, alpha, entropy, cosine similarity
#' @export
#'
#' @examples
#' \dontrun{
#' ldas <- lda.rep(dtm=dtm, max.k=10, method='VEM', SEED=123)
#' plotindex(ldas)
#' }
#'
#'
lda.rep<-function(dtm, max.k=5, method='VEM', SEED=1234){
  cat('Calculating Metrics to choose the best number of topics.\n')
  # model<-vector('list')
  perp <- loglik <- alpha <- entropy <- cosine <- c()
  for(i in 2:max.k){
    cat('---------------- k=',i,' --------------\n')
    # model[[i-1]] <- LDA(dtm, k = i, method = method, control = list(seed = SEED, verbose=0))
    modell <- topicmodels::LDA(dtm, k = i, method = method, control = list(seed = SEED, verbose=0))
    perp[i-1] <- topicmodels::perplexity(modell)  # perplexity
    loglik[i-1] <- logLik(modell)[1]  # loglikelihood
    alpha[i-1] <- slot(modell, "alpha") # alpha
    entropy[i-1] <- mean(apply(posterior(modell)$topics, 1, function(z) - sum(z * log(z))))  # entropy
    ccc <- t(combn(1:i, 2))
    cosine[i-1] <- purrr::map2(ccc[,1], ccc[,2], .f = function(x,y){
      mattt <- posterior(modell)$terms
      x=mattt[x,]
      y=mattt[y,]
      sum(x*y)/sqrt(sum(x^2)*sum(y^2))
    }) %>% unlist %>% mean
  }
  return(list(#model=model,
    perp=perp, loglik=loglik, alpha=alpha, entropy=entropy, cosine=cosine))
}


#' plot metrics to find the best number of topics
#'
#' @param ldas list of 4 metrics(perplexity, loglikelihood, alpha, entropy, cosine similarity), result of \code{lda.rep()}.
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' ldas <- lda.rep(dtm=dtm, max.k=10, method='VEM', SEED=123)
#' plotindex(ldas)
#' }
#'
#'
plotindex <- function(ldas) {
  n <- length(ldas$perp)+1
  est=data.frame(x=rep(2:n, 3), y=c(ldas$perp, ldas$loglik, ldas$cosine),
                 type=rep(c('perplexity', 'loglikelihood', 'cosine similarity'), each=n-1))
  ggplot(est, aes(x, y, color=type)) +
    geom_line() + geom_point()+ scale_x_continuous(breaks = 0:n) +
    facet_wrap(~ type, scales = "free", ncol=3)+
    labs(x = "number of topic")
}


#' Run topic modeling with GSEA result
#'
#' @param SeuratObj Object of class "Seurat"
#' @param by using which GSEA result to run topic modeling, one of "GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"
#' @param k number of topics. Can be a specific k, if NULL, calculate best k automatically(time-consuming)
#' @param method method used for fitting a LDA model; currently "VEM" or "Gibbs" are supported.
#' @param SEED seed
#' @param plot whether or not to plot metrics used to decide the number of topics
#'
#' @return a Seurat object with topic modeling result stored in the \code{SeuratObj@misc$ldaOut}
#' @export
#'
#' @examples
#' \dontrun{
#' k <- 10
#' SeuratObj <- runLDA(SeuratObj, by = "GO", k = k, method = "VEM", SEED = 1234, plot = T)
#' SeuratObj@misc$ldaOut # where the topic modeling result is stored
#' }
#'
#'
runLDA <- function(SeuratObj, by = "GO", k=NULL, method='VEM', SEED=1234, plot=TRUE) {
  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  #create a DocumentTermMatrix
  dd <- GSEAresult %>% dplyr::mutate(score = as.integer(-log2(p.adjust)))
  dtm <- tidytext::cast_dtm(data = dd, document = cluster, term = ID, value = score)
  rowTotals <- apply(dtm, 1, sum)
  colTotals <- apply(dtm, 2, sum)
  dtm <- dtm[rowTotals > 0, colTotals > 0] # remove pathways and clusters that shows no enrichment (word frequency as 0)
  # dtm <- dtm[slam::row_sums(dtm) > 0,] # code of same effect
  if (is.null(k)) {
    max.k <- length(unique(GSEAresult$cluster))
    ldas <- lda.rep(dtm=dtm, max.k=max.k, method=method, SEED=SEED)
    slot(object = SeuratObj, name = 'misc')[["ldas"]] <- ldas
    bestK1 <- seq(2,max.k)[which.min(ldas$perp)]
    bestK2 <- seq(2,max.k)[which.max(ldas$loglik)]
    k <- min(bestK1, bestK2)
    cat('The number of topics is set to', k, '.\n')
    if (plot) {
      plotindex(ldas)
    }
  }
  ldaOut <- topicmodels::LDA(dtm, k = k, method = method, control = list(seed = SEED))
  slot(object = SeuratObj, name = 'misc')[["ldaOut"]] <- ldaOut
  return(SeuratObj)
}


#' barplot to show top terms of every topic
#'
#' @param betaDF topic~term probability data.frame
#' @param topics topics to show
#' @param topn top n terms of topic
#' @param axis.text.y.size size of y axis text
#'
#' @importFrom dplyr filter group_by slice_max ungroup
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' betaDF <- tidytext::tidy(ldaOut, matrix = "beta")
#' pws <- SeuratObj@misc$GSEAresult_GO %>% dplyr::select(ID, Description) %>% unique %>% tibble::deframe()
#' betaDF %<>% dplyr::mutate(descrip=unname(pws[term]))
#' hist_topic_term(betaDF, topics=1, topn=20, axis.text.y.size=5)
#' hist_topic_term(betaDF, topics=c(1,2,3,4), topn=20, axis.text.y.size=5)
#' }
#'
#'
hist_topic_term <- function(betaDF, topics, topn=5, axis.text.y.size=5) {
  if (is.null(topics)) topics <- head(unique(betaDF$topic), 4)
  top_terms <- betaDF %>% dplyr::filter(topic %in% topics) %>%
    group_by(topic) %>% slice_max(order_by = beta, n=topn, with_ties=F) %>% ungroup()
  pp <- top_terms %>%
    mutate(descrip = reorder(descrip, beta)) %>%
    ggplot(aes(descrip,beta,fill=factor(topic))) +
    geom_bar(stat="identity") +
    coord_flip(expand = TRUE) +
    labs(title = paste0("Top ", topn, " terms of Topics ", paste(topics, collapse=",")), fill="topic", x="terms", y="probability") +
    facet_wrap(~ topic, scales = "free", ncol=2) +
    theme(axis.text.y = element_text(size=axis.text.y.size),
          # axis.text.x = element_text(size=5),
          panel.background = element_blank(),
          panel.border=element_rect(fill="transparent",color="black"),
          strip.background = element_rect(fill = "light gray", color = "black"),
          plot.title = element_text(lineheight = 610,colour = "black",size = 15))
  return(pp)
}

#' barplot to show probability of assigned topics in clusters
#'
#' @param gammaDF cluster~topic probability data.frame
#' @param clusters clusters to show
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' gammaDF <- tidytext::tidy(ldaOut, matrix = "gamma")
#' hist_cluster_topic(gammaDF, clusters=head(Seurat::Idents(SeuratObj)))
#' }
#'
hist_cluster_topic <- function(gammaDF, clusters) {
  if (is.null(clusters)) clusters <- head(unique(gammaDF$document), 3)
  df <- gammaDF %>% dplyr::filter(document %in% clusters)
  pp <- df %>% ggplot(aes(factor(topic), gamma, fill=factor(topic))) +
    geom_bar(stat="identity") +
    coord_flip() +
    facet_wrap(~ document, ncol=3) +
    labs(fill="topic", x="topic") +
    theme(axis.text.y = element_text(size=11),
          axis.text.x = element_text(size=10),
          panel.background = element_blank(),
          panel.border=element_rect(fill="transparent",color="black"),
          strip.background = element_rect(fill = "light gray", color = "black"),
          plot.title = element_text(lineheight = 610,colour = "black",size = 15))
  return(pp)
}


#' network showing top terms of topics
#'
#' @param betaDF topic~term probability data.frame
#' @param topn number of top terms of each topic
#'
#' @importFrom igraph graph_from_data_frame
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' betaDF <- tidytext::tidy(ldaOut, matrix = "beta")
#' topicNW1(betaDF, topn=10)
#' }
#'
topicNW1 <- function(betaDF, topn) {
  df <- betaDF %>% dplyr::group_by(topic) %>% dplyr::slice_max(order_by = beta, n=topn, with_ties=F)
  ll <- split(df$term, df$topic) %>% lapply(embed, dimension = 2)
  edgedf <- do.call(rbind, ll)[, 2:1] %>% as.data.frame(stringsAsFactors = F) %>%
    cbind(rep(unique(df$topic), each=topn-1)) %>% setNames(c("source", "target", "topic"))
  vv <- unique(c(edgedf$source, edgedf$target))
  gg <- graph_from_data_frame(edgedf[,1:2], directed = F, vertices = vv)
  plot(gg, layout=layout.fruchterman.reingold,
       vertex.size=1.5,
       vertex.color='gray90',
       vertex.label.cex=0.5,
       edge.color=as.integer(edgedf$topic),
       edge.width = 2,
       edge.curved=TRUE,
       main=paste0('Top ', topn, ' terms of ', length(unique(edgedf$topic)), ' Topics'))
}

#' network showing top terms of topics
#'
#' @param betaDF topic~term probability data.frame
#' @param topn number of top terms of each topic
#'
#' @importFrom igraph graph_from_data_frame
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' betaDF <- tidytext::tidy(ldaOut, matrix = "beta")
#' topicNW2(betaDF, topn=10)
#' }
#'
topicNW2 <- function(betaDF, topn) {
  df <- betaDF %>% dplyr::group_by(topic) %>% dplyr::slice_max(order_by = beta, n=topn, with_ties=F) %>% dplyr::ungroup()
  edgedf <- df %>% dplyr::mutate(topic=paste0("Topic ", topic)) %>% dplyr::select(topic, term, beta) %>%
    setNames(c("source", "target", "weight"))# %>%
  # dplyr::mutate(width=2*(weight-min(weight))/(max(weight)-min(weight))+1)   # 这样好像不对，会出现NaN
  ww <- edgedf$weight
  edgedf$width <- 2*(ww-min(ww))/(max(ww)-min(ww))+1
  colorss <- setNames(scPalette2(length(unique(edgedf$source))), unique(edgedf$source))
  edgedf$color <- colorss[edgedf$source]
  nodes <- c(unique(edgedf$source), unique(edgedf$target)) %>%
    cbind(c(unname(colorss), rep("#D3D3D3", length(unique(edgedf$target))))) %>%
    as.data.frame() %>%  setNames(c("node","color")) %>%
    dplyr::mutate(size=c(rep(8, length(unique(edgedf$source))), rep(5, length(unique(edgedf$target)))),
                  label.cex=c(rep(1, length(unique(edgedf$source))), rep(0.5, length(unique(edgedf$target)))))
  gg <- graph_from_data_frame(edgedf, directed = T, vertices = nodes)
  # ww <- E(gg)$weight
  # E(gg)$width <- 2*(ww-min(ww))/(max(ww)-min(ww))+1
  V(gg)$frame.color <- NA
  plot(gg, layout=layout_nicely,#layout_with_kk,
       edge.arrow.size=0.3,
       main=paste0('Top ', topn, ' terms of ', length(unique(edgedf$source)), ' Topics'))
}



#' echart4r network showing top terms of topics
#'
#' echart4r network showing top terms of topics with tooltips showing description of pathway ID
#'
#' @param betaDF topic~term probability data.frame
#' @param topn top n terms of topic
#' @param pws named vector, pathway description with ID as names
#'
#' @importFrom echarts4r e_charts e_graph e_graph_nodes e_graph_edges e_tooltip e_title e_legend e_toolbox e_toolbox_feature
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' betaDF <- tidytext::tidy(ldaOut, matrix = "beta")
#' pws <- SeuratObj@misc$GSEAresult_GO %>% dplyr::select(ID, Description) %>% unique %>% tibble::deframe()
#' topicNW3(betaDF, topn=10, pws)
#' }
#'
topicNW3 <- function(betaDF, topn, pws) {
  df <- betaDF %>% dplyr::group_by(topic) %>% dplyr::slice_max(order_by = beta, n=topn, with_ties=F) %>% dplyr::ungroup()
  edgedf <- df %>% dplyr::mutate(topic=paste0("Topic ", topic)) %>% dplyr::select(topic, term, beta)) %>%
    setNames(c("source", "target", "weight")) %>%
    dplyr::mutate(target=paste(target, unname(pws[target]), sep=": "))  # 适用于多种情况 GO/KEGG
                  # target=paste(target, GOfuncR::get_names(target)$go_name, sep=": ")) # 只适用于GO的情况
  ww <- edgedf$weight
  edgedf$width <- 2*(ww-min(ww))/(max(ww)-min(ww))+1
  aa <- unique(edgedf$source)
  bb <- unique(edgedf$target)
  nodes <- data.frame(node=c(aa,bb),
                      size=c(rep(20, length(aa)), rep(5, length(bb)))) %>%
    dplyr::mutate(value=node,
                  label.cex=c(rep(1, length(aa)), rep(0.5, length(bb))))
  ddd <- edgedf %>% dplyr::group_by(target) %>% slice_max(order_by = weight, n=1, with_ties=F) %>% dplyr::select(target, source) %>%
    tibble::deframe()
  nodes$cat <- c(aa, ddd[bb])
  e_charts() %>%
    e_graph(layout = "force",
            roam = T,
            draggable = T,
            symbolKeepAspect = T,
            focusNodeAdjacency = T,
            force = list(repulsion = 100),
            edgeSymbol = list('none', 'arrow'),
            edgeSymbolSize = 5,
            lineStyle = list(color = 'source'),
            label = list(show = T, color = "#000000", fontWeight = "normal", fontSize = 5)
    ) %>%
    e_graph_nodes(nodes=nodes, names=node, value=value, size=size, category=cat) %>%
    e_graph_edges(edges=edgedf, source, target) %>%
    e_labels(
      formatter = htmlwidgets::JS(
        paste0(
          '
          function(params) {
          if (params.value) {
          const text = params.value.split(": ")
          const status = (text[0].substring(0,5) == "Topic")? "ff" : "ee"
          return(`{${status}|${text[0]}}`)
          }
          }
          '
        )
      ),
      rich = list(
        ff = list(borderWidth = 0, padding = 0, fontSize = 13),
        ee = list(borderWidth = 0, padding = 0, fontSize = 8)
      )
    ) %>%
    e_tooltip(
      formatter = htmlwidgets::JS(
        paste0(
          '
          function(params) {
          if (params.value) {
          const text = params.value.split(": ")
          const label = (text[0].substring(0,5) == "Topic")? text[0] : `${text[0]}<br>${text[1]}`
          return(label)
          }
          }
          '
        )
      )
    ) %>%
    e_title(
      text = paste0('Top ', topn, ' terms of ', length(aa), ' Topics')
    ) %>%
    e_legend(right=0, top=20) %>%
    e_toolbox() %>%
    e_toolbox_feature(feature = c("saveAsImage", "dataView"))
}

# topicNW3 <- function(betaDF, topn) {
#   df <- betaDF %>% dplyr::group_by(topic) %>% dplyr::slice_max(order_by = beta, n=topn, with_ties=F) %>% dplyr::ungroup()
#   edgedf <- df %>% dplyr::mutate(topic=paste0("Topic ", topic)) %>% dplyr::select(topic, term, beta) %>%
#     setNames(c("source", "target", "weight"))
#   ww <- edgedf$weight
#   edgedf$width <- 2*(ww-min(ww))/(max(ww)-min(ww))+1
#   aa <- unique(edgedf$source)
#   bb <- unique(edgedf$target)
#   nodes <- data.frame(node=c(aa,bb),
#                       size=c(rep(20, length(aa)), rep(5, length(bb)))) %>%
#     dplyr::mutate(value=node,
#                   label.cex=c(rep(1, length(aa)), rep(0.5, length(bb))))
#   ddd <- edgedf %>% dplyr::group_by(target) %>% slice_max(order_by = weight, n=1, with_ties=F) %>% dplyr::select(target, source) %>%
#     tibble::deframe()
#   nodes$cat <- c(aa, ddd[bb])
#   e_charts() %>%
#     e_graph(layout = "force",
#             roam = T,
#             draggable = T,
#             symbolKeepAspect = T,
#             focusNodeAdjacency = T,
#             force = list(repulsion = 100),
#             edgeSymbol = list('none', 'arrow'),
#             edgeSymbolSize = 5,
#             lineStyle = list(color = 'source'),
#             label = list(show = T, color = "#000000", fontWeight = "normal", fontSize = 8)
#     ) %>%
#     e_graph_nodes(nodes=nodes, names=node, value=value, size=size, category=cat) %>%
#     e_graph_edges(edges=edgedf, source, target) %>%
#     e_tooltip() %>%
#     e_title(
#       text = paste0('Top ', topn, ' terms of ', length(aa), ' Topics')
#     ) %>%
#     e_legend(right=0, top=10) %>%
#     e_toolbox() %>%
#     e_toolbox_feature(feature = c("saveAsImage", "dataView"))
# }


#' word cloud of top terms of topic
#'
#' @param betaDF topic~term probability data.frame
#' @param pws named vector, pathway description with ID as names
#' @param topic topic to draw
#' @param topn top n terms of topic
#'
#' @importFrom ggwordcloud geom_text_wordcloud_area
#'
#' @rdname wordcloud_topic
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' betaDF <- tidytext::tidy(ldaOut, matrix = "beta")
#' pws <- SeuratObj@misc$GSEAresult_GO %>% dplyr::select(ID, Description) %>% unique %>% tibble::deframe()
#' wordcloud_topic(betaDF, pws, topic=1, topn=20)
#' wordcloud_topic_3D(betaDF, pws, topic=1, topn=20)
#' }
#'
wordcloud_topic <- function(betaDF, pws, topic, topn=20) {
  pp <- betaDF %>% filter(.data$topic == .env$topic) %>% slice_max(order_by = beta, n=topn, with_ties=F) %>% mutate(term=pws[term]) %>%
    ggplot(aes(label=term, size=beta, color=beta)) +
    geom_text_wordcloud_area(show.legend = T) +
    scale_color_gradientn(colours = rev(scPalette2(20))) + # scale_size(range = c(1,6)) +
    theme_minimal() + #ggtitle(paste0("Topic ", topic)) +
    theme(legend.position = "none")# +
    # theme(plot.margin = margin(t = 0,
    #                            r = 0,
    #                            b = 0,
    #                            l = 0,
    #                            unit = "cm"))
  return(pp)
}

#' @rdname wordcloud_topic
#' @export
#'
#'
wordcloud_topic_3D <- function(betaDF, pws, topic, topn=20) {
  pal2 <- colorRampPalette(brewer.pal(9,"Set1"))(20)
  dd <- betaDF %>% filter(.data$topic == .env$topic) %>% dplyr::slice_max(order_by = beta, n=topn, with_ties=F)
  # suppressWarnings(wordcloud(pws[dd$term],dd$beta, scale=c(2,.1), min.freq=0,
  #           max.words=20, random.order=F, rot.per=.15, colors=pal2, ordered.colors=T))
  # # JS 3D wordcloud
  wordcloud2::wordcloud2(dd %>% mutate(descrip=unname(pws[term])) %>% dplyr::select(descrip, beta), size = 0.1,color= "random-dark")
}


#' Sankey diagram showing the best assigned topic of each cluster
#'
#' Create a D3 JavaScript Sankey diagram which shows the best assigned topic of each cluster.
#'
#' @param gammaDF cluster~topic probability data.frame
#' @param topn top n terms of each topic
#' @param plotHeight height of plot
#'
#' @importFrom networkD3 sankeyNetwork
#' @importFrom dplyr group_by slice_max ungroup mutate n
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' gammaDF <- tidytext::tidy(ldaOut, matrix = "gamma")
#' plot_sankey(gammaDF, topn=1)
#' plot_sankey(gammaDF, topn=2) %>% networkD3::saveNetwork(file = "plot_sankey.html")
#' }
#'
plot_sankey <- function(gammaDF, topn=1, plotHeight=600) {
  links <- gammaDF %>% group_by(document) %>% slice_max(order_by = gamma, n=topn) %>% ungroup() %>%
    mutate(document=as.character(document))
  nodes <- data.frame(name=c(unique(links$document), paste0("Topic_", sort(unique(links$topic))))) %>% mutate(index=0:(n()-1))
  links <- links %>% mutate(topic=paste0("Topic_", topic)) %>% as.data.frame()
  links <- merge(links, nodes, by.x="document", by.y="name")
  links <- merge(links, nodes, by.x="topic", by.y="name")
  links <- links[,c(4,5,3,1)]
  colnames(links) <- c("Source","Target","Value", "topic")
  sankeyNetwork(Links = links, Nodes = nodes, Source = "Source",
                Target = "Target", NodeID = "name", Value = "Value", NodeGroup="name", LinkGroup = "topic",
                fontSize = 12, nodeWidth = 30, height = plotHeight, width = 800, sinksRight=F)
}


#' UMAP on cluster-topic probability matrix
#'
#' Each point represents clusters. Point color indicates the topic with the highest probability of each cluster.
#'
#' @param ldaOut topic modeling result
#'
#' @importFrom modeltools posterior
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' umap_cluster(ldaOut)
#' }
#'
#'
umap_cluster <- function(ldaOut) {
  mat <- posterior(ldaOut)$topics  # cluster-topic probability matrix
  em <- uwot::umap(mat, n_neighbors = 3)
  colnames(x = em) <- paste0("UMAP_", 1:ncol(x = em))
  rownames(x = em) <- rownames(mat)
  emdf <- em %>% as.data.frame() %>% tibble::rownames_to_column("cluster") %>%
    mutate(topic=apply(mat, 1, function(x){names(which.max(x))})) # 每个cluster概率最高的topic
  emdf$cluster <- factor(emdf$cluster, levels = unique(emdf$cluster))
  emdf$topic <- factor(emdf$topic, levels = sort(as.numeric(unique(emdf$topic))))
  # umap散点图按topic着色
  pp <- ggplot(emdf, aes(UMAP_1, UMAP_2, color=topic)) +
    geom_point(size=1) + theme_classic() +
    scale_color_manual(name="topic", values=setNames(scPalette2(ncol(mat)), colnames(mat))) +
    guides(color=guide_legend(title = "topic", override.aes = list(size=3))) +
    # geom_text(inherit.aes = T, aes(label=cluster), size=3) + #coord_fixed() +
    ggrepel::geom_text_repel(aes(label=cluster), max.overlaps = 20, segment.colour="black") +
    ggtitle("UMAP on cluster-topic probability matrix")
  return(pp)
}


#' Heatmap showing probability between topics and clusters
#'
#' @param ldaOut topic modeling result
#'
#' @importFrom modeltools posterior
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' cluster_topic_hmp(ldaOut)
#' cluster_topic_hmp(ldaOut,  cluster_cols = F)
#' }
#'
cluster_topic_hmp <- function(ldaOut, cluster_rows = T, cluster_cols = T) {
  mm <- posterior(ldaOut)$topics
  colnames(mm) <- paste0('Topic ', colnames(mm))
  pheatmap::pheatmap(mm, cluster_rows = cluster_rows, cluster_cols = cluster_cols, angle_col = "0", main = "probability between topic and cluster",
                     color = colorRampPalette(c('white', brewer.pal(n=9,name="OrRd")))(100))
}


#' Heatmap showing cosine similarity between topics
#'
#' calculated by topic~term probability matrix
#'
#' @param ldaOut topic modeling result
#'
#' @importFrom modeltools posterior
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' cosineheatmap(ldaOut)
#' }
#'
cosineheatmap <- function(ldaOut) {
  mm <- posterior(ldaOut)$terms  # topic~term probability matrix
  cc <- suppressMessages(philentropy::distance(mm, method = "cosine"))
  rownames(cc) <- colnames(cc) <- rownames(mm)
  pheatmap::pheatmap(cc, color = colorRampPalette(rev(c("#000000", "#9932CC", "#EF3B2C","#FFFFCC")))(100), angle_col = "0",
                     fontsize = 9, main = "cosine similarity between topics")#, cellwidth = 30, cellheight = 30)
}


#' Network showing cosine similarity between clusters
#'
#' Network of cosine similarity between clusters, calculated by cluster~topic probability matrix. Width of edge shows cosine similarity between clusters, node pie shows topic probability distribution of each cluster.
#'
#' @param ldaOut topic modeling result
#' @param layout network layout
#' @param cos_sim_thresh only shows edge of cosine similarity greater than specific threshold
#' @param SEED seed
#' @param radius node pie radius
#' @param width_range range of width of edges
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph ggraph geom_edge_link scale_edge_width_continuous
#' @importFrom scatterpie geom_scatterpie
#' @importFrom modeltools posterior
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' ldaOut <- SeuratObj@misc$ldaOut
#' cosine_network_cluster(ldaOut, layout="circle", cos_sim_thresh=0.2, SEED=123, radius=0.12, width_range=c(0.1, 0.8))
#' cosine_network_cluster(ldaOut, layout="fr", cos_sim_thresh=0.2, SEED=123, radius=0.12, width_range=c(0.1, 0.8))
#' }
#'
cosine_network_cluster <- function(ldaOut,
                                   layout="circle",
                                   cos_sim_thresh=0.2,
                                   SEED=123,
                                   radius=0.12,
                                   width_range=c(0.1, 0.8),
                                   text_size = 3)
{
  mm <- posterior(ldaOut)$topics  # cluster~topic probability matrix
  ccc <- t(combn(rownames(mm), 2))
  ccc <- plyr::mdply(ccc, .fun = function(x,y){
    x=mm[x, ]
    y=mm[y, ]
    sum(x*y)/(sqrt(sum(x^2))*sqrt(sum(y^2)))
  })
  links <- ccc %>% dplyr::filter(V1 >= cos_sim_thresh) %>% setNames(c('source', 'target', 'weight')) %>% dplyr::mutate(width=weight)
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


#' Network of cosine similarity between terms
#'
#' Links showing cosine similarity between terms, calculated by topic~term probability matrix or GSEA result.
#' Node pie shows topic distribution of each term or enrichment score distribution of each term.
#'
#' @param SeuratObj Object of class "Seurat"
#' @param cosine_cal_by calculate cosine similarity by topic~term probability matrix or GSEA result. "Topic modeling", "GSEA result".
#' @param pie_by what node pie shows. topic distribution or enrichment score distribution of each term. "Topic modeling", "GSEA result".
#' @param GSEA_by using which GSEA result to run topic modeling, one of "GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"
#' @param topn top n terms of each topic or each cluster
#' @param layout network layout
#' @param cos_sim_thresh only shows edge of cosine similarity greater than specific threshold
#' @param SEED seed
#' @param radius node pie radius
#' @param width_range range of width of edges
#' @param text_size size of node labels
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph ggraph geom_edge_link scale_edge_width_continuous
#' @importFrom scatterpie geom_scatterpie
#' @importFrom modeltools posterior
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' cosine_network_term(SeuratObj, cosine_cal_by = "Topic modeling", pie_by = "Topic modeling", GSEA_by = "GO",
#'    topn = 10, layout = "fr", cos_sim_thresh = 0.8, radius = 0.1, text_size = 2)
#' cosine_network_term(SeuratObj, cosine_cal_by = "GSEA result", pie_by = "GSEA result", GSEA_by = "GO")
#' # edge(cosine similarity calculation) and node pie information could be different.
#' cosine_network_term(SeuratObj, cosine_cal_by = "GSEA result", pie_by = "Topic modeling", GSEA_by = "GO")
#' }
#'
#'
cosine_network_term <- function(SeuratObj,
                                cosine_cal_by = "Topic modeling",
                                pie_by = "Topic modeling",
                                GSEA_by = "GO",
                                topn = 10,
                                layout = "fr",
                                cos_sim_thresh = 0.8,
                                SEED = 123,
                                radius = 0.1,
                                width_range = c(0.05, 0.55),
                                text_size = 2)
{
  cosine_cal_by <- match.arg(cosine_cal_by, choices = c("Topic modeling", "GSEA result"))
  pie_by <- match.arg(pie_by, choices = c("Topic modeling", "GSEA result"))

  # calculate cosine similarity
  if (cosine_cal_by == "Topic modeling") { # calculate cosine similarity with topic~term probability matrix
    ldaOut <- slot(object = SeuratObj, name = 'misc')[["ldaOut"]]
    mmc <- posterior(ldaOut)$terms  # topic~term probability matrix
    # show top terms of every topic
    betaDF <- tidytext::tidy(ldaOut, matrix = "beta")
    tt <- betaDF %>% dplyr::group_by(topic) %>% slice_max(order_by = beta, n=topn, with_ties=F) %>% dplyr::pull(term) %>% unique()
  } else { # calculate cosine similarity with GSEA result
    # show top terms of every cluster
    GSEA_by <- match.arg(GSEA_by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
    GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", GSEA_by)]]
    GSEAresult %<>% dplyr::mutate(FDR=-log10(p.adjust))
    tt <- GSEAresult %>% dplyr::group_by(cluster) %>% slice_max(order_by = FDR, n=topn, with_ties=F) %>% dplyr::pull(ID) %>% unique()
    mmc <- reshape2::acast(GSEAresult, cluster~ID, value.var = "FDR")
    mmc[is.na(mmc)] <- 0
  }
  mmc <- mmc[, tt]
  ccc <- t(combn(colnames(mmc), 2))
  ccc <- plyr::mdply(ccc, .fun = function(x,y){
    x=mmc[, x]
    y=mmc[, y]
    sum(x*y)/(sqrt(sum(x^2))*sqrt(sum(y^2)))
  })
  links <- ccc %>% dplyr::filter(V1 >= cos_sim_thresh) %>% setNames(c('source', 'target', 'weight')) %>% dplyr::mutate(width=weight)
  nodes <- data.frame(name=unique(c(links$source, links$target)), stringsAsFactors = F)
  nodes$descrip <- GOfuncR::get_names(nodes$name)$go_name

  # node pie information
  if (pie_by == "Topic modeling") {
    ldaOut <- slot(object = SeuratObj, name = 'misc')[["ldaOut"]]
    mmp <- posterior(ldaOut)$terms  # topic~term probability matrix
  } else if (pie_by == "GSEA result") {
    GSEA_by <- match.arg(GSEA_by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
    GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", GSEA_by)]]
    GSEAresult %<>% dplyr::mutate(FDR=-log10(p.adjust))
    mmp <- reshape2::acast(GSEAresult, cluster~ID, value.var = "FDR")
    mmp[is.na(mmp)] <- 0
  }
  mmp <- mmp[, tt]
  clusColor <- setNames(scPalette2(nrow(mmp)), rownames(mmp))
  if (cosine_cal_by == pie_by) { # 饼图展示和余弦相似度计算一致的话则显示edge的颜色
    links$color <- apply(links, 1, function(x){names(which.max(c(mmp[, x[1]], mmp[, x[2]])))})
    links$color <- grDevices::adjustcolor(clusColor[links$color], 0.5)
  }
  gg <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  set.seed(SEED)
  p <- ggraph(gg, layout = layout)
  if (cosine_cal_by == pie_by) {
    p <- p + geom_edge_link(alpha=0.5, aes_(width=~width, edge_colour=~I(color))) + scale_edge_width_continuous(range = width_range)
  } else {
    p <- p + geom_edge_link(alpha=0.5, aes_(width=~width)) + scale_edge_width_continuous(range = width_range)
  }
  ## Get the matrix data for the pie plot
  ID_Cluster_mat <- as.data.frame(t(mmp[, nodes$name]))
  ID_Cluster_mat <- cbind(ID_Cluster_mat[p$data$name, ], p$data[, c("x", "y")])
  ID_Cluster_mat$radius <- radius
  p <- p + ggnewscale::new_scale_fill() + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                                                          cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA, alpha=1, sorted_by_radius=T) +
    scale_fill_manual(name=ifelse(pie_by == "Topic modeling", "topic", "cluster"), values=clusColor) +
    coord_fixed() + theme_void()# + theme(legend.position = "none")
  p <- p + ggrepel::geom_text_repel(aes_(x =~ x, y =~ y, label =~ descrip), size = text_size, show.legend = FALSE,
                                    segment.size = 0.3, force=3,force_pull=2, max.overlaps=30)
  return(p)
}


#' Show topic activity on UMAP/TSNE cell map
#'
#' @param SeuratObj Seurat object
#' @param reduction umap or tsne
#' @param topic which topic to show
#' @param pointSize 0.1 default.
#'
#' @importFrom modeltools posterior
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' topicProb(SeuratObj, reduction="umap", topic=3, pointSize=0.1)
#' }
#'
topicProb <- function(SeuratObj, reduction="umap", topic=1, pointSize=0.1) {
  if (!reduction %in% names(SeuratObj@reductions)) {
    stop("No '", reduction, "' reduction in Seurat object!")
  }
  ldaOut <- SeuratObj@misc$ldaOut
  em <- SeuratObj@meta.data
  em <- cbind(em, magrittr::set_colnames(Seurat::Embeddings(SeuratObj, reduction = reduction)[rownames(em), 1:2],
                                         c('Component_1', 'Component_2')))
  mm <- posterior(ldaOut)$topics  # cluster~topic probability matrix
  em$gamma <- mm[, as.character(topic)][as.character(Seurat::Idents(SeuratObj))]
  # scatter plot
  pp <- ggplot(em, aes(x=Component_1, y=Component_2, color=gamma)) +
    geom_point(size=pointSize) + theme_classic() +
    labs(x='Component_1', y='Component_2', color="topic probability", title = paste0("Topic ", topic)) +
    scale_color_gradientn(colours = rev(pals::brewer.piyg(10)))
  return(pp)
}



