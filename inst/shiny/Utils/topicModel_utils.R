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


plotindex <- function(ldas) {
  n <- length(ldas$perp)+1
  est=data.frame(x=rep(2:n, 3), y=c(ldas$perp, ldas$loglik, ldas$cosine),
                 type=rep(c('perplexity', 'loglikelihood', 'cosine similarity'), each=n-1))
  ggplot(est, aes(x, y, color=type)) +
    geom_line() + geom_point()+ scale_x_continuous(breaks = 0:n) +
    facet_wrap(~ type, scales = "free", ncol=3)+
    labs(x = "number of topic")
}


runLDA <- function(SeuratObj, by = "GO", k=NULL, method='VEM', SEED=1234, plot=TRUE) {
  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  #create a DocumentTermMatrix
  dd <- GSEAresult %>% dplyr::mutate(score = as.integer(-log2(p.adjust)))
  dtm <- tidytext::cast_dtm(data = dd, term = ID, document = cluster, value = score)
  rowTotals<- apply(dtm,1,sum) #Find the sum of words in each Document
  dtm <- dtm[rowTotals>0,]
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

# show top terms of every topic
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


# # echart4r network
# topicNW3 <- function(betaDF, topn) {
#   df <- betaDF %>% dplyr::group_by(topic) %>% slice_max(order_by = beta, n=topn, with_ties=F)
#   edgedf <- df %>% dplyr::mutate(topic=paste0("Topic ", topic)) %>% dplyr::select(topic, everything()) %>%
#     setNames(c("source", "target", "weight")) %>%
#     dplyr::mutate(width=2*(weight-min(weight))/(max(weight)-min(weight))+1)
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

# # echart4r network
# tooltips show description of ID
topicNW3 <- function(betaDF, topn, pws) {
  df <- betaDF %>% dplyr::group_by(topic) %>% slice_max(order_by = beta, n=topn, with_ties=F)
  edgedf <- df %>% dplyr::mutate(topic=paste0("Topic ", topic)) %>% dplyr::select(topic, everything()) %>%
    setNames(c("source", "target", "weight")) %>%
    dplyr::mutate(# width=2*(weight-min(weight))/(max(weight)-min(weight))+1, # 这样好像不对，会出现NaN
                  target=paste(target, unname(pws[target]), sep=": "))  # 适用于多种情况 GO/KEGG
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


wordcloud_topic <- function(betaDF, pws, topic, topn=20) {
  pp <- betaDF %>% filter(.data$topic == .env$topic) %>% slice_max(order_by = beta, n=topn, with_ties=F) %>% mutate(term=pws[term]) %>%
    ggplot(aes(label=term, size=beta, color=beta)) +
    geom_text_wordcloud_area(show.legend = T) +
    scale_color_gradientn(colours = rev(scPalette2(20))) + # scale_size(range = c(1,6)) +
    theme_minimal() + #ggtitle(paste0("Topic ", topic)) +
    theme(legend.position = "none")# +
    # theme(plot.margin = margin(t = 0,  # 顶部边缘距离
    #                            r = 0,  # 右边边缘距离
    #                            b = 0,  # 底部边缘距离
    #                            l = 0,  # 左边边缘距离
    #                            unit = "cm"))#设置单位为cm
  return(pp)
}

# wordcloud_topic <- function(betaDF, pws, topic, topn=20) {
#   pal2 <- colorRampPalette(brewer.pal(9,"Set1"))(20)
#   dd <- betaDF %>% filter(.data$topic == .env$topic) %>% dplyr::slice_max(order_by = beta, n=topn, with_ties=F)
#   suppressWarnings(wordcloud(pws[dd$term],dd$beta, scale=c(2,.1), min.freq=0,
#             max.words=20, random.order=F, rot.per=.15, colors=pal2, ordered.colors=T))
#   # # JS 3D wordcloud
#   # wordcloud2(dd %>% mutate(descrip=unname(pws[term])) %>% dplyr::select(descrip, beta), size = 0.1,color= "random-dark")
# }


# sankey
plot_sankey <- function(gammaDF, topn=1) {
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
                fontSize = 12, nodeWidth = 30, height = 600, width = 800, sinksRight=F)
}


# cluster做umap降维
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
    guides(color=guide_legend(title = "topic", override.aes = list(size=2))) +
    # geom_text(inherit.aes = T, aes(label=cluster), size=3) + #coord_fixed() +
    ggrepel::geom_text_repel(aes(label=cluster), max.overlaps = 20, segment.colour="black") +
    ggtitle("UMAP on cluster-topic probability matrix")
  return(pp)
}

# heatmap of cosine similarity between topics
cosineheatmap <- function(ldaOut) {
  mm <- posterior(ldaOut)$terms  # topic~term probability matrix
  cc <- suppressMessages(philentropy::distance(mm, method = "cosine"))
  rownames(cc) <- colnames(cc) <- rownames(mm)
  pheatmap::pheatmap(cc, color = colorRampPalette(rev(c("#000000", "#9932CC", "#EF3B2C","#FFFFCC")))(100), angle_col = "0",
                     fontsize = 9, main = "cosine similarity between topics")#, cellwidth = 30, cellheight = 30)
}

# network of cosine similarity between clusters
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

# cosine similarity between terms can be calculated by  topic~term probability matrix or GSEA result
# network of cosine similarity between terms
cosine_network_term <- function(SeuratObj,
                                cosine_cal_by = "Topic modeling",
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
  ldaOut <- slot(object = SeuratObj, name = 'misc')[["ldaOut"]]
  mm <- posterior(ldaOut)$terms  # topic~term probability matrix
  if (cosine_cal_by == "Topic modeling") { # calculate cosine similarity with topic~term probability matrix
    # show top terms of every topic
    betaDF <- tidytext::tidy(ldaOut, matrix = "beta")
    tt <- betaDF %>% dplyr::group_by(topic) %>% slice_max(order_by = beta, n=topn, with_ties=F) %>% dplyr::pull(term) %>% unique()
    mmc <- mm
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

  mm <- mm[, tt]
  clusColor <- setNames(scPalette2(nrow(mm)), rownames(mm))
  if (cosine_cal_by == "Topic modeling") {
    links$color <- apply(links, 1, function(x){names(which.max(c(mm[, x[1]], mm[, x[2]])))})
    links$color <- grDevices::adjustcolor(clusColor[links$color], 0.5)
  }
  gg <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  set.seed(SEED)
  p <- ggraph(gg, layout = layout)
  if (cosine_cal_by == "Topic modeling") {
    p <- p + geom_edge_link(alpha=0.5, aes_(width=~width, edge_colour=~I(color))) + scale_edge_width_continuous(range = width_range)
  } else {
    p <- p + geom_edge_link(alpha=0.5, aes_(width=~width)) + scale_edge_width_continuous(range = width_range)
  }
  ## Get the matrix data for the pie plot
  ID_Cluster_mat <- as.data.frame(t(mm[, nodes$name]))
  ID_Cluster_mat <- cbind(ID_Cluster_mat[p$data$name, ], p$data[, c("x", "y")])
  ID_Cluster_mat$radius <- radius
  p <- p + ggnewscale::new_scale_fill() + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                                                          cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA, alpha=1, sorted_by_radius=T) +
    scale_fill_manual(name="topic", values=clusColor) +
    coord_fixed() + theme_void()# + theme(legend.position = "none")
  p <- p + ggrepel::geom_text_repel(aes_(x =~ x, y =~ y, label =~ descrip), size = text_size, show.legend = FALSE,
                                    segment.size = 0.3, force=3,force_pull=2, max.overlaps=30)
  return(p)
}





