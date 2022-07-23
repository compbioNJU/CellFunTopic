

# observeEvent(input$sideBarTab, {
observe({
  if ("ldaOut" %in% names(ldaDF()$SeuratObj@misc)) {
    # shinyjs::hide(id = "tm_hide1")
    shinyjs::show(id = "tm_hide2")
  } else {
    # shinyjs::show(id = "tm_hide1")
    shinyjs::hide(id = "tm_hide2")
  }

  if ("ldas" %in% names(ldaDF()$SeuratObj@misc)) {
    shinyjs::show(id = "tm_hide3")
  } else {
    shinyjs::hide(id = "tm_hide3")
  }
})


ldaDF <- reactive({
  shinyjs::hide(id = "tm_hide1")
  if (!("ldaOut" %in% names(SeuratObj@misc))) {
    cat('No topic modeling result were found in SeuratObj@misc, topic modeling is running now.\n')
    withProgress(value=0.5,
                 message="Running topic modeling",
                 detail="please wait for a while ...",
                 {
                   SeuratObj <<- runLDA(SeuratObj, by = "GO", k=NULL, method='VEM', SEED=1234, plot=F)
                   incProgress(0.5)
                   Sys.sleep(0.1)
                 })
    shinyjs::show(id = "tm_hide1")
  }
  if ("ldas" %in% names(SeuratObj@misc)) {
    ldas <- SeuratObj@misc$ldas
  } else {ldas <- NULL}
  ldaOut <- SeuratObj@misc$ldaOut
  betaDF <- tidytext::tidy(ldaOut, matrix = "beta")
  gammaDF <- tidytext::tidy(ldaOut, matrix = "gamma")
  dd <- SeuratObj@misc$GSEAresult_GO %>% dplyr::mutate(score = as.integer(-log2(p.adjust)))
  dtm <- tidytext::cast_dtm(data = dd, term = ID, document = cluster, value = score)
  rowTotals <- apply(dtm, 1, sum)
  colTotals <- apply(dtm, 2, sum)
  dtm <- dtm[rowTotals > 0, colTotals > 0]
  assignws <- tidytext::augment(ldaOut, data = dtm)
  return(list(SeuratObj=SeuratObj, ldas=ldas, ldaOut=ldaOut, betaDF=betaDF, gammaDF=gammaDF, assignws=assignws))
})


output$tm_echt <- renderEcharts4r({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  betaDF <- ldaDF()$betaDF
  pws <- ID2Description(SeuratObj, by = "GO")
  topicNW3(betaDF, topn=input$tm_topn1, pws=pws)
})


tm_hist3 <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  betaDF <- ldaDF()$betaDF
  pws <- ID2Description(SeuratObj, by = "GO")
  betaDF %<>% dplyr::mutate(descrip=unname(pws[term]))
  # betaDF %<>% dplyr::mutate(descrip=GOfuncR::get_names(term)$go_name) # 只适用于GO的情况
  tt <- gsub("Topic ", "", input$tm_echt_clicked_data$source)
  if (length(tt)<1) tt <- unique(betaDF$topic)[1]
  hist_topic_term(betaDF, topics=tt, topn=input$tm_topn1, axis.text.y.size=10) + theme(legend.position = "none")
})

output$tm_hist3 <- renderPlot({
  tm_hist3()
})

tm_wordcloud <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  tt <- gsub("Topic ", "", input$tm_echt_clicked_data$source)
  betaDF <- ldaDF()$betaDF
  if (length(tt)<1) tt <- unique(betaDF$topic)[1]
  pws <- ID2Description(SeuratObj, by = "GO")
  wordcloud_topic(betaDF, pws, topic=tt, topn=input$tm_topn1)
})

output$tm_wordcloud <- renderPlot({
  tm_wordcloud()
})

# output$tm_wordcloud <- renderWordcloud2({
#   tm_wordcloud()
# }) # 3D词云，但是会导致其他的图比如网络图加载不出来


betatable <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  betaDF <- ldaDF()$betaDF
  pws <- ID2Description(SeuratObj, by = "GO") # 这样做适用于GO、KEGG等等多种情况
  betaDF %<>% dplyr::mutate(descrip=unname(pws[term])) %>% dplyr::rename(probability=beta)
  # GOfuncR::get_names(term)$go_name # 只适用于GO的情况

  DT::datatable(
    data = betaDF,
    class = 'display',
    rownames = FALSE,
    filter = "top",
    plugins = "natural",
    extensions = 'Buttons',
    options = list(
      scrollCollapse = T,
      autoWidth = T,
      paging = T,
      pageLength = 10,
      fixedHeader = T,
      scrollX = T,
      scrollY = T,
      searchHighlight = TRUE,
      dom = "Bfrtip",
      fnDrawCallback = htmlwidgets::JS("
            function() {
              HTMLWidgets.staticRender();
            }
          "),
      buttons = list(
        list(
          extend = "collection",
          text = "Download",
          buttons = list(
            list(
              extend = "csv",
              filename = paste0(input$childparent, "Terms"),
              title = paste0(input$childparent, "Terms")
            ),
            list(
              extend = "excel",
              filename = paste0(input$childparent, "Terms"),
              title = paste0(input$childparent, "Terms")
            )
          )
        )
      ),
      columnDefs = list(
        list(
          className = "dt-center",
          targets = c(0,1)
        ),
        list(
          className = "dt-left",
          targets = c(2,3)
        ),
        list(width = "200px", targets = c(3))
      )
    )
  ) %>%
    formatStyle(
      columns = "topic",
      color = styleEqual(
        unique(betaDF$topic),
        scPalette(length(unique(betaDF$topic)))
      ),
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "probability",
      color = styleInterval(seq(min(betaDF$probability), max(betaDF$probability), length.out = 49),
                            colorRampPalette(c("#807DBA", "#3F007D"))(50)),
      background = styleColorBar(betaDF$probability,
                                 "#6A51A3", angle = -90),
      backgroundSize = "98% 18%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "bottom",
      fontWeight = "bold"
    )
})

output$betatable <- DT::renderDataTable(
  betatable()
)

topicsChoiceui <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  betaDF <- ldaDF()$betaDF
  selectInput(
    inputId = "tm_topics",
    label = "topics to show",
    choices = unique(betaDF$topic),
    selected = head(unique(betaDF$topic), 4),
    multiple = T
  )
})
output$topicsChoiceui <- renderUI(topicsChoiceui())


tm_hist1 <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  betaDF <- ldaDF()$betaDF
  pws <- ID2Description(SeuratObj, by = "GO")
  betaDF %<>% dplyr::mutate(descrip=unname(pws[term]))
  # betaDF %<>% dplyr::mutate(descrip=GOfuncR::get_names(term)$go_name) # 只适用于GO的情况
  hist_topic_term(betaDF, topics=input$tm_topics, topn=input$tm_topn2, axis.text.y.size=input$tm_size)
})

output$tm_hist1 <- renderPlot({
  tm_hist1()
})







gammatable <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  gammaDF <- ldaDF()$gammaDF %>% dplyr::rename(probability=gamma, cluster=document)

  DT::datatable(
    data = gammaDF,
    class = 'display',
    rownames = FALSE,
    filter = "top",
    plugins = "natural",
    extensions = 'Buttons',
    options = list(
      scrollCollapse = T,
      autoWidth = T,
      paging = T,
      pageLength = 10,
      fixedHeader = T,
      scrollX = T,
      scrollY = T,
      searchHighlight = TRUE,
      dom = "Bfrtip",
      fnDrawCallback = htmlwidgets::JS("
            function() {
              HTMLWidgets.staticRender();
            }
          "),
      buttons = list(
        list(
          extend = "collection",
          text = "Download",
          buttons = list(
            list(
              extend = "csv",
              filename = paste0(input$childparent, "Terms"),
              title = paste0(input$childparent, "Terms")
            ),
            list(
              extend = "excel",
              filename = paste0(input$childparent, "Terms"),
              title = paste0(input$childparent, "Terms")
            )
          )
        )
      ),
      columnDefs = list(
        list(
          className = "dt-center",
          targets = c(0,1)
        ),
        list(
          width = "100px",
          className = "dt-left",
          targets = 2
        )
      )
    )
  ) %>%
    formatStyle(
      columns = "topic",
      color = styleEqual(
        unique(gammaDF$topic),
        scPalette(length(unique(gammaDF$topic)))
      ),
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "cluster",
      color = styleEqual(
        unique(gammaDF$cluster),
        scPalette(length(unique(gammaDF$cluster)))
      ),
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "probability",
      color = styleInterval(seq(min(gammaDF$probability), max(gammaDF$probability), length.out = 49),
                            colorRampPalette(c("#807DBA", "#3F007D"))(50)),
      background = styleColorBar(gammaDF$probability,
                                 "#6A51A3", angle = -90),
      backgroundSize = "98% 18%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "bottom",
      fontWeight = "bold"
    )
})

output$gammatable <- DT::renderDataTable(
  gammatable()
)

tm_clustersChoiceui <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  gammaDF <- ldaDF()$gammaDF
  selectInput(
    inputId = "tm_clusters",
    label = "clusters to show",
    choices = unique(gammaDF$document),
    selected = head(unique(gammaDF$document), 3),
    multiple = T
  )
})
output$tm_clustersChoiceui <- renderUI(tm_clustersChoiceui())


tm_hist2 <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  gammaDF <- ldaDF()$gammaDF
  hist_cluster_topic(gammaDF, clusters=input$tm_clusters)
})

output$tm_hist2 <- renderPlot({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  tm_hist2()
})




assignwstable <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  pws <- ID2Description(SeuratObj, by = "GO")
  assignws <- ldaDF()$assignws %>% dplyr::mutate(description=unname(pws[term])) %>%
    dplyr::rename(cluster=document, topic=.topic) %>% dplyr::select(cluster, term, description, everything())

  DT::datatable(
    data = assignws,
    class = 'display',
    rownames = FALSE,
    filter = "top",
    plugins = "natural",
    extensions = 'Buttons',
    options = list(
      scrollCollapse = T,
      autoWidth = T,
      paging = T,
      pageLength = 10,
      fixedHeader = T,
      scrollX = T,
      scrollY = T,
      searchHighlight = TRUE,
      dom = "Bfrtip",
      fnDrawCallback = htmlwidgets::JS("
            function() {
              HTMLWidgets.staticRender();
            }
          "),
      buttons = list(
        list(
          extend = "collection",
          text = "Download",
          buttons = list(
            list(
              extend = "csv",
              filename = paste0(input$childparent, "Terms"),
              title = paste0(input$childparent, "Terms")
            ),
            list(
              extend = "excel",
              filename = paste0(input$childparent, "Terms"),
              title = paste0(input$childparent, "Terms")
            )
          )
        )
      ),
      columnDefs = list(
        list(
          className = "dt-center",
          targets = c(0,1,3,4)
        ),
        list(
          width = "200px",
          className = "dt-left",
          targets = 2
        )
      )
    )
  ) %>%
    formatStyle(
      columns = "topic",
      color = styleEqual(
        unique(assignws$topic),
        scPalette(length(unique(assignws$topic)))
      ),
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "cluster",
      color = styleEqual(
        unique(assignws$cluster),
        scPalette(length(unique(assignws$cluster)))
      ),
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "count",
      color = styleInterval(seq(min(assignws$count), max(assignws$count), length.out = 49),
                            colorRampPalette(c("#807DBA", "#3F007D"))(50)),
      background = styleColorBar(assignws$count,
                                 "#6A51A3", angle = -90),
      backgroundSize = "98% 18%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "bottom",
      fontWeight = "bold"
    )
})

output$assignwstable <- DT::renderDataTable(
  assignwstable()
)


tm_sankey <- reactive({
  gammaDF <- ldaDF()$gammaDF
  plot_sankey(gammaDF, topn=input$tm_topn3)
})

output$tm_sankey <- renderSankeyNetwork(
  tm_sankey()
)

tm_umap_cluster <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  ldaOut <- ldaDF()$ldaOut
  umap_cluster(ldaOut)
})

output$tm_umap_cluster <- renderPlot(
  tm_umap_cluster()
)

tm_heatmap <- reactive({
  # # probability between topic and cluster, heatmap
  ldaOut <- ldaDF()$ldaOut
  cluster_topic_hmp(ldaOut)
})

output$tm_heatmap <- renderPlot(
  tm_heatmap()
)


tm_index <- reactive({
  # if (!("ldas" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  # ldas <- slot(object = SeuratObj, name = 'misc')[["ldas"]]
  if (is.null(ldaDF()$ldas)) return(NULL)
  ldas <- ldaDF()$ldas
  plotindex(ldas)
})

output$tm_index <- renderPlot(
  tm_index()
)


tm_cosineheatmap <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  ldaOut <- ldaDF()$ldaOut
  cosineheatmap(ldaOut)
})

output$tm_cosineheatmap <- renderPlot(
  tm_cosineheatmap()
)


tm_cosine_network_cluster <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  ldaOut <- ldaDF()$ldaOut
  cosine_network_cluster(ldaOut,
                         layout=input$tm_layout,
                         cos_sim_thresh=input$cos_sim_thresh,
                         SEED=123,
                         radius=input$radiuscosnet,
                         width_range=input$width_range)
})

output$tm_cosine_network_cluster <- renderPlot(
  tm_cosine_network_cluster()
)


tm_cosine_network_term <- reactive({
  if (!("ldaOut" %in% names(ldaDF()$SeuratObj@misc))) return(NULL)
  cosine_network_term(ldaDF()$SeuratObj,
                      cosine_cal_by = input$cosine_cal_by,
                      topn = input$tm_cosnet_topn,
                      layout = input$tm_layout2,
                      cos_sim_thresh = input$cos_sim_thresh2,
                      SEED = 123,
                      radius = input$radiuscosnet2,
                      width_range = input$width_range2,
                      text_size = input$tm_text_size)
})

output$tm_cosine_network_term <- renderPlot(
  tm_cosine_network_term()
)




######## download
output$tm_sankeydl <- downloadHandler(
  filename = "sankey.html",
  content = function(file) {
    gammaDF <- ldaDF()$gammaDF
    plot_sankey(gammaDF) %>% saveNetwork(file = file)
  },
  contentType = "html"
)

output$tm_umap_clusterdl <- downloadHandler(
  filename = "umap_cluster.pdf",
  content = function(file) {
    pdf(file, width = 10, height = 10)
    ldaOut <- ldaDF()$ldaOut
    pp <- umap_cluster(ldaOut)
    print(pp)
    dev.off()
  },
  contentType = "image/pdf"
)

output$tm_heatmapdl <- downloadHandler(
  filename = "cluster_topic_heatmap.pdf",
  content = function(file) {
    pdf(file, width = 10, height = 15)
    ldaOut <- ldaDF()$ldaOut
    cluster_topic_hmp(ldaOut)
    dev.off()
  },
  contentType = "image/pdf"
)

output$tm_indexdl <- downloadHandler(
  filename = "metrics.pdf",
  content = function(file) {
    pdf(file, width = 15, height = 6)
    ldas <- ldaDF()$ldas
    pp <- plotindex(ldas)
    print(pp)
    dev.off()
  },
  contentType = "image/pdf"
)

output$tm_cosineheatmapdl <- downloadHandler(
  filename = "cosineheatmap.pdf",
  content = function(file) {
    pdf(file, width = 10, height = 10)
    ldaOut <- ldaDF()$ldaOut
    cosineheatmap(ldaOut)
    dev.off()
  },
  contentType = "image/pdf"
)

output$tm_cosine_network_clusterdl <- downloadHandler(
  filename = "cosine_network_cluster.pdf",
  content = function(file) {
    pdf(file, width = 13, height = 13)
    ldaOut <- ldaDF()$ldaOut
    pp <- cosine_network_cluster(ldaOut,
                                 layout=input$tm_layout,
                                 cos_sim_thresh=input$cos_sim_thresh,
                                 SEED=123,
                                 radius=input$radiuscosnet,
                                 width_range=input$width_range)
    print(pp)
    dev.off()
  },
  contentType = "image/pdf"
)

output$tm_cosine_network_termdl <- downloadHandler(
  filename = "cosine_network_term.pdf",
  content = function(file) {
    pdf(file, width = 13, height = 13)
    pp <- cosine_network_term(ldaDF()$SeuratObj,
                              cosine_cal_by = input$cosine_cal_by,
                              topn = input$tm_cosnet_topn,
                              layout = input$tm_layout2,
                              cos_sim_thresh = input$cos_sim_thresh2,
                              SEED = 123,
                              radius = input$radiuscosnet2,
                              width_range = input$width_range2,
                              text_size = input$tm_text_size)
    print(pp)
    dev.off()
  },
  contentType = "image/pdf"
)



