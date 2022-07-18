
gseaHeatmap1 <- reactive({
  gseaHeatmap(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, toshow = input$toshow, topPath = input$topPath,
              colour = input$hmpcol, scale = input$scale, fontsize_row = input$fontsize_row)
})
# output$gseaheatmap <- renderPlot(gseaHeatmap1())
output$gseaheatmap <- renderPlot({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if (nn %in% names(SeuratObj@misc)) {
    gseaHeatmap1()
  } else {
    return(NULL)
  }
})


# output$gseaheatmap1 <- renderUI({
#   plotOutput('gseaheatmap', height = input$plotheight) %>% withSpinner()
# })



clustercorplot1 <- reactive({
  if (input$similarity == "Jaccard index") {
    clustercorplot_jaccard(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, vertex.size.cex=input$vertex.size.cex,
                           vertex.label.cex=input$vertex.label.cex, edge.max.width=input$edge.max.width,
                           vertex.label.dist=input$vertex.label.dist, link_threshold=input$link_threshold_clustercorplot)
  } else {
    clustercorplot(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, vertex.size.cex=input$vertex.size.cex,
                   vertex.label.cex=input$vertex.label.cex, edge.max.width=input$edge.max.width,
                   vertex.label.dist=input$vertex.label.dist, link_threshold=input$link_threshold_clustercorplot)
  }
})

# output$clustercorplot <- renderPlot(clustercorplot1())
output$clustercorplot <- renderPlot({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if (nn %in% names(SeuratObj@misc)) {
    clustercorplot1()
  } else {
    return(NULL)
  }
})


circleplot1 <- reactive({
  circleplot(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, pvaluecutoff = input$pvaluecutoff,
             link_threshold = input$link_threshold_circleplot)
})

output$circleplot <- renderPlot({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if (nn %in% names(SeuratObj@misc)) {
    circleplot1()
  } else {
    return(NULL)
  }
})

embeddedplot1 <- reactive({
  embeddedplot(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, topaths = input$topaths,
               reduction = input$reduction, type = input$type, pie.size.cex = input$pie.size.cex)
})

output$embeddedplot <- renderPlot({
  if (paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) {
    embeddedplot1()
  } else {
    return(NULL)
  }
})


hierarchyplot_tree1 <- reactive({
  hierarchyplot_tree(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, topaths = input$topath1,
                     # cluster_cutree_k = input$cluster_cutree_k, pathway_cutree_k = input$pathway_cutree_k,
                     vertex.size.cex = input$vertex.size.cex2,
                     vertex.label.cex = input$vertex.label.cex2,
                     edge.max.width = input$edge.max.width2, alpha.edge = input$alpha.edge)
})

# output$hierarchyplot_tree <- renderPlot(hierarchyplot_tree1())
output$hierarchyplot_tree <- renderPlot({
  if (paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) {
    hierarchyplot_tree1()
  } else {
    return(NULL)
  }
})

emapplotPie1 <- reactive({
  emapplotPie(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, showCategory = input$showCategory,
              node_label_cex = input$node_label_cex,
              node_size_cex = input$node_size_cex, layout = input$layout,
              cex_line = input$cex_line, pie = input$pie)
})

# output$emapplotPie <- renderPlot(emapplotPie1())
output$emapplotPie <- renderPlot({
  if (paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) {
    emapplotPie1()
  } else {
    return(NULL)
  }
})

emapplot1 <- reactive({
  emapplot2(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, showCategory = input$showCategory2,
            node_label_cex = input$node_label_cex2, cluster = input$cluster,
            node_size_cex = input$node_size_cex2, layout = input$layout2,
            cex_line = input$cex_line2)
})

# output$emapplot <- renderPlot(emapplot1())
output$emapplot <- renderPlot({
  if (paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) {
    emapplot1()
  } else {
    return(NULL)
  }
})


goplot1 <- reactive({
  goplot2(SeuratObj, cluster = input$cluster2, pathwayIDs = inter$pwss, showCategory = input$showCategory3,
          ont = input$ont, label_size = input$label_size)
})
# output$goplot <- renderPlot(goplot1())
output$goplot <- renderPlot({
  if ((paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) & (input$switchTerm == "GO")) {
    goplot1()
  } else {
    return(NULL)
  }
})


# simplifyEnrichmentplot1 <- reactive({
#   simplifyEnrichmentplot(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, showCategory = input$showCategory4, GO_ont = input$ont2)
# })
# # output$goplot <- renderPlot(goplot1())
# output$simplifyEnrichmentplot <- renderPlot({
#   if ((paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) & (input$switchTerm %in% c("GO", "KEGG", "Reactome", "MSigDb"))) {
#     simplifyEnrichmentplot1()
#   } else {
#     return(NULL)
#   }
# })


