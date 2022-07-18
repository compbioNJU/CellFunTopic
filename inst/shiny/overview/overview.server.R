source(
  file = "overview/overview.plot.R",
  local = TRUE,
  encoding = "UTF-8"
)

GSEAresulttable <- reactive({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if (nn %in% names(SeuratObj@misc)) {
    dt2 <- slot(object = SeuratObj, name = 'misc')[[nn]]
    if (nn == "GSEAresult_GO") {
      dt2$QuickGO <- sprintf('<a href="https://www.ebi.ac.uk/QuickGO/term/%s"><i class=\"fa fa-external-link-alt\"></i> %s</a>', dt2$ID, dt2$ID)
      # dt2$QuickGO <- sprintf('<a href="https://www.ebi.ac.uk/QuickGO/term/%s" class=\"fa fa-external-link-alt\"> %s</a>', dt2$ID, dt2$ID)
    }
    dt2 %<>% dplyr::mutate(`-log10(FDR)`=-log10(p.adjust)) %>% dplyr::select("cluster","ID", "Description", "-log10(FDR)", everything())
  } else {
    return(NULL)
    }
  
  DT::datatable(
    data = dt2,
    class = 'display',
    escape = FALSE,
    rownames = FALSE,
    filter = "top",
    plugins = "natural",
    extensions = c("RowGroup","Buttons"),
    callback = htmlwidgets::JS(paste0(
      "table.rowGroup().",
      ifelse(input$tableShowSetting, "enable()", "disable()"),
      ".draw();"
    )),
    options = list(
      scrollCollapse = T,
      autoWidth = T,
      orderClasses = T,
      paging = T,
      pageLength = 10,
      rowGroup = list(dataSrc = 0),
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
          extend = "colvis",
          columns = c(3:14),
          text = "Columns to show"
        ),
        list(
          extend = "collection",
          text = "Download",
          buttons = list(
            list(
              extend = "csv",
              filename = "enriched_pathways",
              title = "Enriched pathways"
            ),
            list(
              extend = "excel",
              filename = "enriched_pathways",
              title = "Enriched pathways"
            )
          )
        )
      ),
      columnDefs = list(
        list(
          visible = F,
          targets = c(5,12,13)
        ),
        list(
          className = "dt-center", 
          targets = c(0,1,4:7,11,14)
        ),
        list(
          className = "dt-left", 
          targets = c(3,8:10)
        ),
        list(
          width = "220px",
          className = "dt-left",
          targets = 2
        ),
        list(width = "20px", targets = c(4,5)),
        list(width = "15px", targets = c(11)),
        list(width = "100px", targets = c(14))
      )
    )
  ) %>% 
    formatStyle(
    columns = "cluster",
    color = styleEqual(
      unique(dt2$cluster), 
      scPalette(length(unique(dt2$cluster)))
    ),
    fontWeight = "bold"
  ) %>% 
    formatStyle(
      columns = "ONTOLOGY",
      backgroundColor = styleEqual(
        unique(dt2$ONTOLOGY), 
        grDevices::adjustcolor(scPalette(length(unique(dt2$ONTOLOGY))),0.3)
      ),
      fontWeight = "bold"
    ) %>% 
    formatStyle(
    columns = 'p.adjust',
    background = styleColorBar(dt2$p.adjust, 
                               '#FC4E07', angle = -90),
    backgroundSize = '100% 50%',
    backgroundRepeat = 'no-repeat',
    backgroundPosition = 'center',
    fontWeight = "bold"
  ) %>% 
    formatStyle(
    columns = 'pvalue',
    background = styleColorBar(dt2$pvalue, 
                               '#00AFBB', angle = -90),
    backgroundSize = '100% 50%',
    backgroundRepeat = 'no-repeat',
    backgroundPosition = 'center',
    fontWeight = "bold"
  ) %>% 
    formatStyle(
      columns = 'qvalues',
      background = styleColorBar(dt2$qvalues,
                                 "#F39C11", angle = -90),
      backgroundSize = '100% 50%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center',
      fontWeight = "bold"
    ) %>%
    # formatStyle(
    #   columns = "qvalues",
    #   color = styleInterval(seq(min(dt2$qvalues), max(dt2$qvalues), length.out = 49),
    #                         colorRampPalette(c("#F8BF76", "#DB8B0A"))(50)),
    #   background = styleColorBar(dt2$qvalues,
    #                              "#F39C11", angle = -90),
    #   backgroundSize = "98% 18%",
    #   backgroundRepeat = "no-repeat",
    #   backgroundPosition = "bottom",
    #   fontWeight = "bold"
    # ) %>%
    formatStyle(
      columns = "-log10(FDR)",
      color = styleInterval(seq(min(dt2$`-log10(FDR)`), max(dt2$`-log10(FDR)`), length.out = 49),
                            colorRampPalette(c("#807DBA", "#3F007D"))(50)),
      background = styleColorBar(dt2$`-log10(FDR)`,
                                 "#6A51A3", angle = -90),
      backgroundSize = "98% 18%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "bottom",
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "NES",
      backgroundColor = styleInterval(seq(min(dt2$NES), max(dt2$NES), length.out = 49),
                                      colorRampPalette(c("#FFFFFF", "#B03C2D"))(50)), 
      fontWeight = "bold"
      ) %>%
    formatStyle(
      columns = "enrichmentScore",
      backgroundColor = styleInterval(seq(min(dt2$enrichmentScore), max(dt2$enrichmentScore), length.out = 49),
                                      colorRampPalette(c("#FFFFFF", "#DB8B0A"))(50)), 
      fontWeight = "bold"
    ) %>% 
    formatSignif(
      columns = c("-log10(FDR)","enrichmentScore", "NES", "pvalue","p.adjust","qvalues"),
      digits = 3
    )
})

output$GSEAtable <- DT::renderDataTable(
  GSEAresulttable(),
  server = TRUE
)






output$grouping <- renderUI({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if (nn %in% names(SeuratObj@misc)) {
    tags$span(
      dropdownButton(
        tags$h4(icon("eye"), "Options"),
        tags$hr(),
        materialSwitch(
          inputId = "tableShowSetting",
          label = tagList(icon("object-group"), "Grouping"),
          status = "danger",
          value = T
        ),
        circle = F,
        right = T,
        inline = T,
        status = "danger",
        icon = icon("gear"),
        size = "sm",
        width = "300px",
        tooltip = tooltipOptions(title = "Options", placement = "top")
      ),
      style = "float:right;"
    )
  } else {
    return(NULL)
  }
})


clusterChoiceui <- reactive({
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", input$switchTerm)]]
  clus <- unique(GSEAresult$cluster)
  selectInput(
    inputId = "cluster",
    label = "cluster", 
    choices = clus
  )
})

output$clusterChoiceui <- renderUI(clusterChoiceui())


clusterChoiceui2 <- reactive({
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[["GSEAresult_GO"]]
  clus <- unique(GSEAresult$cluster)
  selectInput(
    inputId = "cluster2",
    label = "cluster", 
    choices = clus
  )
})
output$clusterChoiceui2 <- renderUI(clusterChoiceui2())



# 
# observeEvent(input$GSEAtable_rows_selected, {
#   output$pathwayList <- renderUI({
#     rids <- input$GSEAtable_rows_selected
#     GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", input$switchTerm)]]
#     nn <- paste(GSEAresult[rids, "ID"], GSEAresult[rids, "Description"], sep = " ")
#     dat <- GSEAresult[rids, "ID"]
#     if (!is.null(dat)) {
#       dat <- setNames(object = dat, nn)
#     }
#     tags$span(
#       pickerInput(
#         inputId = "selectedpw", 
#         label = "Selected:",
#         choices = dat, 
#         selected = GSEAresult[rids, "ID"],
#         options = list(
#           `actions-box` = TRUE, 
#           size = 10,
#           `selected-text-format`= "count",
#           width = "300px"
#         ),
#         multiple = TRUE,
#         inline = TRUE
#       ),
#       style = "float:right;"
#     )
#   })
# })

output$pathwayList <- renderUI({
  if (paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) {
    rids <- input$GSEAtable_rows_selected
    GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", input$switchTerm)]]
    nn <- paste(GSEAresult[rids, "ID"], GSEAresult[rids, "Description"], sep = " ")
    dat <- GSEAresult[rids, "ID"]
    if (!is.null(dat)) {
      dat <- setNames(object = dat, nn)
    }
    tags$span(
      pickerInput(
        inputId = "selectedpw",
        label = "Selected:",
        choices = dat,
        selected = GSEAresult[rids, "ID"],
        options = list(
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format`= "count",
          width = "300px"
        ),
        multiple = TRUE,
        inline = TRUE
      ),
      style = "float:right;"
    )
  } else {
    return(NULL)
  }
})






output$select_plotui2 <- renderUI({
  tagList(tags$b("No GSEA result of "), dashboardLabel(input$switchTerm, status = "info", style = "square"), tags$b("were found."))
})

observeEvent(input$switchTerm, {
  if (!(paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc))) {
    shinyjs::hide(id = "pa")
    shinyjs::hide(id = "tohide1")
    shinyjs::show(id = "tohide2")
  } else {
    shinyjs::show(id = "pa")
    shinyjs::show(id = "tohide1")
    shinyjs::hide(id = "tohide2")
  }
  if (input$switchTerm == "GO") {
    shinyjs::show(id = "tohide3")
  } else {
    shinyjs::hide(id = "tohide3")
  }
})



observeEvent(input$button1, {
  dataTableProxy('GSEAtable') %>% selectRows(input$GSEAtable_rows_all)
})
observeEvent(input$button2, {
  dataTableProxy('GSEAtable') %>% selectRows(NULL)
})


# plotui <- reactive({
#   if (paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) {
#     fluidRow(
#       fluidRow(
#         column(width = 12,
#                boxPlus(
#                  title = tagList(
#                    icon("expand-arrows-alt"),
#                    "Hierarchyplot"
#                  ),
#                  closable = F,
#                  collapsible = T,
#                  width = 12,
#                  hierarchyplot_treeui(),
#                  plotOutput('hierarchyplot_tree', height = "800px") %>% withSpinner()
#                )
#         )
#       ),
#       fluidRow(
#         column(
#           width = 6,
#           boxPlus(
#             title = tagList(
#               icon("th"),
#               "Heatmap"
#             ),
#             closable = F,
#             collapsible = T,
#             width = 12,
#             gseaheatmapui(),
#             plotOutput('gseaheatmap', height = "980px") %>% withSpinner()
#             # uiOutput('gseaheatmap1') %>% withSpinner()  # Warning较多
#           )
#         ),
#         column(
#           width = 6,
#           boxPlus(
#             title = tagList(
#               icon("chart-bar"),
#               "Embeddedplot"
#             ),
#             closable = F,
#             collapsible = T,
#             width = 12,
#             embeddedplotui(),
#             plotOutput('embeddedplot', height = "400px") %>% withSpinner()
#           ),
#           boxPlus(
#             title = tagList(
#               icon("dharmachakra"),
#               # icon("vector-path-circle"),
#               "Circleplot"
#             ),
#             closable = F,
#             collapsible = T,
#             footer = tags$small(icon("lightbulb"), "Link width represents number of intersection pathways between clusters,
#                             link color is the same as the cluster with higher-score intersection pathways."),
#             width = 12,
#             circleplotui(),
#             plotOutput('circleplot', height = "400px") %>% withSpinner()
#           )
#         )
#       ),
#       fluidRow(
#         column(
#           width = 6,
#           boxPlus(
#             title = tagList(
#               icon("connectdevelop"),
#               "Clustercorplot"
#             ),
#             closable = F,
#             collapsible = T,
#             footer = tags$small(icon("lightbulb"), "Helps to infer the relationship and similarity between clusters."),
#             width = 12,
#             clustercorplotui(),
#             plotOutput('clustercorplot', height = "400px") %>% withSpinner()
#           )
#         ),
#         column(
#           width = 6,
#           boxPlus(
#             title = tagList(
#               icon("bezier-curve"),
#               "EmapplotPie"
#             ),
#             closable = F,
#             collapsible = T,
#             width = 12,
#             emapplotPieui(),
#             plotOutput('emapplotPie', height = "400px") %>% withSpinner()
#           )
#         )
#       ),
#       fluidRow(
#         column(
#           width = 6,
#           boxPlus(
#             title = tagList(
#               icon("bezier-curve"),
#               "Emapplot"
#             ),
#             closable = F,
#             collapsible = T,
#             width = 12,
#             emapplotui(),
#             plotOutput('emapplot', height = "400px") %>% withSpinner()
#           )
#         ),
#         column(
#           width = 6,
#           uiOutput("goplot_ui") %>% withSpinner()
#         )
#       )#,
#       # fluidRow(
#       #   column(
#       #     width = 6,
#       #     uiOutput("simplifyEnrichmentplot_ui") %>% withSpinner()
#       #   ),
#       #   column(
#       #     width = 6
#       #   )
#       # )
#     )
#   } else {
#     return(NULL)
#   }
# })
# 
# output$plotui <- renderUI(plotui())

# output$plotui <- renderUI({
#   fluidPage(
#     fluidRow(
#       column(width = 12,
#              boxPlus(
#                title = tagList(
#                  icon("expand-arrows-alt"),
#                  "Hierarchyplot"
#                ),
#                closable = F,
#                collapsible = T,
#                width = 12,
#                hierarchyplot_treeui(),
#                plotOutput('hierarchyplot_tree', height = "800px") %>% withSpinner()
#              )
#       )
#     ),
#     fluidRow(
#       column(
#         width = 6,
#         boxPlus(
#           title = tagList(
#             icon("th"),
#             "Heatmap"
#           ),
#           closable = F,
#           collapsible = T,
#           width = 12,
#           gseaheatmapui(),
#           plotOutput('gseaheatmap', height = "980px") %>% withSpinner()
#           # uiOutput('gseaheatmap1') %>% withSpinner()  # Warning较多
#         )
#       ),
#       column(
#         width = 6,
#         boxPlus(
#           title = tagList(
#             icon("chart-bar"),
#             "Embeddedplot"
#           ),
#           closable = F,
#           collapsible = T,
#           width = 12,
#           embeddedplotui(),
#           plotOutput('embeddedplot', height = "400px") %>% withSpinner()
#         ),
#         boxPlus(
#           title = tagList(
#             icon("dharmachakra"),
#             # icon("vector-path-circle"),
#             "Circleplot"
#           ),
#           closable = F,
#           collapsible = T,
#           footer = tags$small(icon("lightbulb"), "Link width represents number of intersection pathways between clusters,
#                             link color is the same as the cluster with higher-score intersection pathways."),
#           width = 12,
#           circleplotui(),
#           plotOutput('circleplot', height = "400px") %>% withSpinner()
#         )
#       )
#     ),
#     fluidRow(
#       column(
#         width = 6,
#         boxPlus(
#           title = tagList(
#             icon("connectdevelop"),
#             "Clustercorplot"
#           ),
#           closable = F,
#           collapsible = T,
#           footer = tags$small(icon("lightbulb"), "Helps to infer the relationship and similarity between clusters."),
#           width = 12,
#           clustercorplotui(),
#           plotOutput('clustercorplot', height = "400px") %>% withSpinner()
#         )
#       ),
#       column(
#         width = 6,
#         boxPlus(
#           title = tagList(
#             icon("bezier-curve"),
#             "EmapplotPie"
#           ),
#           closable = F,
#           collapsible = T,
#           width = 12,
#           emapplotPieui(),
#           plotOutput('emapplotPie', height = "400px") %>% withSpinner()
#         )
#       )
#     ),
#     fluidRow(
#       column(
#         width = 6,
#         boxPlus(
#           title = tagList(
#             icon("bezier-curve"),
#             "Emapplot"
#           ),
#           closable = F,
#           collapsible = T,
#           width = 12,
#           emapplotui(),
#           plotOutput('emapplot', height = "400px") %>% withSpinner()
#         )
#       ),
#       column(
#         width = 6,
#         uiOutput("goplot_ui") %>% withSpinner()
#       )
#     )#,
#     # fluidRow(
#     #   column(
#     #     width = 6,
#     #     uiOutput("simplifyEnrichmentplot_ui") %>% withSpinner()
#     #   ),
#     #   column(
#     #     width = 6
#     #   )
#     # )
#   )
# })



# goplot_ui <- reactive({
#   if (input$switchTerm == "GO") {
#     boxPlus(
#       title = tagList(
#         icon("sitemap"),
#         "GOplot"
#       ),
#       closable = F,
#       collapsible = T,
#       width = 12,
#       goplotui(),
#       plotOutput('goplot', height = "400px") %>% withSpinner()
#     )
#   } else {
#     return(NULL)
#   }
# })
# output$goplot_ui <- renderUI(goplot_ui())

# 
# simplifyEnrichmentplot_ui <- reactive({
#   if (input$switchTerm %in% c("GO", "KEGG", "Reactome", "MSigDb")) {
#     boxPlus(
#       title = tagList(
#         icon("cloud"),
#         "simplifyEnrichment"
#       ),
#       closable = F,
#       collapsible = T,
#       width = 12,
#       simplifyEnrichmentplotui(),
#       plotOutput('simplifyEnrichmentplot') %>% withSpinner()
#     )
#   } else {
#     return(NULL)
#   }
# })
# output$simplifyEnrichmentplot_ui <- renderUI(simplifyEnrichmentplot_ui())



pws <- eventReactive(input$button3, input$selectedpw)

inter <- reactiveValues(pwss=NULL)


observeEvent({
  input$button3
}, {
  if (length(input$selectedpw) > 1) {
    inter$pwss <- pws()
  }else{
    show_alert(title="Warnings", 
               type="warning",
               text="Please select >=2 pathways to generate plots.")
  }
})






dlp <- function(ppp) {
  output[[paste0(ppp,"dl")]] <- downloadHandler(
    filename = function(x){paste0(ppp, ".pdf")},
    content = function(file) {
      pdf(file, width = 13, height = 10)
      p <- eval(parse(text = paste0(ppp, "1()")))
      print(p)
      dev.off()
    },
    contentType = "image/pdf"
  )
}
for (ppp in c("gseaHeatmap", "embeddedplot", "emapplotPie", "goplot", "emapplot")) {
  dlp(ppp)
}

output$clustercorplotdl <- downloadHandler(
  filename = function(x){
    ifelse(input$similarity == "Jaccard index", "clustercorplot_jaccard.pdf", "clustercorplot_pearson.pdf")
  },
  content = function(file) {
    pdf(file, width = 13, height = 10)
    if (input$similarity == "Jaccard index") {
      clustercorplot_jaccard(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, vertex.size.cex=input$vertex.size.cex,
                             vertex.label.cex=input$vertex.label.cex, edge.max.width=input$edge.max.width,
                             vertex.label.dist=input$vertex.label.dist)
    } else {
      clustercorplot(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, vertex.size.cex=input$vertex.size.cex,
                     vertex.label.cex=input$vertex.label.cex, edge.max.width=input$edge.max.width,
                     vertex.label.dist=input$vertex.label.dist)
    }
    dev.off()
  },
  contentType = "image/pdf"
)

output$circleplotdl <- downloadHandler(
  filename = "circleplot.pdf",
  content = function(file) {
    pdf(file, width = 13, height = 10)
    circleplot(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, pvaluecutoff = input$pvaluecutoff)
    dev.off()
  },
  contentType = "image/pdf"
)

output$hierarchyplot_treedl <- downloadHandler(
  filename = "hierarchyplot.pdf",
  content = function(file) {
    pdf(file, width = 13, height = 10)
    hierarchyplot_tree(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, topaths = input$topath1, 
                       vertex.size.cex = input$vertex.size.cex2,
                       vertex.label.cex = input$vertex.label.cex2,
                       edge.max.width = input$edge.max.width2, alpha.edge = input$alpha.edge)
    dev.off()
  },
  contentType = "image/pdf"
)


# output$simplifyEnrichmentplotdl <- downloadHandler(
#   filename = "simplifyEnrichmentplot.pdf",
#   content = function(file) {
#     pdf(file, width = 13, height = 10)
#     simplifyEnrichmentplot(SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, showCategory = input$showCategory4, GO_ont = input$ont2)
#     dev.off()
#   },
#   contentType = "image/pdf"
# )








# 
# 
# plotui <- fluidRow(
#   fluidRow(
#     column(width = 12,
#            boxPlus(
#              title = tagList(
#                icon("expand-arrows-alt"),
#                "Hierarchyplot"
#              ),
#              closable = F,
#              collapsible = T,
#              width = 12,
#              hierarchyplot_treeui(),
#              plotOutput('hierarchyplot_tree', height = "800px") %>% withSpinner()
#            )
#     )
#   ),
#   fluidRow(
#     column(
#       width = 6,
#       boxPlus(
#         title = tagList(
#           icon("th"),
#           "Heatmap"
#         ),
#         closable = F,
#         collapsible = T,
#         width = 12,
#         gseaheatmapui(),
#         plotOutput('gseaheatmap', height = "800px") %>% withSpinner()
#         # uiOutput('gseaheatmap1') %>% withSpinner()  # Warning较多
#       ),
#       boxPlus(
#         title = tagList(
#           icon("connectdevelop"),
#           "Clustercorplot"
#         ),
#         closable = F,
#         collapsible = T,
#         footer = tags$small(icon("lightbulb"), "Helps to infer the relationship and similarity between clusters."),
#         width = 12,
#         clustercorplotui(),
#         plotOutput('clustercorplot') %>% withSpinner()
#       ),
#       boxPlus(
#         title = tagList(
#           icon("bezier-curve"),
#           "Emapplot"
#         ),
#         closable = F,
#         collapsible = T,
#         width = 12,
#         emapplotui(),
#         plotOutput('emapplot') %>% withSpinner()
#       )#,
#       # uiOutput("simplifyEnrichmentplot_ui") %>% withSpinner()
#     ),
#     column(
#       width = 6,
#       boxPlus(
#         title = tagList(
#           icon("chart-bar"),
#           "Embeddedplot"
#         ),
#         closable = F,
#         collapsible = T,
#         width = 12,
#         embeddedplotui(),
#         plotOutput('embeddedplot') %>% withSpinner()
#       ),
#       boxPlus(
#         title = tagList(
#           icon("dharmachakra"),
#           # icon("vector-path-circle"),
#           "Circleplot"
#         ),
#         closable = F,
#         collapsible = T,
#         footer = tags$small(icon("lightbulb"), "Link width represents number of intersection pathways between clusters, 
#                             link color is the same as the cluster with higher-score intersection pathways."),
#         width = 12,
#         circleplotui(),
#         plotOutput('circleplot') %>% withSpinner()
#       ),
#       boxPlus(
#         title = tagList(
#           icon("bezier-curve"),
#           "EmapplotPie"
#         ),
#         closable = F,
#         collapsible = T,
#         width = 12,
#         emapplotPieui(),
#         plotOutput('emapplotPie') %>% withSpinner()
#       ),
#       uiOutput("goplot_ui") %>% withSpinner()
#     )
#   )
# )
# 
# 
# totalui <- renderUI({
#   # if (paste0("GSEAresult_", "GO") %in% names(SeuratObj@misc)) {
#   if (paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) {
#     fluidPage(
#       fluidRow(
#         boxPlus(
#           title = tagList(
#             icon("table"),
#             "GSEA result"
#           ),
#           closable = F,
#           collapsible = T,
#           width = 12,
#           DT::dataTableOutput("GSEAtable") %>% withSpinner(),
#           select_plotui
#         )
#       ),
#       plotui
#     )
#   } else {
#     tagList(tags$b("No GSEA result of "), dashboardLabel(input$switchTerm, status = "info", style = "square"), tags$b("were found."))
#   }
# })










































