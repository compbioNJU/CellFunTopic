fluidPage(
  fluidRow(
    id = "tm_hide1",
    br(),
    column(width = 12,
           boxPlus(
             width = 12,
             closable = F,
             # enable_label = T,
             # label_text = "New",
             # label_status = "warning",
             solidHeader = T,
             status = "danger",
             title = tagList(icon("bullhorn"), "warning"),
             collapsible = T,
             collapsed = F,
             tags$b("The number of topics is calculated automatically. To get ideal results, 
                    please run topic modeling with a specific ", tags$code("k"), "in advance. 
                    For example: \n"),
             tags$br(),
             tags$code("k <- 10"),
             tags$br(),
             tags$code('SeuratObj <- runLDA(SeuratObj, by = "GO", k = k, method = "VEM", SEED = 1234, plot = T)'),
             tags$br()
           )
    )
  ),
  fluidRow(
    id = "tm_hide2",
    fluidRow(
      column(
        width = 12,
        boxPlus(
          title = tagList(icon("project-diagram"), "Topic modeling network"),
          footer = tagList(icon("lightbulb"), "Click network to view wordcloud and barplot of specific topic."),
          closable = F,
          collapsible = T,
          collapsed = F,
          width = 12,
          fluidRow(
            column(width = 8,
                   numericInput(inputId = "tm_topn1", "show top terms:", 10, min = 1, max = 100, width = "300px"),
                   echarts4rOutput("tm_echt", height = "900px") %>% withSpinner()
                   ),
            column(width = 4,
                   plotOutput('tm_hist3') %>% withSpinner(),
                   plotOutput('tm_wordcloud', height = "500px") %>% withSpinner()
                   # wordcloud2Output('tm_wordcloud', width = "100%", height = "400px") %>% withSpinner()  # 会导致其他的图比如网络图加载不出来
            )
            )
        )
      )
      ),
    fluidRow(
      column(width = 12,
             boxPlus(
               title = tagList(icon("table"), "Topic modeling results"),
               closable = F,
               collapsible = T,
               collapsed = F,
               width = 12,
               tabsetPanel(
                 tabPanel(
                   title = "topic-term probabilities",
                   icon = icon("lightbulb"),
                   fluidRow(
                     column(width = 5,
                            boxPad(color = "gray",
                                   DT::dataTableOutput("betatable") %>% withSpinner()
                            )
                     ),
                     column(width = 7,
                            fluidRow(
                              column(width = 4, 
                                     uiOutput("topicsChoiceui")
                              ),
                              column(width = 4, 
                                     numericInput(inputId = "tm_topn2", "top terms:", 5, min = 1, max = 100)
                              ),
                              column(width = 4, 
                                     sliderInput("tm_size", label = "axis.text.y.size", value = 5,
                                                 min = 1, max = 30, step = 1)
                              )
                            ),
                            plotOutput('tm_hist1') %>% withSpinner()
                     )
                   )
                 ),
                 tabPanel(
                   title = "cluster-topic probabilities",
                   icon = icon("lightbulb"),
                   fluidRow(
                     column(width = 5,
                            boxPad(color = "gray",
                                   DT::dataTableOutput("gammatable") %>% withSpinner()
                            )
                     ),
                     column(width = 7,
                            fluidRow(
                              column(width = 12, 
                                     uiOutput("tm_clustersChoiceui")
                            ),
                            plotOutput('tm_hist2') %>% withSpinner()
                     )
                   )
                 )
               ),
               tabPanel(
                 title = "document-term-count-topic probabilities",
                 icon = icon("lightbulb"),
                 fluidRow(
                   column(width = 12,
                          DT::dataTableOutput("assignwstable") %>% withSpinner()
                   )
                 )
               )
             )
      )
    )
    ),
    fluidRow(
      column(width = 12,
             boxPlus(
               footer = fluidRow(
                 column(width = 4,
                        downloadButton("tm_sankeydl", label="Download html", style = "background-color:#1DA1F2;color:white;")
                 ),
                 column(width = 4,
                        downloadButton("tm_umap_clusterdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;")
                 ),
                 column(width = 4,
                        downloadButton("tm_heatmapdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;")
                 )
               ),
               closable = F,
               collapsible = T,
               collapsed = F,
               width = 12,
               fluidRow(
                 column(width = 4,
                        numericInput(inputId = "tm_topn3", "show top topics of every cluster:", 1, min = 1, max = 10, width = "300px"),
                        sankeyNetworkOutput("tm_sankey", height = "600px") %>% withSpinner()
                        ),
                 column(width = 4,
                        plotOutput('tm_umap_cluster') %>% withSpinner()
                 ),
                 column(width = 4,
                        plotOutput('tm_heatmap', height = "600px") %>% withSpinner()
                 )
               )
               )
             )),
    fluidRow(
      column(width = 5,
             boxPlus(
               footer = tagList(
                 icon("lightbulb"), "Show cosine similarity between topics.",
                 downloadButton("tm_cosineheatmapdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;")
               ),
               closable = F,
               collapsible = T,
               collapsed = F,
               width = 12,
               plotOutput('tm_cosineheatmap') %>% withSpinner()
             )
             ),
      column(width = 7,
             boxPlus(
               footer = tagList(icon("lightbulb"), "Show cosine similarity between clusters."),
               closable = F,
               collapsible = T,
               collapsed = F,
               width = 12,
               fluidRow(
                 column(width = 4,
                        boxPad(color = "gray",
                               radioGroupButtons(
                                 inputId = "tm_layout",
                                 label = "layout", 
                                 choices = c("force-directed"="fr", "circle"="circle"),
                                 status = "info",
                                 individual = T
                               ),
                               sliderInput("cos_sim_thresh", label = "cosine similarity threshold", value = 0.2,
                                           min = 0.01, max = 1, step = 0.01),
                               sliderInput("radiuscosnet", label = "pie radius", value = 0.12,
                                           min = 0.01, max = 1, step = 0.01),
                               sliderTextInput(
                                 inputId = "width_range",
                                 label = "edge width range:", 
                                 choices = seq(from = 0, to = 2, by = 0.1),
                                 selected = c(0.1, 0.8),
                                 grid = TRUE
                               ),
                               downloadButton("tm_cosine_network_clusterdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;")
                        )
                 ),
                 column(width = 8,
                        plotOutput('tm_cosine_network_cluster') %>% withSpinner()
                        # flipBox(
                        #   id = 1,
                        #   width = 12,
                        #   header_img = NULL,
                        #   main_img = "project-diagram-solid.svg",
                        #   front_btn_text = "circle layout",
                        #   back_btn_text = "force-directed layout",
                        #   plotOutput('tm_cosine_network_cluster_fr') %>% withSpinner(),
                        #   back_content = plotOutput('tm_cosine_network_cluster_circle') %>% withSpinner()
                        # )
                 )
               )
             )
      )
    ),
    fluidRow(
      column(width = 12,
             boxPlus(
               footer = tagList(icon("lightbulb"), "Show cosine similarity between terms."),
               closable = F,
               collapsible = T,
               collapsed = F,
               width = 12,
               fluidRow(
                 column(width = 3,
                        boxPad(color = "gray",
                               radioGroupButtons(
                                 inputId = "cosine_cal_by",
                                 label = "cosine similarity between terms calculated by", 
                                 choices = c("Topic modeling", "GSEA result"),
                                 status = "primary",
                                 individual = T
                               ),
                               sliderInput("tm_cosnet_topn", label = "show top terms of every topic/cluster", value = 10,
                                           min = 1, max = 100, step = 1),
                               selectInput(
                                 inputId = "tm_layout2",
                                 label = "layout", 
                                 selected = "fr",
                                 choices = c("fr", "kk", "star", "circle", "dh", "graphopt", "graphopt", "sphere", "drl", "nicely", "lgl")
                               ),
                               sliderInput("cos_sim_thresh2", label = "cosine similarity threshold", value = 0.8,
                                           min = 0.01, max = 1, step = 0.01),
                               sliderInput("radiuscosnet2", label = "pie radius", value = 0.2,
                                           min = 0.01, max = 1, step = 0.01),
                               sliderTextInput(
                                 inputId = "width_range2",
                                 label = "edge width range:", 
                                 choices = seq(from = 0, to = 2, by = 0.05),
                                 selected = c(0.05, 0.55),
                                 grid = TRUE
                               ),
                               sliderInput("tm_text_size", label = "text size", value = 4,
                                           min = 0.5, max = 20, step = 0.5),
                               downloadButton("tm_cosine_network_termdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;")
                        )
                 ),
                 column(width = 9,
                        plotOutput('tm_cosine_network_term', height = "800px") %>% withSpinner()
                 )
               )
             )
             )
    ),
    fluidRow(
      id = "tm_hide3",
      column(width = 12,
             boxPlus(
               title = tagList(icon("chart-line"), "Metrics for the best number of topics"),
               footer = tagList(
                 icon("lightbulb"), "Helps to choose the best number of topics.",
                 downloadButton("tm_indexdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;")
               ),
               # footer = tagList(icon("lightbulb"), "Helps to choose the best number of topics."),
               closable = F,
               collapsible = T,
               collapsed = F,
               width = 12,
               plotOutput('tm_index') %>% withSpinner()
             )
      ))
  )
)