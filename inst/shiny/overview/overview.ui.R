fluidPage(
  fluidRow(
    boxPlus(
      title = tagList(
        icon("table"),
        "GSEA result"
      ),
      closable = F,
      collapsible = T,
      width = 12,
      fluidRow(
        column(
          width = 11,
          radioGroupButtons(
            justified = F,
            inputId = "switchTerm",
            label = NULL,
            choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN", "customized"),
            # selected = "GO",
            checkIcon = list(
              yes = icon("ok", lib = "glyphicon")
            ),
            status = "danger",
            individual = T
          )
        ),
        column(
          width = 1,
          uiOutput("grouping")
        )
      ),
      DT::dataTableOutput("GSEAtable") %>% withSpinner(),
      # uiOutput("select_plotui") %>% withSpinner()
      fluidRow(
        id = "tohide2",
        column(width = 12,
               uiOutput("select_plotui2") %>% withSpinner()
               )
      ),
      fluidRow(
        id = "tohide1",
        br(),
        column(width = 3,
               actionBttn(
                 inputId = "button1",
                 label = "Select All",
                 style = "fill",
                 color = "success",
                 icon=icon("check square"),
                 size = "sm"
               ),
               actionBttn(
                 inputId = "button2",
                 label = "Clear Selection",
                 style = "simple",
                 color = "primary",
                 icon=icon("broom"),
                 size = "sm"
               )
        ),
        column(width = 5,
               uiOutput("pathwayList")
        ),
        column(width = 4,
               tags$b("Send selected pathway(s) to:"),
               actionBttn(
                 inputId = "button3",
                 label = "Overview",
                 style = "jelly",
                 icon=icon("arrow-circle-right"),
                 color = "danger",
                 size = "sm"
               ),
               actionBttn(
                 inputId = "button4",
                 label = "Focus",
                 icon=icon("arrow-circle-right"),
                 style = "unite",
                 color = "warning",
                 size = "sm"
               ),
        )
      )
    )
  ),
  # fluidRow(
  #   id = "pa",
  #   boxPlus(
  #     title = tagList(
  #       icon("connectdevelop"),
  #       "Pathway Analysis"
  #     ),
  #     closable = F,
  #     collapsible = T,
  #     width = 12,
  #     uiOutput("plotui") %>% withSpinner()
  #   )
  #   
  # ),
  fluidRow(
    id = "pa",
    fluidRow(
      column(width = 12,
             boxPlus(
               title = tagList(
                 icon("expand-arrows-alt"),
                 "Hierarchyplot"
               ),
               closable = F,
               collapsible = T,
               width = 12,
               hierarchyplot_treeui(),
               plotOutput('hierarchyplot_tree', height = "800px") %>% withSpinner()
             )
      )
    ),
    fluidRow(
      column(
        width = 6,
        boxPlus(
          title = tagList(
            icon("th"),
            "Heatmap"
          ),
          closable = F,
          collapsible = T,
          width = 12,
          gseaheatmapui(),
          plotOutput('gseaheatmap', height = "980px") %>% withSpinner()
          # uiOutput('gseaheatmap1') %>% withSpinner()  # Warning较多
        )
      ),
      column(
        width = 6,
        boxPlus(
          title = tagList(
            icon("chart-bar"),
            "Embeddedplot"
          ),
          closable = F,
          collapsible = T,
          width = 12,
          embeddedplotui(),
          plotOutput('embeddedplot', height = "400px") %>% withSpinner()
        ),
        boxPlus(
          title = tagList(
            icon("dharmachakra"),
            # icon("vector-path-circle"),
            "Circleplot"
          ),
          closable = F,
          collapsible = T,
          footer = tags$small(icon("lightbulb"), "Link width represents number of intersection pathways between clusters,
                            link color is the same as the cluster with higher-score intersection pathways."),
          width = 12,
          circleplotui(),
          plotOutput('circleplot', height = "400px") %>% withSpinner()
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        boxPlus(
          title = tagList(
            icon("connectdevelop"),
            "Clustercorplot"
          ),
          closable = F,
          collapsible = T,
          footer = tags$small(icon("lightbulb"), "Helps to infer the relationship and similarity between clusters."),
          width = 12,
          clustercorplotui(),
          plotOutput('clustercorplot', height = "400px") %>% withSpinner()
        )
      ),
      column(
        width = 6,
        boxPlus(
          title = tagList(
            icon("bezier-curve"),
            "EmapplotPie"
          ),
          closable = F,
          collapsible = T,
          width = 12,
          emapplotPieui(),
          plotOutput('emapplotPie', height = "400px") %>% withSpinner()
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        boxPlus(
          title = tagList(
            icon("bezier-curve"),
            "Emapplot"
          ),
          closable = F,
          collapsible = T,
          width = 12,
          emapplotui(),
          plotOutput('emapplot', height = "400px") %>% withSpinner()
        )
      ),
      column(
        width = 6,
        # uiOutput("goplot_ui") %>% withSpinner()
        boxPlus(
          id = "tohide3",
          title = tagList(
            icon("sitemap"),
            "GOplot"
          ),
          closable = F,
          collapsible = T,
          width = 12,
          goplotui(),
          plotOutput('goplot', height = "400px") %>% withSpinner()
        )
      )
    )#,
    # fluidRow(
    #   column(
    #     width = 6,
    #     uiOutput("simplifyEnrichmentplot_ui") %>% withSpinner()
    #   ),
    #   column(
    #     width = 6
    #   )
    # )
  )
  
)

# fluidPage(
#   fluidRow(
#     column(
#       width = 11,
#       radioGroupButtons(
#         justified = F,
#         inputId = "switchTerm",
#         label = NULL,
#         choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN", "customized"),
#         # selected = "GO",
#         checkIcon = list(
#           yes = icon("ok", lib = "glyphicon")
#         ),
#         status = "danger",
#         individual = T
#       )
#     ),
#     column(
#       width = 1#,
#       # uiOutput("grouping")
#     )
#   ),
#   uiOutput("totalui") %>% withSpinner()
# )





