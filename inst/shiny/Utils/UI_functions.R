circleplotui <- function() {
  fluidRow(
    column(width = 10),
    column(width = 2,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               sliderInput("pvaluecutoff", label = "pvalue cutoff", value = 0.01,
                           min = 0, max = 0.05, step = 0.001),
               sliderInput("link_threshold_circleplot", label = "threshold number of intersection", value = 10,
                           min = 1, max = 100, step = 2),
               downloadButton("circleplotdl", label="Download PDF", status = "success", style = "background-color:#1DA1F2;color:white;"),
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
    )
  )
}

clustercorplotui <- function() {
  fluidRow(
    column(width = 10),
    column(width = 2,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               radioGroupButtons(
                 inputId = "similarity",
                 label = "Similarity measurement",
                 choices = c("Jaccard index", "Pearson's correlation"),
                 status = "danger",
                 checkIcon = list(
                   yes = icon("check-square"),
                   no = icon("square-o")
                 ),
                 selected = "Jaccard index",
                 individual = F,
                 direction = "vertical"
               ),
               sliderInput("link_threshold_clustercorplot", label = "only show links whose correlation/Jaccard-index bigger than", value = 0.5,
                           min = 0, max = 1, step = 0.05),
               sliderInput("vertex.size.cex", label = "vertex.size.cex", value = 1,
                           min = 0, max = 5, step = 0.1),
               sliderInput("vertex.label.cex", label = "vertex.label.cex", value = 1.5,
                           min = 0, max = 5, step = 0.1),
               sliderInput("edge.max.width", label = "edge.max.width", value = 8,
                           min = 1, max = 30, step = 1),
               sliderInput("vertex.label.dist", label = "vertex.label.dist", value = 2,
                           min = 0, max = 10, step = 0.5),
               downloadButton("clustercorplotdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;"),
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
    )
  )
}


gseaheatmapui <- function() {
  fluidRow(
    column(width = 7),
    column(width = 4,
           pickerInput(
             inputId = "hmpcol",
             label = "color",
             choices = palettes,
             choicesOpt = color_choicesOpt(palettes, id=3),
             inline = T
           )
    ),
    column(width = 1,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               sliderInput("topPath", label = "top pathways of every cluster to show", value = 10,
                           min = 1, max = 30, step = 1),
               sliderInput("fontsize_row", label = "row fontsize", value = 10,
                           min = 1, max = 30, step = 1),
               selectInput("toshow", label = "heatmap showing",
                           choices = list("-log10(p.adjust)" = "-logFDR", "normalized enrichment score(NES)" = "NES"),
                           selected = "-logFDR"),
               radioGroupButtons(
                 inputId = "scale",
                 label = "scale",
                 choices = c("none", "row", "column"),
                 status = "primary",
                 individual = T
               ),
               downloadButton("gseaHeatmapdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;"),
               # sliderInput("plotheight", label = "height of plot", value = 800,
               #             min = 100, max = 2000, step = 50),
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
    )
  )
}


embeddedplotui <- function() {
  fluidRow(
    column(width = 10),
    column(width = 2,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               sliderInput("topaths", label = "top pathways of every cluster to show", value = 1,
                           min = 1, max = 5, step = 1),
               sliderInput("pie.size.cex", label = "pie.size.cex", value = 1,
                           min = 0, max = 5, step = 0.1),
               radioGroupButtons(
                 inputId = "reduction",
                 label = "reduction",
                 choices = c("umap", "tsne", "pca"),
                 status = "primary",
                 individual = T
               ),
               radioGroupButtons(
                 inputId = "type",
                 label = "child plot",
                 choiceNames = c(
                   paste(icon("chart-pie"), "pie"),
                   paste(icon("chart-bar"), "hist")
                 ),
                 choiceValues = c("pie", "hist"),
                 status = "warning",
                 individual = T
               ),
               downloadButton("embeddedplotdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;"),
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
    )
  )
}



hierarchyplot_treeui <- function() {
  fluidRow(
    column(width = 10),
    column(width = 2,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               sliderInput("topath1", label = "top pathways of every cluster to show", value = 5,
                           min = 1, max = 50, step = 1),
               # selectInput("cluster_cutree_k", label = "cluster_cutree_k",  # NULL无法成为其中的选项，选项默认为1
               #             choices = append(list('NULL'=NULL), 1:100),
               #             selected = 'NULL'),
               # selectInput("pathway_cutree_k", label = "pathway_cutree_k",
               #             choices = append(list('NULL'=NULL), 1:100),
               #             selected = 'NULL'),
               sliderInput("edge.max.width2", label = "edge.max.width", value = 6,
                           min = 0, max = 20, step = 0.5),
               sliderInput("vertex.size.cex2", label = "vertex.size.cex", value = 1.2,
                           min = 0, max = 10, step = 0.1),
               sliderInput("vertex.label.cex2", label = "vertex.label.cex", value = 1.2,
                           min = 0, max = 10, step = 0.1),
               sliderInput("alpha.edge", label = "transparency of edge", value = 0.6,
                           min = 0, max = 1, step = 0.05),
               downloadButton("hierarchyplot_treedl", label="Download PDF", style = "background-color:#1DA1F2;color:white;"),
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
    )
  )
}

emapplotPieui <- function() {
  fluidRow(
    column(width = 10),
    column(width = 2,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               sliderInput("showCategory", label = "top pathways of every cluster to show", value = 5,
                           min = 1, max = 50, step = 1),
               sliderInput("cex_line", label = "width of link", value = 1,
                           min = 0, max = 5, step = 0.1),
               sliderInput("node_size_cex", label = "node.size.cex", value = 1.5,
                           min = 0, max = 10, step = 0.1),
               sliderInput("node_label_cex", label = "node.label.cex", value = 1,
                           min = 0, max = 10, step = 0.1),
               radioGroupButtons(
                 inputId = "pie",
                 label = "draw pie by",
                 choices = c("-log10FDR", "count"),
                 status = "primary",
                 individual = T
               ),
               selectInput(
                 inputId = "layout",
                 label = "layout",
                 choices = c("nicely", "kk", "star", "circle", "dh", "graphopt", "graphopt", "sphere", "drl", "fr", "lgl")
               ),
               downloadButton("emapplotPiedl", label="Download PDF", style = "background-color:#1DA1F2;color:white;"),
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
    )
  )
}



emapplotui <- function() {
  fluidRow(
    column(width = 10),
    column(width = 2,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               uiOutput("clusterChoiceui"),
               sliderInput("showCategory2", label = "top pathways of every cluster to show", value = 20,
                           min = 1, max = 50, step = 1),
               sliderInput("cex_line2", label = "width of link", value = 1,
                           min = 0, max = 5, step = 0.1),
               sliderInput("node_size_cex2", label = "node.size.cex", value = 1.5,
                           min = 0, max = 10, step = 0.1),
               sliderInput("node_label_cex2", label = "node.label.cex", value = 1,
                           min = 0, max = 10, step = 0.1),
               selectInput(
                 inputId = "layout2",
                 label = "layout",
                 choices = c("nicely", "kk", "star", "circle", "dh", "graphopt", "graphopt", "sphere", "drl", "fr", "lgl")
               ),
               downloadButton("emapplotdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;"),
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
    )
  )
}


goplotui <- function() {
  fluidRow(
    column(width = 10),
    column(width = 2,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               uiOutput("clusterChoiceui2"),
               sliderInput("showCategory3", label = "top pathways of every cluster to show", value = 10,
                           min = 1, max = 50, step = 1),
               sliderInput("label_size", label = "label size", value = 3,
                           min = 0.1, max = 10, step = 0.5),
               radioGroupButtons(
                 inputId = "ont",
                 label = "GO ontology",
                 choices = c("BP", "MF", "CC"),
                 status = "primary",
                 individual = T
               ),
               downloadButton("goplotdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;"),
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
    )
  )
}



pathwayScatterplotui <- function() {
  fluidRow(
    column(width = 7),
    column(width = 4,
           pickerInput(
             inputId = "sctcol",
             label = "color",
             choices = palettes,
             selected = "OrRd",
             choicesOpt = color_choicesOpt(palettes, id=3),
             inline = T
           )
    ),
    column(width = 1,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               materialSwitch(
                 inputId = "labelornot",
                 label = tags$b("Label the clusters"),
                 status = "primary",
                 value = TRUE
               ),
               sliderInput("pointsize", label = "point size", value = 1,
                           min = 0, max = 10, step = 0.1),
               sliderInput("label.size", label = "label size", value = 4,
                           min = 1, max = 10, step = 0.2),
               radioGroupButtons(
                 inputId = "reduction1",
                 label = "reduction",
                 choices = c("umap", "tsne", "pca"),
                 status = "primary",
                 individual = T
               ),
               downloadButton("pathwayScatterplotdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;"),
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
    )
  )
}




GOboxplotui <- function() {
  fluidRow(
    column(width = 11),
    column(width = 1,
           tags$span(
             dropdownButton(
               tags$h4(icon("sliders-h"), "Options"),
               tags$hr(),
               materialSwitch(
                 inputId = "flip",
                 label = tagList(icon("sync-alt"), tags$b("Flip the coordinate")),
                 status = "primary",
                 value = FALSE
               ),
               sliderInput("pointsize2", label = "point size", value = 1,
                           min = 0, max = 10, step = 0.1),
               downloadButton("GOboxplotdl", label="Download PDF", style = "background-color:#1DA1F2;color:white;"),
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
    )
  )
}


# simplifyEnrichmentplotui <- function() {
#   fluidRow(
#     column(width = 11),
#     column(width = 1,
#            tags$span(
#              dropdownButton(
#                tags$h4(icon("sliders-h"), "Options"),
#                tags$hr(),
#                sliderInput("showCategory4", label = "top pathways of every cluster to show", value = 10,
#                            min = 1, max = 100, step = 1),
#                radioGroupButtons(
#                  inputId = "ont2",
#                  label = "GO ontology",
#                  choices = c("BP", "MF", "CC"),
#                  status = "primary",
#                  individual = T
#                ),
#                downloadButton("simplifyEnrichmentplotdl", label="Download PDF"),
#                circle = F,
#                right = T,
#                inline = T,
#                status = "danger",
#                icon = icon("gear"),
#                size = "sm",
#                width = "300px",
#                tooltip = tooltipOptions(title = "Options", placement = "top")
#              ),
#              style = "float:right;"
#            )
#     )
#   )
# }






