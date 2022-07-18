

pathwayScatterplot1 <- reactive({
  pathwayScatterplot(SeuratObj, by = input$switchTerm, pathwayID = interr$pwidd,
                     reduction = input$reduction1, colour = input$sctcol, pointsize = input$pointsize, 
                     label = input$labelornot, label.size = input$label.size)
})

output$pathwayScatterplot <- renderPlot({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if (nn %in% names(SeuratObj@misc)) {
    pathwayScatterplot1()
  } else {
    return(NULL)
  }
})


GOboxplot1 <- reactive({
  GOboxplot(SeuratObj, goid = interr$pwidd, type = input$childparent, pointsize = input$pointsize2, flip = input$flip)
})

output$GOboxplot <- renderPlot({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if ((nn %in% names(SeuratObj@misc)) & (input$switchTerm == "GO")) {
    GOboxplot1()
  } else {
    return(NULL)
  }
})

# 
# gseaHeatmap2 <- reactive({
#   if (input$childparent == "child") {
#     nodes <- get_child_nodes(interr$pwidd)$child_go_id
#   } else if (input$childparent == "parent") {
#     nodes <- get_parent_nodes(interr$pwidd)$parent_go_id
#   }
#   gseaHeatmap(SeuratObj, by = "GO", pathwayIDs = nodes, toshow = input$toshow2,
#               colour = input$hmpcol2, scale = input$scale2, fontsize_row = input$fontsize_row2)
# })
# # output$gseaheatmap <- renderPlot(gseaHeatmap1())
# output$gseaheatmap2 <- renderPlot({
#   nn <- paste0("GSEAresult_", input$switchTerm)
#   if ((nn %in% names(SeuratObj@misc)) & (input$switchTerm == "GO")) {
#     gseaHeatmap2()
#   } else {
#     return(NULL)
#   }
# })
# 
# 
# gseaheatmap2ui <- function() {
#   fluidRow(
#     column(width = 7),
#     column(width = 4,
#            pickerInput(
#              inputId = "hmpcol2",
#              label = "color", 
#              choices = palettes,
#              choicesOpt = color_choicesOpt(palettes, id=3),
#              inline = T
#            )
#     ),
#     column(width = 1,
#            tags$span(
#              dropdownButton(
#                tags$h4(icon("sliders-h"), "Options"),
#                tags$hr(),
#                sliderInput("fontsize_row2", label = "row fontsize", value = 10,
#                            min = 1, max = 30, step = 1),
#                selectInput("toshow2", label = "heatmap showing",
#                            choices = list("-log10(p.adjust)" = "-logFDR", "normalized enrichment score(NES)" = "NES"),
#                            selected = "-logFDR"),
#                radioGroupButtons(
#                  inputId = "scale2",
#                  label = "scale", 
#                  choices = c("none", "row", "column"),
#                  status = "primary",
#                  individual = T
#                ),
#                downloadButton("gseaHeatmap2dl", label="Download PDF"),
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

termtable <- reactive({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if ((nn %in% names(SeuratObj@misc)) & (input$switchTerm == "GO")) {
    if (input$childparent == "child") {
      df <- get_child_nodes(interr$pwidd)
    } else if (input$childparent == "parent") {
      df <- get_parent_nodes(interr$pwidd)
    }
  } else {
    return(NULL)
  }
  
  DT::datatable(
    data = df,
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
          targets = c(0,1,3)
        ),
        list(
          width = "270px",
          className = "dt-left",
          targets = 2
        ),
        list(width = "50px", targets = c(3))
      )
    )
  ) %>%
    formatStyle(
      columns = "distance",
      backgroundColor = styleEqual(
        sort(unique(df$distance)),
        colorRampPalette(c("#FFFFFF", "#B03C2D"))(length(unique(df$distance)))
      ),
      fontWeight = "bold"
    )
})

output$termtable <- DT::renderDataTable(
  termtable()
)


output$termsui <- renderUI({
  fluidRow(
    column(width = 7,
           boxPlus(
             title = tagList(tags$b(paste0(input$childparent, " terms"))),
             closable = F,
             collapsible = T,
             status = "warning",
             solidHeader = T,
             width = 12,
             DT::dataTableOutput("termtable") %>% withSpinner()
           )
    ),
    column(width = 5,
           boxPlus(
             title = tagList(tags$b(paste0(input$childparent, " terms"))),
             # title = dashboardLabel(paste0(input$childparent, " terms"), status = "info", style = "square"),
             closable = F,
             status = "warning", 
             solidHeader = T,
             collapsible = T,
             width = 12,
             GOboxplotui(),
             plotOutput('GOboxplot')
           )),
  )
})

output$focusui_sub <- renderUI({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if ((nn %in% names(SeuratObj@misc)) & (input$switchTerm == "GO")) {
    fluidRow(
      column(width = 12,
             boxPlus(
               closable = F,
               collapsible = T,
               width = 12,
               fluidRow(
                 column(width = 8,
                        radioGroupButtons(
                          justified = F,
                          inputId = "childparent",
                          label = NULL,
                          choices = c("Child Terms" = "child", "Parent Terms" = "parent"),
                          checkIcon = list(
                            yes = icon("ok", lib = "glyphicon")
                          ),
                          status = "danger",
                          individual = T
                        )
                 )
               ),
               uiOutput("termsui") %>% withSpinner()
             )
             )
    )
  } else {
    return(NULL)
  }
})



pwid <- eventReactive(input$button4, input$selectedpw)

interr <- reactiveValues(pwidd=NULL)


output$focusui <- renderUI({
  if ((paste0("GSEAresult_", input$switchTerm) %in% names(SeuratObj@misc)) & (length(pwid()) == 1)) {
    fluidPage(
      fluidRow(
        column(width = 6,
               boxPlus(
                 closable = F,
                 collapsible = T,
                 width = 12,
                 pathwayScatterplotui(),
                 plotOutput('pathwayScatterplot') %>% withSpinner()
               ))
      ),
      uiOutput("focusui_sub") %>% withSpinner()
    )
  } else {
    return(NULL)
  }
})


observeEvent({
  input$button4
}, {
  if (length(input$selectedpw) == 1) {
    interr$pwidd <- pwid()
    
    updateTabItems(
      session = session,
      inputId = "sideBarTab",
      selected = "focus"
    )
  }else{
    show_alert(title="Warnings", 
               type="warning",
               text="Please select only one pathway.")
  }
})



output$pathwayScatterplotdl <- downloadHandler(
  filename = "pathwayScatterplot.pdf",
  content = function(file) {
    pdf(file, width = 13, height = 10)
    p <- pathwayScatterplot1()
    print(p)
    dev.off()
  },
  contentType = "image/pdf"
)

output$GOboxplotdl <- downloadHandler(
  filename = "GOboxplot.pdf",
  content = function(file) {
    pdf(file, width = 13, height = 10)
    p <- GOboxplot1()
    print(p)
    dev.off()
  },
  contentType = "image/pdf"
)











