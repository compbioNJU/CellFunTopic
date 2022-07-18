source(
  file = "global.R",
  local = TRUE,
  encoding = "UTF-8"
)


ui <- dashboardPagePlus(
  skin = "green",
  title = "cellFunMap",
  # options = list(sidebarExpandOnHover = TRUE),
  header = dashboardHeaderPlus(
    title = "cellFunMap",
    left_menu = tagList(
      tags$span(
        "Functional Annotation for Single-Cell Transcriptomics",
        style = 'font-family:"Helvetica Neue",Helvetica,Arial,sans-serif;color:#fff;line-height:34px;font-size:20px;font-weight:300;overflow:hidden;'
      )
    ),
    userOutput("user")
  ),
  sidebar = dashboardSidebar(
    # collapsed = TRUE,
    sidebarMenu(
      id = "sideBarTab",
      menuItem(
        "pathway overview",
        tabName = "overview",
        icon = icon("globe")
      ),
      menuItem(
        "pathway of focus",
        tabName = "focus",
        icon = icon("eye")
      ),
      menuItem(
        "Topic modeling",
        tabName = "topicModel",
        icon = icon("project-diagram")
        # icon = icon("cluster", lib="glyphicon")
      )
  )),
  body = dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    shinyjs::useShinyjs(),
    tabItems(
      tabItem(
        tabName = "overview",
        source(
          file = "overview/overview.ui.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value
      ),
      tabItem(
        tabName = "focus",
        source(
          file = "focus/focus.ui.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value
      ),
      tabItem(
        tabName = "topicModel",
        source(
          file = "topicModel/topicModel.ui.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value
      )
    )
  )
)
