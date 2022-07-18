source(
  file = "global.R",
  local = TRUE,
  encoding = "UTF-8"
)

server <- function(input, output, session) {
  
  source(
    file = "Utils/Developer.server.R",
    local = TRUE,
    encoding = "UTF-8"
  )
  source(
    file = "overview/overview.server.R",
    local = TRUE,
    encoding = "UTF-8"
  )
  source(
    file = "focus/focus.server.R",
    local = TRUE,
    encoding = "UTF-8"
  )
  source(
    file = "topicModel/topicModel.server.R",
    local = TRUE,
    encoding = "UTF-8"
  )
  

  

  
  
  
  
}





