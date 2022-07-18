output$user <- renderUser({
  dashboardUser(
    name = "Shanni Cao",
    src = "profile.png",
    title = "Developer",
    subtitle = "https://github.com/xiaobaicainiao666",
    fluidRow(
      tags$div(
        style = "text-align:center",
        HTML('<p>GitHub: <a href="https://github.com/xiaobaicainiao666">https://github.com/xiaobaicainiao666</a></p>'),
      )
    )
  )
})


