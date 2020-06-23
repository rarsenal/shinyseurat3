fluidPage(theme = shinytheme("superhero"),

  titlePanel(titlename),
  # First row containing: input field for user gene selection
  fluidRow(
    wellPanel(
      selectizeInput(
      'geneNames', 'Select Genes',
      choices=NULL, selected=sample(top100,3), multiple=TRUE, options = list(placeholder ='Start typing gene name'), width = '3000px'
    ))
  ), # Second row containing the configuration button and its popup menu panel, and then the main reduction map
  fluidRow(style='height:50vh',
    dropdownButton(
      tags$h3("Plot Options"),
           wellPanel(
             radioButtons(inputId = "pt1Type", label="Reduction for Plot 1", choices = reductionChoices, selected = reductionChoices[[2]]),
             radioButtons(inputId = "pt2Type", label="Reduction for Plot 2", choices = reductionChoices, selected = reductionChoices[[3]]),
             radioButtons(inputId = "labelBoolean", label="Display labels?", choices = c("TRUE","FALSE"), selected="TRUE"),
             radioButtons(inputId = "legendBoolean", label="Legend Position", choices = c("none","right","bottom"), selected="right"),
             tags$hr(),
             sliderInput(inputId = "dotSize", label = "Set point size", value=0.1, min=0.01, max=10),
             sliderInput(inputId = "labelSize", label = "Set label size", value=6, min=0.5, max=10)
           ),shinythemes::themeSelector(),
      circle = FALSE, status = "info", icon = icon("gear"), width = "300px",
      tooltip = tooltipOptions(title = "Click to see Options!")
    ),
    column(12,
           
           align="center",
           plotOutput("Main", width="80%", height="50vh")
           
    )
  ), # Thrid row containing the tabbed panel with the various plots
  fluidRow(
    column(12,
           
           align="center",
           
           tags$hr(),
           tabsetPanel(
             tabPanel("Feature 1",plotOutput("feature1")),
             tabPanel("Feature 2",plotOutput("feature2")),
             tabPanel("Violins",plotOutput("vlnPlot")),
             tabPanel("Ridges",plotOutput("ridgePlot")),
             tabPanel("Dot Plots",plotOutput("dotPlot")),
             tabPanel("HeatMap",plotOutput("heatPlot"))
           )
    )
  )
)
