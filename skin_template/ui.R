library(shinycssloaders)


fluidPage(


  HTML('<meta name="viewport" content="width=1024">'),


  theme = shinytheme("united"),

  fluidRow(style = "padding-top:5px;",
    column(4,
        style='padding:1x;',
        offset = 0,
        align="center",
        #imageOutput("Logo", height="10vh")
        tags$a(img(src='skingenes_logo.png', style="width: 66%"), href="../index_inner.html")    
    ),
    column(4,
        style='padding:0px;',
           offset = 0,

	align="center",

	wellPanel(align="left",
	      
        fixedRow(style='padding-left: 5px; padding-right: 5px; padding-top: 0px; padding-bottom: 0px;',
		column(10,
		        style='padding-left: 5px; padding-right: 0px; padding-top: 0px; padding-bottom: 0px;',
           		offset = 0,

			align="left",
		
#          		selectizeInput('geneNames', "Select up to 3 genes",
#               			choices=NULL, selected=sample(top100,3), multiple=TRUE, options = list(placeholder ='Start typing gene name', maxItems = 3, plugins=list("remove_button")))
          		selectizeInput('geneNames', "Select genes to explore the dataset",
               			choices=NULL, multiple=TRUE, selected=NULL, options = list(placeholder ='Type up to 3 genes and click Go! (leave blank to sample from top 100 genes)', maxItems = 3, plugins=list("remove_button")))
		),
		column(2,
		        style='margin-top: 22px; padding-right:5px; padding-left: 0px; padding-top: 0px; padding-bottom: 0px;',
           		offset = 0,

			align="center",

			actionButton(
        			inputId = "submit_loc",
        			label = "Go",
				width = "90%"
			)
		)
	),

	style="padding-left: 2%; padding-right: 2%; padding-top: 0px; padding-bottom: 0px; width: 80%, height: 100%"
        )

    ),
    column(4,
        style='padding:0px;',
           offset = 0,

        align="center"
    )
  ),

  fluidRow(style='padding:1px;',
    #column(4,
    #       offset = 0,
    #       align="left",
    #       plotOutput("Main") %>% withSpinner(color="#66cc66")
    #),
    column(4,
           style='padding-left: 5px;',
           offset = 0,
           align="left",
           includeHTML("description.html")
    ),
#    column(4, style='padding:0px;', offset = 0, align="center",
#           #img(src='tsne.svg',style="width: 100%")
#           plotOutput("Main")
#    ),
    column(8,
           style='padding:0px;',
           offset = 0,
           align="center",
           #plotOutput("Draw", height="46vh")
           img(src='Mouse_wound_day12.jpg',style="width: 90%")
    )
  ),
  fluidRow(style='padding:0px',

           tags$hr(),

    column(12,
           style='padding:0px:',
           align="center",
           
#           uiOutput("feature1")
           uiOutput("feature1")
    )
  ),

  fluidRow(style='padding:0px',
    column(12,
           
           align="center",
           
           tags$hr(),
           tabsetPanel(id = "miscTab",
             tabPanel("Violins",plotOutput("vlnPlot", height="auto") %>% withSpinner(color="#cc66cc")),
             tabPanel("Ridges",plotOutput("ridgePlot", height="auto") %>% withSpinner(color="#21b700")),
             tabPanel("HeatMap",plotOutput("heatPlot") %>% withSpinner(color="#ff65ab"))
           )
    )

#  ),

  
  #fluidRow(
  #   column(1,

  #             dropdownButton(
  #    tags$h3("Plot Options"),
  #         wellPanel(
  #           radioButtons(inputId = "pt1Type", label="Reduction for Plot 1", choices = reductionChoices, selected = reductionChoices[[2]]),
  #           radioButtons(inputId = "labelBoolean", label="Display labels?", choices = c("TRUE","FALSE"), selected="FALSE"),
  #           tags$hr(),
  #           sliderInput(inputId = "dotSize", label = "Set point size", value=0.1, min=0.01, max=10),
  #           sliderInput(inputId = "labelSize", label = "Set label size", value=6, min=0.5, max=10)
  #         ),
  #    circle = FALSE, status = "info", icon = icon("gear"), width = "300px",
  #    circle = FALSE, status = "default", icon = NULL, width = "0px",
  #    tooltip = tooltipOptions(title = "Click to see Options!"), inputId = "dropBottom"
  #  )

  # )
  ,
  tags$script('
        $(document).on("keydown", function (e) {
        if (e.keyCode == "13") {
  		$("#submit_loc").click();
  	}
        });



  ')

  )
)
