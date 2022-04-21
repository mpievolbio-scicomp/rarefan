suppressMessages(library(shiny))
suppressMessages(library(shinycssloaders))

fluidPage(
		# App title
		titlePanel("REPIN and RAYT analysis"),
		sidebarLayout(
				sidebarPanel(
					textOutput("text"),
					tableOutput("run_setup"),
					selectInput(inputId = 'rayt',
								label="Select REPIN group",
								choices = list(
										"0" = 0,
										"1" = 1,
										"2" = 2,
										"3" = 3,
										"4" = 4,
										"5" = 5
								),
								selected = 0),
						htmlOutput("back_to_results")
						),
				mainPanel(
				  tabsetPanel(
				    tabPanel("RAYT tree",
    					withSpinner(plotOutput(outputId = 'rayt_tree', inline=TRUE))
				  ),
					tabPanel("REPINs",
  					withSpinner(plotOutput(outputId = 'repin_tree', inline=TRUE))
					),
					tabPanel("Correlation",
					         withSpinner(plotOutput(outputId = 'correlations', inline=TRUE))
					  )
				  )
				)
		)
)
