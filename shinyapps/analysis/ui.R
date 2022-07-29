suppressMessages(library(shiny))
suppressMessages(library(shinycssloaders))

fluidPage(
		# App title
		titlePanel("REPIN and RAYT analysis"),
		sidebarLayout(
				sidebarPanel(
					textOutput("text"),
					tableOutput("run_setup"),
					numericInput(inputId = 'rayt',
            								label="Select REPIN group",
            								value=0,
            								min=0,
            								step=1
					),
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
