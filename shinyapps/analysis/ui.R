suppressMessages(library(shiny))
suppressMessages(library(shinycssloaders))

fluidPage(
		# App title
		titlePanel("REPIN and RAYT analysis"),
		sidebarLayout(
				sidebarPanel(
					textOutput("text"),
					selectInput(inputId = 'rayt',
								label="Select REP/RAYT group",
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
					h4("RAYT tree"),
					withSpinner(plotOutput(outputId = 'rayt_tree', width="100%")),
					hr(),
					h4("REPINs"),
					withSpinner(plotOutput(outputId = 'repin_tree', width="100%")),
					hr(),
					h4("Correlation"),
				  withSpinner(plotOutput(outputId = 'correlations', width="100%"))
				)
		)
)
