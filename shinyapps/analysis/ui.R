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
										"RAYT 0" = 0,
										"RAYT 1" = 1,
										"RAYT 2" = 2,
										"RAYT 3" = 3,
										"RAYT 4" = 4,
										"RAYT 5" = 5
								),
								selected = 0),
						textOutput("plot_instructions")
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
