suppressMessages(library(shiny))
# suppressMessages(library(plotly))

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
					plotOutput(outputId = 'rayt_tree', width="100%"),
					hr(),
					h4("REPINs"),
					plotOutput(outputId = 'repin_tree', width="100%"),
					hr(),
					h4("Correlation"),
					plotOutput(outputId = 'correlations', width="100%")
				)
		)
)
