suppressMessages(library(shiny))
suppressMessages(library(plotly))

fluidPage(
		
		# App title
		titlePanel("REPIN and RAYT analysis"),
		sidebarLayout(
				sidebarPanel(
						textOutput("text"),
					hr(),
					selectInput(inputId = 'rayt',
								label="Select RAYT",
								choices = list(
										"RAYT 1" = 0,
										"RAYT 2" = 1,
										"RAYT 3" = 2,
										"RAYT 4" = 3,
										"RAYT 5" = 4,
										"RAYT 6" = 5
								),
								selected = 0)
							),
				mainPanel(
					h4("RAYT tree"),
					plotlyOutput(outputId = 'rayt_tree'),
					hr(),
					h4("REPINs"),
					plotlyOutput(outputId = 'repin_tree'),
					hr(),
					h4("Correlation"),
					plotlyOutput(outputId = 'correlations')
				)
		)
)
