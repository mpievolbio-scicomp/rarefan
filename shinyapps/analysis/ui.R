suppressMessages(library(shiny))

fluidPage(
		
		# App title
		titlePanel("REPIN and RAYT analysis"),
		sidebarLayout(
				sidebarPanel(
					textOutput("text"),
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
								selected = 0),
						textOutput("plot_instructions")
						),
				mainPanel(
					h4("RAYT tree"),
					plotOutput(outputId = 'rayt_tree', width="60%", height="400px" ),
					hr(),
					h4("REPINs"),
					plotOutput(outputId = 'repin_tree', width="60%", height="400px" ),
					hr(),
					h4("Correlation"),
					plotOutput(outputId = 'correlations', width="60%", height="400px" )
				)
		)
)
