suppressMessages(library(shiny))

fluidPage(
		
		# App title
		titlePanel("REPIN and RAYT analysis"),
		sidebarLayout(
				sidebarPanel(
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
								selected = 0),
					hr(),
					textOutput("These plots can be reproduced by running the R scripts in the results zip file downloaded from the results page.
					The syntax is \n
					\n
						Rscript run_analysis.R [-r RAYT] [-t TREEFILE] DIR OUTFILE \n
					\n
					where DIR is the path to the output directory in the unzipped results dataset and OUTFILE is the file where to save the figure.")
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
