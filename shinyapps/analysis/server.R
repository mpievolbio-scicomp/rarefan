#!/usr/bin/env Rscript
# Required libraries

suppressMessages(library(shiny))
suppressMessages(library(here))

uploads_path <- here('app', 'static', 'uploads')
print(uploads_path)

# Include local definitions
source("analysis.R")

logging::addHandler(writeToFile, file="/tmp/R.log", level='DEBUG')

# Font size
fontsize=16

function(input, output, session) {
    logging::logdebug("Entering shiny app main()")
    observe(
            {
              query <- parseQueryString(session$clientData$url_search)
              output$text <- renderText({
          				paste("Run ID ", query$run_id, sep=" ")
          			}
              )
              output$back_to_results <- renderText({
                               a(h4("Back to results summary page"),
                                 href=paste("http://rarefan.evolbio.mpg.de/results/?run_id=",
                                            query$run_id
                                            )
                                 )
                                 }
                                 )
#			  output$plot_instructions <- renderText({
#          				paste("  ", "The plots can be reproduced with the R script 'run_analysis.R' which is part of the zip archive on the results page.", sep="")
#          			}
#			  )
              logging::logdebug(session$clientData$url_search)
              logging::logdebug("Still alive")
              run_dir <- paste0(uploads_path, '/', query$run_id)
              logging::logdebug(paste0("run_dir = ", run_dir))
              out_dir <- paste0(run_dir, "/out")
              logging::logdebug(paste0("out_dir = ", out_dir))
              treefile <- 'tmptree.nwk'
              logging::logdebug(paste0("treefile = ", treefile))

              logging::logdebug("Calling 'drawRAYTphylogeny'.")
              output$rayt_tree <-  renderPlot({
                drawRAYTphylogeny(out_dir)
              })

              logging::logdebug("Calling 'plotREPINs'.")
              output$repin_tree <-  renderPlot({
                plotREPINs(out_dir,
                           treefile,
                           input$rayt
                )
              })

              logging::logdebug("Calling 'plotCorrelationSingle'.")
              output$correlations <-  renderPlot({
                plotCorrelationSingle(out_dir,
                                      input$rayt
                )
              })
            })
}
