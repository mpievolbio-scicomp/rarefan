#!/usr/bin/env Rscript
# Required libraries

suppressMessages(library(shiny))

# Include local definitions
source("analysis.R")

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
              logging::logdebug(session$clientData$url_search)
              logging::logdebug("Still alive")
              run_dir <- paste0("/home/rarefan/repinpop/app/static/uploads/", query$run_id)
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
                           input$rayt,
                           "#40e0d0",
                           2,
                           fontsize
                )
              })

              logging::logdebug("Calling 'plotCorrelationSingle'.")
              output$correlations <-  renderPlot({
                plotCorrelationSingle(out_dir,
                                      input$rayt,
                                      c(0,1),
                                      c(0,320),
                                      theme,
                                      fontsize,
                                      "left",
                                      "bottom"
                )
              })
            })
}
