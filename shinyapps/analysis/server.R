#!/usr/bin/env Rscript
# Required libraries

suppressMessages(library(shiny))
suppressMessages(library(here))
suppressMessages(library(mongolite))

uploads_path <- here('app', 'static', 'uploads')
mongo_username <- "rarefan"
mongo_password <- ""
mongo_uri = sprintf("mongodb://%s:%s@localhost:27017",
                    URLencode(mongo_username, reserved=TRUE),
                    URLencode(mongo_password, reserved=TRUE)
)

db <- mongo(db='rarefan', collection='job', url="mongodb://localhost:27017")

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
              db_query <- paste0('{"run_id":"', query$run_id, '"}')

              # Get run setup from database, filter relevant columns and transpose for nicer display.
              db_run_setup <- data.frame(db$find(db_query)$setup)
              setup <- db_run_setup %>% select(reference_strain,
                                               query_rayt,
                                               nmer_length,
                                               min_nmer_occurrence,
                                               e_value_cutoff
                                               )
              if('distance_group_seeds' %in% colnames(db_run_setup)) {
                setup$distance_group_seeds <- db_run_setup$distance_group_seeds
              }
              else {
                setup$distance_group_seeds <- 30
              }

              # Change order
              setup <- setup %>% select(reference_strain,
                                               query_rayt,
                                               nmer_length,
                                               min_nmer_occurrence,
                                               distance_group_seeds,
                                               e_value_cutoff
              )
              setup <- t(setup)
              rownames(setup) <- c("Ref. strain",
                                         "Query RAYT",
                                         "Seed length",
                                         "Min. seed count",
                                         "Dist. group seeds",
                                         "e-value cutoff"
              )
              output$text <- renderText({
          				paste("Run ID ", query$run_id, sep=" ")
          			}
              )
              output$run_setup <- renderTable(setup, rownames = T, colnames = F)
              output$back_to_results <- renderUI({
                                                  tags$a("Back to results summary page",
                                                         href=paste("http://rarefan.evolbio.mpg.de/results?run_id=",
                                                         query$run_id,
                                                         sep=""
                                                         )
                                                         )
                                                         }
                                                         )
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
              rayt_phylogeny <- drawRAYTphylogeny(out_dir)
              number_of_strains <- dim(rayt_phylogeny$data)[[1]]
              plot_height <- max(600, ceiling(number_of_strains * fontsize/4 * 1.2))
              plot_width = ceiling(1.41 * plot_height)
              logging::logdebug(paste0("plot width: ", plot_width))
              logging::logdebug(paste0("plot height: ", plot_height))
              output$rayt_tree <-  renderPlot({rayt_phylogeny},
                                              height=plot_height,
                                              width=plot_width
              )

              logging::logdebug("Calling 'plotREPINs'.")
              output$repin_tree <-  renderPlot({
                plotREPINs(out_dir,
                           treefile,
                           input$rayt
                )
              },
              height=plot_height, width=plot_width
              )

              logging::logdebug("Calling 'plotCorrelationSingle'.")
              output$correlations <-  renderPlot({
                plotCorrelationSingle(out_dir,
                                      input$rayt
                )
              },
              height=600, width=846
              )
            })
}
