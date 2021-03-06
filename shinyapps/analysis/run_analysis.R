#! /usr/bin/env R

suppressMessages(library(here))
r_src_dir <- here('shinyapps/analysis')
source(paste0(r_src_dir, "/analysis.R"))

packages=c("ggpubr", "argparse", "filenamer")

package_check = lapply(
                      packages,
                      FUN = function(x) {
                          if(!require(x, character.only = TRUE)) {
                              install.packages(x, dependencies=TRUE, repos="https://ftp.gwdg.de/pub/misc/cran")
                              suppressMessages(library(x, character.only = TRUE))
                          }
                      }
                  )

parser <- ArgumentParser(description='Create plots from rarefan results.')

# Input data dir
parser$add_argument('-d', '--data_dir',
                    metavar='DIR',
                    type="character",
					default='.',
                    help='The output data dir of the rarefan run to analyse. Defaults to current working directory.'
)

# RAYT index
parser$add_argument('-r', '--rayt',
                    metavar='RAYT',
                    type="integer",
                    dest='rayt_type',
                    default=0,
                    choices=c(0,1,2,3,4,5,6,7),
                    help='The RAYT index to calculate results for.'
)

# tree file
parser$add_argument('-t', '--tree',
                    metavar='TREEFILE',
                    type="character",
                    dest='treefile',
                    default='tmptree.nwk',
                    help='The name of the treefile to use for plotting the genome phylogenies.'
)

# tree file
parser$add_argument('-o', '--outfile',
                    metavar='OUTFILE',
					default='analysis.pdf',
                    type="character",
                    help='Save figures to OUTFILE'
)

args = parser$parse_args(commandArgs(TRUE))

data_dir = args$data_dir
rayt_type = args$rayt_type
treefile = args$treefile
outfile = args$outfile

# Add rayt type to outfile basename
fname = as.filename(outfile)
outfile = as.character(insert(fname, paste0("rayt",rayt_type, "_")))

logging::loginfo(paste0("Reading data from ", data_dir))
logging::loginfo(paste0("RAYT index = ", rayt_type))
logging::loginfo(paste0("treefile = ", treefile))
logging::loginfo(paste0("outfile = ", outfile))
# Parse command line args

correlation_plot = plotCorrelationSingle(
		data_dir,
		rayt_type
		)

repin_facet_plot = plotREPINs(
       data_dir,
		   treefile,
		   rayt_type
	)


phylogeny_plot = drawRAYTphylogeny(data_dir)

figure = ggarrange(phylogeny_plot, repin_facet_plot, correlation_plot,  ncol=1)

ggexport(figure, filename=outfile ,height=30)

