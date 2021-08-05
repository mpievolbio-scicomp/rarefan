#! /usr/bin/env R

source("analysis.R")
library(ggpubr)
library(argparse)

parser <- ArgumentParser(description='Create plots from rarefan results.')

# Input data dir
parser$add_argument('data_dir',
                    metavar='DIR',
                    type="character", 
                    help='The output data dir of the rarefan run to analyse.'
)

# RAYT index
parser$add_argument('-r', '--rayt',
                    metavar='RAYT',
                    type="integer", 
                    dest='rayt_type',
                    default=0,
                    choices=c(0,1,2,3,4,5),
                    help='The RAYT index to calculate results for.'
)

# tree file
parser$add_argument('-t', '--tree',
                    metavar='TREEFILE',
                    type="character", 
                    dest='treefile',
                    default='tmptree.nwk',
                    help='The treefile to use (default "DIR/tmptree.nwk")'
)

args = parser$parse_args(commandArgs(TRUE))

data_dir = args$data_dir
rayt_type = args$rayt_type
treefile = args$treefile

logging::loginfo(paste0("Reading data from ", data_dir))
logging::loginfo(paste0("RAYT index = ", rayt_type))
logging::loginfo(paste0("treefile = ", treefile))
# Parse command line args

outfile = paste0("rarefan_", rayt_type, ".pdf")
barcolor = "#40e0d0"
barsize = 2
fontsize = 12

correlation_plot = plotCorrelationSingle(
		folder=data_dir,
		type=rayt_type,
		xlim=c(0, 1),
		ylim=c(0, 320),
		theme=theme,
		fontsize=fontsize,
		pvLabelX='left',
		pvLabelY='bottom'
		)

repin_facet_plot = plotREPINs(
           data_dir,
		   treefile,
		   rayt_type,
		   barcolor,
		   barsize,
		   fontsize 
)

phylogeny_plot = drawRAYTphylogeny(data_dir)

figure = ggarrange(phylogeny_plot, repin_facet_plot, correlation_plot,  ncol=1, nrow=1)

ggexport(figure, filename=outfile)




