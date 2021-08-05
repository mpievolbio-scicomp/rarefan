# TODO: Add comment
# 
# Author: Carsten Fortmann-Grote
###############################################################################

source("analysis.R")
library(ggpubr)
# Parse command line args
args = commandArgs(trailingOnly=TRUE)
#folder: folder that contains presAbs_* files
#treeFile: name of newick tree file 
#type: REPIN type that is supposed to be viewed (_*, * is the type)

max_number_of_expected_args = 3
min_number_of_expected_args = 1
if (length(args)<min_number_of_expected_args ) { 
    stop("Usage: Rscript run_analysis.R DIR [TYPE [TREEFILE TYPE]]", call.=FALSE)
}
if (length(args)>max_number_of_expected_args ) { 
    stop("Usage: Rscript run_analysis.R DIR [TYPE [TREEFILE TYPE]]", call.=FALSE)
}
if (length(args) == 1) {
    data_dir=args[1]
    treefile="tmptree.nwk"
    rayt_type=0
}
if (length(args) == 2) {
    data_dir=args[1]
    treefile=args[2]
    rayt_type=0
}
if (length(args) == 3) {
    data_dir=args[1]
    treefile=args[2]
    rayt_type=args[3]
}


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

repin_facet_plot = plotREPINs(data_dir,
		   treefile,
		   rayt_type,
		   barcolor,
		   barsize,
		   fontsize 
)

phylogeny_plot = drawRAYTphylogeny(data_dir)

figure = ggarrange(phylogeny_plot, repin_facet_plot, correlation_plot,  ncol=1, nrow=2)

ggexport(figure, filename='rarefan_analysis.pdf')




