# TODO: Add comment
# 
# Author: Carsten Fortmann-Grote
###############################################################################

# Load libraries silently.
suppressMessages(library(Biostrings))
suppressMessages(library(ape))
suppressMessages(library(muscle))
suppressMessages(library(ggpmisc))
suppressMessages(library(ggtree))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(logging))
suppressMessages(library(roxygen2))


# Set log level
logging::basicConfig()
logging::setLevel(10) # 10: debug, 20: info, 30: warning, 40: error 

# 6 Colors for plots (corresponding to 6 RAYT types)
colors=c("#45BA55", "#5545BA", "#BA5545", "#B6BD42", "#42B6BD", "#BD42B6")

# Set theme for all plots
logging::logdebug("defining theme")
theme=theme(axis.line.x = element_line(colour = "black"),
            legend.key = element_rect(fill = "white"),
            axis.line.y = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.justification = c(0, 1),
            legend.position = c(0.5, 1),
            legend.title=element_blank(),
            legend.text = element_text(hjust=0),
            panel.spacing=unit(2,"lines"),
)


######################################################################################
plotREPINs=function(folder,treeFile,type,colorBars,bs,fontsize){
#' @description Plot the REPIN phylogeny, RAYT population size and REPIN population size
#' 
#' @param folder Output directory containing the RAREFAN results (character).
#' @param treeFile Name of the treefile in the output directory (character).
#' @param type The RAYT type to analyse (int in range 0-5)
#' @return A ggplot object holding the plot.
#'
#'  @examples
#'  \dontrun{
#'  plotREPINs('/tmp/rarefan_test/neisseria/out', 
#'               'neisseria.nwk',
#'               0,
#'               'green',
#'               2,
#'               12
#'               )
#' }
        

  logging::logdebug("Enter function 'plotREPINs' with ")
  logging::logdebug(paste0("    folder = ", folder))
  logging::logdebug(paste0("    treeFile = ", treeFile))
  logging::logdebug(paste0("    type = ", type))
  logging::logdebug(paste0("    colorBars = ", colorBars))
  logging::logdebug(paste0("    bs = ", bs))
  logging::logdebug(paste0("    fontsize = ", fontsize))

  # Set theme.
  themeCurr=theme(axis.line.x = element_line(colour = "black"),
                  legend.key = element_rect(fill = "white"),
                  axis.line.y = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  legend.position="none",
                  panel.spacing=unit(2,"lines"),
                  plot.margin = unit(c(5.5,12,5.5,10.5), "pt"),
  )

  ### Process tree file
  tree_file = paste0(folder,"/",treeFile)
  logging::logdebug("Reading tree file %s.", tree_file)

  tree=read.tree(tree_file)
  tips=tree$tip.label
  logging::logdebug(str(tree))

  logging::logdebug("Plotting ggtree...")
  p=ggtree(tree)+
      scale_x_continuous(breaks=scales::pretty_breaks(n=3))+
#      xlim_tree(0.2)+
      geom_tiplab()

  # RAYT population size.
  assoc_file = paste0(folder,"/repin_rayt_association_byREPIN.txt")
  logging::logdebug("Read association data fom %s.", assoc_file)
  
  association=read.table(assoc_file,header=TRUE)
  logging::logdebug(colnames(association))
  logging::logdebug(association)

  d <- association[association$repintype==type,]
  logging::logdebug(str(d))
  logging::logdebug(d)
  
  repin_rayt_assoc_table_file = paste0(folder,"/repin_rayt_association.txt")

  colorDF = determineColor(repin_rayt_assoc_table_file)
  logging::logdebug(str(colorDF))
  logging::logdebug(colorDF)

  rayt_color = colorDF[colorDF$repRAYT==type,2]
  logging::logdebug(paste0("typeof(rayt_color)=", typeof(rayt_color)))
  logging::logdebug(paste0("length(rayt_color)=", length(rayt_color)))

  # In some cases, we get a list of two identical colors but need only one.
  if (length(rayt_color) > 1) {
      rayt_color = rayt_color[[1]]
  }

  logging::logdebug(paste0("typeof(rayt_color)=", typeof(rayt_color)))
  logging::logdebug(paste0("length(rayt_color)=", length(rayt_color)))
  logging::logdebug(paste0("rayt_color=", rayt_color))

  # In other cases, no color is defined, set it to grey.
  if(length(rayt_color) == 0) {
      rayt_color=c("grey")
  }
  logging::logdebug(paste0("rayt_color=", rayt_color))

  logging::logdebug("Tree tip labels = %s", tree$tip.label)
  logging::logdebug("repin_rayt_association.txt['genome'] = %s", d$genome)
  logging::logdebug("repin_rayt_association.txt['rayts'] = %s", d$rayts)

  num_rayts = length(d$rayts)
  logging::logdebug("Number of rayts: %d", num_rayts)
  if(num_rayts > 0){
      logging::logdebug("Plotting RAYTs.")
      p <- facet_plot(p,
                    panel='RAYTs',
                    data=d,
                    geom=geom_segment,
                    aes(x=0,
                        xend=rayts,
                        y=y,
                        yend=y),
                    size=bs,
                    color=rayt_color
          )
  }
  else {
      logging::logwarn("No RAYTs in association table %s for RAYT type %d.", repin_rayt_assoc_table_file, type)
      logging::logwarn("Not plotting RAYT population size.")
  }

  # REPIN population size.
  logging::logdebug("Constructing popSize data.")
  data_file = paste0(folder,"/presAbs_",type,".txt")
  logging::logdebug("Reading table from %s.", data_file)

  # data_file_ok = FALSE
  t=tryCatch(read.table(data_file,sep="\t", skip=1),
		  error=function(e)
		  logging::logwarn("File %s is empty.", data_file)
          )

  logging::logdebug(typeof(t))
  logging::logdebug(t)

  data_file_is_corrupt = typeof(t) == "logical"
  if(!data_file_is_corrupt) {
      popSize=data.frame(name=t[,1],
                         rayts=t[,2],
                         repins=t[,3],
                         prop=(t[,5]/t[,3]),
                         propAll=(t[,3]/t[,6]),
                         numClus=t[,7],
                         diffRAYTCluster=t[,7]-t[,2])

      logging::logdebug(str(popSize))

      logging::logdebug("plotting repin population size.")
      p = facet_plot(p,
                    panel='repin population size',
                    data=popSize,
                    geom=geom_segment,
                    aes(x=0,
                        xend=repins,
                        y=y,
                        yend=y),
                    size=bs
                   ,color=rayt_color
                     )
      
  }
  else {
          logging::logwarn("Not plotting REPIN population size for RAYT type %d", type)
      }
  
  logging::logdebug("adding theme.")
  p = p + theme_tree2()

  logging::logdebug("customizing theme.")
  p = p +
        theme(strip.text.x=element_text(hjust=0),
              strip.background = element_rect(color="black",
                                              fill="white",
                                              linetype="blank"
              )
        )

  logging::logdebug("setting font size.")
  p = p +
          theme(text=element_text(size=fontsize)) +
          themeCurr

  logging::logdebug("Done, returning from function 'plotREPIN'")

  return(p)
}

######################################################################################
determineColor=function(associationFile){
    #' @description Determine the color for plotting the REPIN data based on the RAYT type
    #'
    #' @param associationFile Filename (path) of the data file that associates REPINs and RAYT types.
    #' Typically, this is the file 'repin_rayt_association.txt' or 'repin_rayt_association_byREPIN' 
    #' in the 'out/' directory of the run
    #' @return A string designating the color.

  logging::logdebug("Determine colors from %s.", associationFile)
  ass=read.delim(associationFile,header=TRUE)
  logging::logdebug(str(ass))
  logging::logdebug(ass)

  colorAss=c()

  for(i in 1:length(ass[,1])){
     rayt=paste0(ass[i,1],"_",ass[i,2])
     c="NA"

     if(typeof(ass$REPINgroups) == 'logical') {
		 c='grey'
     }
     else if(nchar(ass[i,3])>0){
        split0=str_split(ass[i,3],",")
        c=colors[as.integer(split0[[1]][1])+1]
     }
	 else{
        c="grey"
     }
     temp=data.frame(repRAYT=rayt,color=c)
     colorAss=rbind(colorAss,temp)

  }
  if(typeof(ass$REPINgroups) == 'logical') {
      return(colorAss)
  }
  groups=unique(ass[,3])
  logging::logdebug(paste0("groups=", groups))
  for(i in groups){
     if(nchar(i)>0){
        split=str_split(i,",")
        for(j in split[[1]]){
           pos=as.integer(split[[1]][1])
           j=as.integer(j)
           temp=data.frame(repRAYT=j,color=colors[pos+1])
           colorAss=rbind(colorAss,temp)
		   logging::logdebug(paste0("i=",i, " j=",j, " pos=", pos, " temp=", temp))
        }
     }
  }
  logging::logdebug(paste0("colorAss=", colorAss))
  return(colorAss)
}

######################################################################################
plotCorrelationSingle=function(folder,type,
                               xlim,
                               ylim,
                               theme,
                               fontsize,
                               pvLabelX,
                               pvLabelY,
                               subsetSmooth=F,
                               from=F,
                               to=F,
                               repinThreshold=0,
                               name=F,
                               labelOdd=""){
    #'
    #' @description Plot the REPIN population size vs. RAYT population size as a scatter plot.
    #'
    #' @param folder Name (path) of the output directory generated by RAREFAN (character).
    #' @param xlim x-axis range to plot (vector).
    #' @param ylim y-axis range to plot (vector).
    #' @param theme The ggplot theme to apply. (ggplot2::theme)
    #' @param fontsize Fontsize for plots (int).
    #' @param pvLabelX X position indicator for plot labels (character).
    #' @param pvLabelX Y position indicator for plot labels (character).
    #' @return A ggplot object holding plot.
    #' @example
    #' \dontrun{
    #' plotCorrelationSingle('/tmp/rarefan_test/neisseria/out',
    #'                       0,
    #'                       c(0, 1),
    #'                       c(0, 320),
    #'                       theme(),
    #'                       12,
    #'                       "left",
    #'                       "top"
    #'                       )
    #' }
                             

	logging::logdebug("Plotting correlation.")

	data_file = paste0(folder,"/presAbs_",type,".txt")
	logging::logdebug("Reading data from %s.", data_file)
	t=tryCatch(read.table(data_file,sep="\t", skip=1),
		  error=function(e)
		  logging::logwarn("File %s is empty.", data_file)
	)

	if(typeof(t) == "logical") {
		return(ggplot())
	}

    t$propMaster=t[,5]/t[,9]
    t$numRepin=t[,9]

	assoc_file = paste0(folder,"/repin_rayt_association_byREPIN.txt")
	logging::logdebug("Reading association data from %s.", assoc_file)
    association=read.table(assoc_file,header=TRUE)

	logging::logdebug("Preparing data structures and colors.")
    association=association[association$repintype==type,]
    t$color=association[match(t[,1],association[,1]),]$rayts
    cols=t$color
    names(cols)=cols
	colorDF = determineColor(paste0(folder,"/repin_rayt_association.txt"))
    cols[cols>0]=colorDF[colorDF$repRAYT==type,]$color
    cols[cols==0]="black"

	logging::logdebug("Setting up ggplots.")
    p=ggplot(t,
             aes(x=propMaster,
                 y=numRepin,
                 col=factor(color)
             )
    ) +
      scale_color_manual(values=cols,
                         labels=c("no RAYT", "RAYT"),
                         guide="legend"
                         ) +
      geom_point()

	logging::logdebug("Adding limits, theme, and axis labels.")
    p=p+xlim(xlim)+
        ylim(ylim)+
        theme +
        xlab("Proportion master sequence (~Replication rate)")+
        ylab("REPIN population size")

	logging::logdebug("Adding theme.")
    p=p+theme(axis.text=element_text(size=fontsize),text=element_text(size=fontsize))

	logging::logdebug("Done, return from 'plotCorrelations'.")
    return(p)
}

######################################################################################
get_rayt_phylogeny=function(data_dir){
    #' 
    #' @description Read or compute the RAYT phylogeny
    #' @details Scans the given directory (data_dir) for a file named 'raytAln.phy' 
    #'          resulting from a muscle alignment run.
    #'          If not found, muscle is run on the RAYT aminoacid sequence data found in the given directory.
    #'          Then checks for files named 'raytAln_phyml_tree.txt' and 'raytAln_phyml_stats.txt'. If either
    #'          of these files is missing, `phyml` will be run on 'raytAln.phy'. A `structure` that holds all
    #'          three filenames will be returned.
    #' @param data_dir Name (path) of the output directory of the RAREFAN run.
    #' @return A `structure` holding the filenames for the phylogeny data.

  logging::loginfo("Getting RAYT phylogeny.")
 
  # Check if phylogeny exists.
  raytAlnFile = paste0(data_dir,"/raytAln.phy")
  
  ### TODO: use filenamer methods here.
  raytPhyTreeFile = paste0(raytAlnFile, "_phyml_tree.txt")
  raytPhyStatsFile = paste0(raytAlnFile, "_phyml_stats.txt")

  files = structure(list(raytAlnFile="", raytPhyTreeFile="", raytPhyStatsFile=""))
  if(!file.exists(raytAlnFile)) {
	  
	  logging::loginfo(paste0("RAYT alignment file '", raytAlnFile, "' not found, will perform alignment now."))
	  # The associated rayt sequence file.
	  raytseqFile=paste0(data_dir,"/repin_rayt_association.txt.fas")
	  
	  # Read sequences
	  raytseqs=readDNAStringSet(raytseqFile,format="fasta")
	  
	  # If no rayt sequences found, return empty.
	  if(length(raytseqs) == 0) {
		  logging::logwarn(paste0("No RAYT sequences found in ", raytseqFile))
		  return(files)
	  }

	  
	  # Run muscle alignment.
	  logging::logdebug("Running muscle...")
	  aln=muscle(raytseqs)
	  logging::logdebug("...done.")
	  
	  # Write alignment to file.
	  write.phylip(aln,raytAlnFile)
  }
  else {
	  logging::loginfo(paste0("RAYT alignment file '", raytAlnFile, "' found."))
  }
	  
  if(!(file.exists(raytPhyTreeFile) && file.exists(raytPhyStatsFile))) {
	  logging::loginfo(paste0("RAYT phylogeny data files not found, will compute phylogeny now."))
	  # Run phyml as system command.
	  system(paste0("phyml --quiet -i ",raytAlnFile," -m GTR"))
  }
  else {
	  logging::loginfo(paste0("RAYT phylogeny data files found."))
  }
  
  files$raytAlnFile = raytAlnFile
  files$raytPhyTreeFile=raytPhyTreeFile
  files$raytPhyStatsFile=raytPhyStatsFile

  return(files)

}


######################################################################################
drawRAYTphylogeny=function(data_dir){
    #'
    #' @description Draw the RAYT phylogeny as a tree using the ggtree R package.
    #'
    #' @param data_dir Directory name (path) from where to read the phylogeny data.
    #'
    #' @return A ggplot object holding the tree graph.

  # Check and get phylogeny data.
  rayt_files = get_rayt_phylogeny(data_dir)
  if(rayt_files$raytAlnFile == "") {
	  logging::logwarn("Alignment is empty, will return empty plot.")
	  return(ggplot())
  }
  
  #
  logging::loginfo("Drawing RAYT phylogeny.")
 	
  # Read tree file.
  raytTreeFile=rayt_files$raytPhyTreeFile
  nwk=read.tree(raytTreeFile)

  # Plot the tree
  p <- ggtree(nwk)
	
  # Get colors
  colorDF = determineColor(paste0(data_dir,"/repin_rayt_association.txt"))
  logging::logdebug(colnames(colorDF))
  logging::logdebug(colorDF)
  
  # Retain only RAYT colors.
  onlyRAYTs <- colorDF[colorDF[,1] %in% nwk$tip.label, ]
  logging::logdebug(colnames(onlyRAYTs))
  logging::logdebug(onlyRAYTs)
  
  # Add tip labels.
  p <- p %<+% onlyRAYTs + geom_tiplab(aes(color=color))
  logging::logdebug("Added color tips")
  cols <- onlyRAYTs$color
  names(cols) <- onlyRAYTs$color
  logging::logdebug(cols)
  
  # Add colors
  p <-  p + scale_color_manual(values=cols,guide="none")
  
  logging::logdebug("Added color scale.")

  return(p)
}


