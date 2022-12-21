# TODO: Add comment
#
# Author: Carsten Fortmann-Grote, Frederic Bertels
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
suppressMessages(library(glue))
suppressMessages(library(ggtext))
suppressMessages(library(RColorBrewer))

# Set log level
logging::basicConfig()
logging::setLevel(30) # 10: debug, 20: info, 30: warning, 40: error

# Colors
colors = c(brewer.pal(8, "Dark2"), brewer.pal(9, 'Set1'))

get_color = function(index) {
  if(index>=0 && index < 9) {
    return(colors[index])
  }
  return("black")
}

# Set theme for all plots
rarefan_theme=theme(axis.line.x = element_line(colour = "black"),
            legend.key = element_rect(fill = "white"),
            axis.line.y = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.justification = c(0, 1),
            legend.position = c(0.80, 1),
            legend.text = element_text(hjust=0),
            panel.spacing=unit(2,"lines"),
	          legend.title=element_blank()
)

# Set fontsize globally.
fontsize=16

blank_theme = theme(axis.text=element_text(size=fontsize),
                  axis.line = element_blank(),
                  axis.line.x.bottom = element_blank(),
                  axis.line.y.left = element_blank(),
                  axis.title = element_blank(),
                  axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  text=element_text(size=fontsize)
            )



######################################################################################
# Plot phylo tree for all input NA sequences, number of RAYTs and number of REPs for each species.
plotREPINs=function(folder,
                    treeFile,
                    rep_rayt_group,
                    highlight_strain="",
                    analyse_repins="y"){

  ### Process tree file
  tree_file = paste0(folder,"/",treeFile)
  logging::logdebug("Attempting to read tree file.")

  # Attempt to open the tree file. If not readable (or not existing), return empty plot.
  tree=tryCatch(read.tree(tree_file),
		  error=function(e)
		  logging::logwarn("Tree file %s does not exist or is not readable.", tree_file)
		  )

  logging::logdebug(tree)
  tree_file_is_corrupt = typeof(tree) == "logical"
  if(tree_file_is_corrupt) {
    p = ggplot() +
            geom_blank() +
            xlim(c(0, 1)) +
            ylim(c(0,1)) +
            ggtitle("No REPINs found.") +
            theme(axis.text=element_text(size=fontsize),text=element_text(size=fontsize)) + rarefan_theme

		return(p)
  }

  # Ok, we have read the tree, now render the plot.
  logging::loginfo("Plotting ggtree...")

  # First facet is the tree itself.
  p <- ggtree(tree)+
      scale_x_continuous(breaks=scales::pretty_breaks(n=3))+
      geom_tiplab(size=fontsize*1/4)

  # Set limits.
  p <- p+xlim_tree(layer_scales(p)$x$get_limits()[2]*2)

  # RAYT population size.
  # Read the REP-RAYT association table.
  assoc_file <-  paste0(folder,"/repin_rayt_association_byREPIN.txt")
  logging::logdebug("Read association data fom %s.", assoc_file)
  association <- read.table(assoc_file,header=TRUE)

  # Take out the rows for the queried rep_rayt_group (alias repintype)
  d <- association[association$repintype==rep_rayt_group,]

  # Setup path to repin_rayt association table file.
  repin_rayt_assoc_table_file = paste0(folder,"/repin_rayt_association.txt")

  # Setup colors
  colorDF = determineColor(repin_rayt_assoc_table_file)

  # Fix the rayt color.
  rayt_color = colorDF[colorDF$repRAYT==rep_rayt_group,2]

  # If no color is defined, set it to grey.
  if(length(rayt_color) == 0) {
      rayt_color=c("grey")
  }
  logging::logdebug(paste0("rayt_color=", rayt_color))

  # Calculate total number of rayts
  num_rayts = sum(d$rayts)
  logging::logdebug("Number of rayts: %d", num_rayts)

  # REPIN population size.
  logging::logdebug("Constructing popSize data.")
  data_file = paste0(folder,"/presAbs_",rep_rayt_group,".txt")
  logging::logdebug("Reading table from %s.", data_file)

  # Attemt to read presAbs file. If file does not exist or cannot be read,
  # issue a warning.
  t=tryCatch(read.table(data_file,sep="\t", skip=1),
		  error=function(e)
		  logging::logwarn("File %s is empty.", data_file)
          )

  # Check if file was read.
  data_file_is_corrupt = typeof(t) == "logical"

  if(!data_file_is_corrupt) {
    if(num_rayts > 0){
        p <- facet_plot(p
                        , panel='RAYTs'
                        , data=d
                        , geom=geom_segment
                        , aes(x=0
                              , xend=rayts
                              , y=y
                              , yend=y)
                        , size=2
                        , color=unique(rayt_color)
                        )
    }

    # If no rayts found
    else {
        logging::logwarn("No RAYTs in association table %s for RAYT rep_rayt_group %s.", repin_rayt_assoc_table_file, rep_rayt_group)

    }

    # Construct new data frame holding the data to plot in facets.
    popSize=data.frame(name=t[,1],
                       rayts=t[,2],
                       repins=t[,9],
                       prop=(t[,5]/t[,9]),
                       propAll=(t[,9]/t[,6]),
                       numClus=t[,7],
                       diffRAYTCluster=t[,7]-t[,2])

    logging::logdebug("popSize=%s", str(popSize))
    # Add repin population size.

    if(sum(popSize$repins) > 0) {
    p = facet_plot(p
                   , panel='Largest REPIN population size'
                   , data=popSize
                   , geom=geom_segment
                   , aes(x=0
                         , xend=repins
                         , y=y
                         , yend=y)
                   , size=2
                   ,color=unique(rayt_color)
                   )

    }
  }
  else {
          if(is.na(analyse_repins)) {
              msg <- paste0("No REP sequences could be identified in Group ", rep_rayt_group, ".")
          }
          else if(analyse_repins == "y") {
              msg <- paste0("No REPIN sequences could be identified in Group ", rep_rayt_group, ".")
          }
          else {
              msg <- paste0("No REP sequences could be identified in Group ", rep_rayt_group, ".")
          }
          p <- p + ggtitle(msg)
          # annotate(x=0.0,
          #                   y=layer_scales(p)$y$get_limits()[2]*0.95,
          #                   geom='text',
          #                   label=msg
          #                    )
  }

  # Apply theme.
  p = p + theme_tree2()

  # Adjust theme.
  p = p +
        theme(strip.text.x=element_text(hjust=0),
              strip.background = element_rect(color="black",
                                              fill="white",
                                              linetype="blank"
              )
        )
  if(highlight_strain!="") {
  p = p +
      geom_highlight(node=which(tree$tip.label==highlight_strain),
                     extend=230,
                     colour='red',
                     fill='white',
                     size=0.3,
                     alpha=0.0)
  }
  # Font size.
  p = p +
          theme(text=element_text(size=fontsize),
                axis.text=element_text(size=fontsize)
                )

  return(p)
}

######################################################################################
determineColor=function(associationFile){
  logging::logdebug("Determine colors from %s.", associationFile)
  ass=read.delim(associationFile,header=TRUE)
  logging::logdebug(str(ass))
  logging::logdebug(ass)

  colorAss=c()

  for(i in 1:length(ass[,1])){
     rayt=paste0(ass[i,1],"_",ass[i,2])
     c="NA"

     if(typeof(ass$REPINgroups) == 'logical') {
         logging::logdebug("REPINgroups is empty, set color to 'grey'.")
    		 c='grey'
     }
     else if(!is.na(ass[i,3])&&nchar(ass[i,3])>0){
        split0=str_split(ass[i,3],",")
        position = as.integer(split0[[1]][1])+1
        c=get_color(position)
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
     if(!is.na(i)&&nchar(i)>0){
        split=str_split(i,",")
        for(j in split[[1]]){
           j=as.integer(j)
           pos=as.integer(split[[1]][1])
           if(length(colorAss[colorAss$repRAYT==j,]$color)==0){
             temp=data.frame(repRAYT=j,color=get_color(pos+1))
             colorAss=rbind(colorAss,temp)
		   logging::logdebug(paste0("i=",i, " j=",j, " pos=", pos, " temp=", temp))
           }
        }
     }
  }
  logging::logdebug(paste0("colorAss=", colorAss))
  return(colorAss)
}

######################################################################################
plotCorrelationSingle=function(folder,
                               rep_rayt_group){

	logging::logdebug("Plotting correlation.")

	data_file = paste0(folder,"/presAbs_",rep_rayt_group,".txt")
	logging::logdebug("Reading data from %s.", data_file)
	t=tryCatch(read.table(data_file,sep="\t", skip=1),
		  error=function(e)
		  logging::logwarn("File %s is empty.", data_file)
	)

	if(typeof(t) == "logical") {
        p = ggplot() +
            geom_blank() +
            xlim(c(0, 1)) +
            ylim(c(0,1)) +
            ggtitle("No data to correlate.") +
            blank_theme

		return(p)
	}

  t$propMaster=t[,5]/t[,9]
  t$numRepin=t[,9]

	assoc_file = paste0(folder,"/repin_rayt_association_byREPIN.txt")
	logging::logdebug("Reading association data from %s.", assoc_file)
  association=read.table(assoc_file,header=TRUE)

	logging::logdebug("Preparing data structures and colors.")
  association=association[association$repintype==rep_rayt_group,]
  t$numRAYT=association[match(t[,1],association[,1]),]$rayts
  t$RAYTs = as.factor(t$numRAYT)

  # Get color dataframe.
	colorDF = determineColor(paste0(folder,"/repin_rayt_association.txt"))

	# Color for current REP/RAYT group.
  rayt_color = unique(colorDF[colorDF$repRAYT==rep_rayt_group,]$color)

  # Sort in descending order of numRAYT.
  t <- t[order(t$numRAYT, decreasing=T),]

  legend_values = unique(t[t$numRAYT>0,]$numRAYT)
  legend_values = legend_values[order(legend_values)]
  # browser()
  # Setup the plot
  p <- ggplot(t) +
    # Plot only observations with numRAYT>0 and map size to number of RAYTs (as factor).
    geom_point2(
           aes(x=propMaster,
               y=numRepin,
               size=RAYTs,
               subset=numRAYT>=1
           ),
           color=rayt_color,
       ) +
    # Now add the obs. with numRAYT==0 and set shape (will be fixed manually).
    geom_point2(
           aes(x=propMaster,
               y=numRepin,
               shape=RAYTs,
               subset=numRAYT==0,
           ),
           size=3
       ) +
    # Set correct size of dots.
    scale_size_manual(
      values=legend_values*2 + 1
    ) +
    # Set the shape of numRAYT=0 observations.
    scale_shape_manual(
      values=c(1),
      labels=c("No RAYT"),
      guide=guide_legend(override.aes = list(shape=1))
      )+
    # Legend titles.
    labs(size="RAYTs", shape=NULL) +
    # Set legend order.
    guides(size=guide_legend(order=1), shape=guide_legend(order=2))

  logging::logdebug("Adding limits, theme, and axis labels.")
  p <- p +
        xlim(c(0,1)) + ylim(c(0,max(t$numRepin+0.1*t$numRepin)))+
        rarefan_theme +
        xlab("Proportion master sequence (~Replication rate)") +
        ylab("Largest REPIN population size")

	logging::logdebug("Adding theme.")
  p <- p + theme(axis.text=element_text(size=fontsize),text=element_text(size=fontsize))

	logging::logdebug("Done, return from 'plotCorrelations'.")

  return(p)
}

######################################################################################
get_rayt_phylogeny=function(data_dir){

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
	  aln=muscle(raytseqs, verbose=T, log='/tmp/muscle.log')
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
    phyml_cmd = sprintf("phyml --quiet -i %s -m GTR", raytAlnFile)
	  system(phyml_cmd)
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
drawRAYTphylogeny=function(data_dir, reference_strain=""){

  # Check and get phylogeny data.
  rayt_files = get_rayt_phylogeny(data_dir)
  if(rayt_files$raytAlnFile == "") {
	 p <- ggplot() +
            geom_blank() +
            xlim(c(0, 1)) +
            ylim(c(0,1)) +
            ggtitle("No RAYT sequences could be identified.") +
            blank_theme

    return(p)
  }

  #
  logging::loginfo("Drawing RAYT phylogeny.")

  # Read tree file.
  raytTreeFile=rayt_files$raytPhyTreeFile

  nwk=read.tree(raytTreeFile)
  # If file is empty, return empty plot.

  if(is.null(nwk)) {
    p <-  ggplot() +
            geom_blank() +
            xlim(c(0, 1)) +
            ylim(c(0,1)) +
            ggtitle("No data in RAYT phylogeny.") +
            blank_theme

		return(p)
  }

  nwk$node.label = c(1:nwk$Nnode)

  # Get tree object.
  p <- ggtree(nwk)

  # Get colors
  colorDF = determineColor(paste0(data_dir,"/repin_rayt_association.txt"))
  logging::logdebug(colnames(colorDF))

  # Retain only RAYT colors.
  onlyRAYTs <- colorDF[colorDF[,1] %in% nwk$tip.label, ]
  logging::logdebug(colnames(onlyRAYTs))
  logging::logdebug(onlyRAYTs)

  # Add tip labels.
  p <- p %<+% onlyRAYTs + geom_tiplab(aes(color=color),size=fontsize*1/4)+theme_tree2()

  # Stretch x axis to accomodate tip labels.
  p=p+xlim(layer_scales(p)$x$get_limits()*1.5)

  ## Annotate reference strain.
  # Extract node ids that belong to the reference strain.
  # Check if reference strain annotation is requested.
  if(reference_strain != "") {
    reference_strain_node_ids = sapply(reference_strain, function(y) grep(y, nwk$tip.label))

    # Check if requested ref. strain is actually present.
    if (isFALSE(reference_strain_node_ids[[1]] == 0)) {
      label=nwk$tip.label[reference_strain_node_ids]

      # Temp data.frame to map tip label to ref. strain annotation (unformatted).
      tmp =  data.frame(label=label, x_label=reference_strain)

      # Now format the ref. strain label using html syntax.
      reference_annotation <- dplyr::mutate(tmp,
                                            label=label,
                                            x_label=glue("<b>&#10229;</b>"),
                                            )

      # Add ref strain annotation to ggtree object.
      p <- p %<+% reference_annotation + geom_richtext(data=td_filter(isTip),
                                                       aes(label=x_label),
                                                       label.color=NA,
                                                       nudge_x=0.8,
                                                       fill='orange',
                                                       alpha=0.7,
                                                      label.padding=grid::unit.c(grid::unit(0.0, "pt"),
                                                                                 grid::unit(8.0, "pt"),
                                                                                 grid::unit(-2.5, "pt"),
                                                                                 grid::unit(2.0, "pt")
                                                      )
      )
    }
  }

  # Add color tips.
  logging::logdebug("Added color tips")
  cols <- onlyRAYTs$color
  names(cols) <- onlyRAYTs$color
  logging::logdebug(cols)
  # Add colors
  p <-  p + scale_color_manual(values=cols,guide="none")

  logging::logdebug("Added color scale.")

  return(p)
}


