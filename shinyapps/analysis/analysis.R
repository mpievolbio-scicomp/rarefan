# TODO: Add comment
# 
# Author: Carsten Fortmann-Grote
###############################################################################

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

colors=c("#45BA55", "#5545BA", "#BA5545", "#B6BD42", "#42B6BD", "#BD42B6")
# Define plot routine

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


plotREPINs=function(folder,treeFile,type,colorBars,bs,fontsize){

  logging::logdebug("Enter function 'plotREPINs' with ")
  logging::logdebug(paste0("    folder = ", folder))
  logging::logdebug(paste0("    treeFile = ", treeFile))
  logging::logdebug(paste0("    type = ", type))
  logging::logdebug(paste0("    colorBars = ", colorBars))
  logging::logdebug(paste0("    bs = ", bs))
  logging::logdebug(paste0("    fontsize = ", fontsize))

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
  assoc_file = paste0(folder,"/repin_rayt_association_byREPIN.txt")
  logging::logdebug("Read association data fom %s.", assoc_file)
  association=read.table(assoc_file,header=TRUE)

  logging::logdebug(colnames(association))
  logging::logdebug(association)

  data_file = paste0(folder,"/presAbs_",type,".txt")
  logging::logdebug("Reading table from %s.", data_file)
  t=tryCatch(read.table(data_file,sep="\t", skip=1),
		  error=function(e)
		  logging::logwarn("File %s is empty.", data_file)
          )

  logging::logdebug(typeof(t))
  logging::logdebug(t)
	if(typeof(t) == "logical") {
        logging::logdebug("Returning empty plot.")
		return(ggplot())
	}

  logging::logdebug("Constructing popSize data.")
  popSize=data.frame(name=t[,1],
                     rayts=t[,2],
                     repins=t[,3],
                     prop=(t[,5]/t[,3]),
                     propAll=(t[,3]/t[,6]),
                     numClus=t[,7],
                     diffRAYTCluster=t[,7]-t[,2])


  tree_file = paste0(folder,"/",treeFile)
  logging::logdebug("Reading tree file %s.", tree_file)
  tree=read.tree(tree_file)
  tips=tree$tip.label

  logging::logdebug("Getting ggtree object.")
  p=ggtree(tree)+
      scale_x_continuous(breaks=scales::pretty_breaks(n=3))+
#      xlim_tree(0.2)+
      geom_tiplab()

  logging::logdebug("Plotting RAYTs.")
  colorDF = determineColor(paste0(folder,"/repin_rayt_association.txt"))
  logging::logdebug(str(colorDF))
  logging::logdebug(colorDF)
  d <- association[association$repintype==type,]
  logging::logdebug(str(d))
  logging::logdebug(d)
  rayt_color = unique(colorDF[colorDF$repRAYT==type,2])

  logging::logdebug(paste0("typeof(rayt_color)=", typeof(rayt_color)))
  logging::logdebug(paste0("length(rayt_color)=", length(rayt_color)))
  
  logging::logdebug(paste0("rayt_color=", rayt_color))
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

  logging::logdebug("Plotting REPIN population size.")
  p = facet_plot(p,
                panel='REPIN population size',
                data=popSize,
                geom=geom_segment,
                aes(x=0,
                    xend=repins,
                    y=y,
                    yend=y),
                size=bs,
                color=rayt_color
  )

  logging::logdebug("Adding theme.")
  p = p + theme_tree2()

  logging::logdebug("Adding labels.")
  p = facet_labeller(p,
                    c(Tree="",
                      RAYTs="Number of RAYTs")
  )

  logging::logdebug("Customizing theme.")
  p = p +
        theme(strip.text.x=element_text(hjust=0),
              strip.background = element_rect(color="black",
                                              fill="white",
                                              linetype="blank"
              )
        )

  logging::logdebug("Setting font size.")
  p = p +
          theme(text=element_text(size=fontsize)) +
          themeCurr

  logging::logdebug("Done, returning from function 'plotREPIN'")
  return(p)
}

#folder: folder that contains presAbs_* files
#type: REPIN type that is supposed to be viewed (_*, * is the type)
determineColor=function(associationFile){
  ass=read.delim(associationFile,header=TRUE)
  logging::logdebug(str(ass))
  logging::logdebug(ass)

  colorAss=c()
  for(i in 1:length(ass[,1])){
     rayt=paste0(ass[i,1],"_",ass[i,2])
     c="NA"
     if(nchar(ass[i,3])>0){
        split0=str_split(ass[i,3],",")
        c=colors[as.integer(split0[[1]][1])+1]
     }
	 else{
        c="grey"
     }
     temp=data.frame(repRAYT=rayt,color=c)
     colorAss=rbind(colorAss,temp)

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

drawRAYTphylogeny=function(data_dir){

  logging::loginfo("Drawing RAYT phylogeny.")

  raytseqFile=paste0(data_dir,"/repin_rayt_association.txt.fas")
  raytseqs=readDNAStringSet(raytseqFile,format="fasta")
  aln=muscle(raytseqs)
  raytAlnFile=paste0(data_dir,"/raytAln.phy")
  write.phylip(aln,raytAlnFile)
  system(paste0("phyml -i ",raytAlnFile," -m GTR"))
  raytTreeFile=paste0(raytAlnFile,"_phyml_tree.txt" )
  nwk=read.tree(raytTreeFile)
  logging::logdebug(nwk)
  colorDF = determineColor(paste0(data_dir,"/repin_rayt_association.txt"))
  logging::logdebug(colnames(colorDF))
  logging::logdebug(colorDF)
  onlyRAYTs <- colorDF[colorDF[,1] %in% nwk$tip.label, ]
  logging::logdebug(colnames(onlyRAYTs))
  logging::logdebug(onlyRAYTs)
  p <- ggtree(nwk)
  logging::logdebug(nwk)
  p <- p %<+% onlyRAYTs + geom_tiplab(aes(color=color))
  logging::logdebug("Added color tips")
  cols <- onlyRAYTs$color
  names(cols) <- onlyRAYTs$color
  logging::logdebug(cols)
  p <-  p + scale_color_manual(values=cols,guide=FALSE)
  logging::logdebug("Added color scale.")
  return(p)
}

logging::basicConfig()
logging::addHandler(writeToConsole)
# logging::addHandler(writeToFile, file="/tmp/shiny.log", level='DEBUG')
logging::logwarn("Starting up")

logging::setLevel(10) # DEBUG
