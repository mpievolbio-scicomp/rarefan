#!/usr/bin/env Rscript
# Required libraries
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

logging::basicConfig()
#logging::setLevel(20) # INFO
logging::setLevel(10) # DEBUG
# NOTE: Requires phyml to be installed in system's PATH.

# Parse command line args
args = commandArgs(trailingOnly=TRUE)
#folder: folder that contains presAbs_* files
#treeFile: name of newick tree file 
#type: REPIN type that is supposed to be viewed (_*, * is the type)

max_number_of_expected_args = 3
min_number_of_expected_args = 1
if (length(args)<min_number_of_expected_args ) { 
    stop("Usage: Rscript displayREPINsAndRAYTs.R DIR [TYPE [TREEFILE TYPE]]", call.=FALSE)
}
if (length(args)>max_number_of_expected_args ) { 
    stop("Usage: Rscript displayREPINsAndRAYTs.R DIR [TYPE [TREEFILE TYPE]]", call.=FALSE)
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

# Set theme for all plotse
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

# Font size
fontsize=14

#determine color for all plots
#returns a table with two columns repintype/rayt and color
determineColor=function(associationFile){
  ass=read.delim(associationFile,header=TRUE)
  colors=c("blue","red","green","purple","teal","orange")
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
  for(i in groups){
     if(nchar(i)>0){
        split=str_split(i,",")
        for(j in split[[1]]){
           pos=as.integer(split[[1]][1])
           j=as.integer(j)
           temp=data.frame(repRAYT=j,color=colors[pos+1])
           colorAss=rbind(colorAss,temp)
        }
     }
  }
  print(colorAss)
  return(colorAss)
}

associationFile=paste0(data_dir,"/repin_rayt_association.txt")
logging::logdebug("Setting colors according to association file %s.", associationFile)
colorDF=determineColor(associationFile)

# Define plot routine
plotREPINs=function(folder,treeFile,type,colorBars,bs,fontsize){
	logging::logdebug("Setting theme.")
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
  
  data_file = paste0(folder,"/presAbs_",type,".txt")
  logging::logdebug("Reading table from %s.", data_file)
  t=tryCatch(read.table(data_file,sep="\t", skip=1),
		  error=function(e)	
		  logging::logwarn("File %s is empty.", data_file)
)
  
	if(typeof(t) == "logical") {
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
      xlim_tree(0.045)+
      geom_tiplab()

  logging::logdebug("Setting up facet plot p2.")
  p2=facet_plot(p,
                panel='RAYTs',
                data=association[association$repintype==type,],
                geom=geom_segment,
                aes(x=0,
                    xend=rayts,
                    y=y,
                    yend=y),
                size=bs#,
                #color=colorDF[colorDF$repRAYT==type,2]
  )

#  p3=facet_plot(p2,
#                panel='Number of \nREPIN clusters',
#                data=popSize,
#                geom=geom_segment,
#                aes(x=0,
#                    xend=numClus,
#                    y=y,
#                    yend=y),
#                size=bs,
#                color=colorBars
#  )

#  p4=facet_plot(p3,
#                panel='Difference \nRAYTs-REPIN \nclusters',
#                data=popSize,
#                geom=geom_segment,
#                aes(x=0,
#                    xend=diffRAYTCluster,
#                    y=y,
#                    yend=y),
#                size=bs,
#                color=colorBars
#  )

  logging::logdebug("Setting up facet plot p5.")
  p5=facet_plot(p2,
                panel='REPIN\npopulation\nsize',
                data=popSize,
                geom=geom_segment,
                aes(x=0,
                    xend=repins,
                    y=y,
                    yend=y),
                size=bs#,
                #color=colorDF[colorDF$repRAYT==type,2]
  )

  logging::logdebug("Setting up facet plot p6.")
  p6=p5+theme_tree2()

  logging::logdebug("Setting up facet plot p7.")
  p7=facet_labeller(p6,
                    c(Tree="",
                      RAYTs="Number of \nRAYTs")
  )

  logging::logdebug("Setting up facet plot p8.")
  p8 = p7+
        theme(strip.text.x=element_text(hjust=0),
              strip.background = element_rect(color="black",
                                              fill="white",
                                              linetype="blank"
              )
        )

  logging::logdebug("Finalizing facet plot p8.")
  p8 = p8 +
          theme(text=element_text(size=fontsize)) +
          themeCurr
  logging::logdebug("Done, returning from function 'plotREPIN'")
  return(p8)
}

#folder: folder that contains presAbs_* files
#type: REPIN type that is supposed to be viewed (_*, * is the type)


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
                               labelOdd){
						  
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
    cols[cols>0]=colorDF[colorDF$repRAYT==type,]$color
    cols[cols==0]="black"
	
	logging::logdebug("Setting up ggplots.")
    p=ggplot(t,
             aes(x=propMaster,
                 y=numRepin,col=factor(color)))+
            scale_color_manual(values=cols,guide=FALSE)+
            geom_point()

 #   p=p+geom_smooth(method=lm,
 #                   se=FALSE,
 #                   fullrange=TRUE,
 #                   formula=y~x,
 #                   show.legend=F
 #   )

 #   p=p+stat_fit_glance(method = 'lm',
 #                       method.args=list(formula=y~x),
 #                       aes(label = paste("P-value = ",
 #                                         signif(..p.value..,
 #                                                digits = 2),
 #                                         sep = "")),
 #                       label.x=pvLabelX,
 #                       label.y=pvLabelY,
 #                       size = fontsize/3)

	logging::logdebug("Adding limits, theme, and axis labels.")
    p=p+xlim(xlim)+
        ylim(ylim)+
        theme+
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
  onlyRAYTs=colorDF[colorDF[,1]%in%nwk$tip.label,]
  p=ggtree(nwk)
  p=p%<+%onlyRAYTs+geom_tiplab(aes(color=color))
  cols=onlyRAYTs$color
  names(cols)=onlyRAYTs$color
  p = p + scale_color_manual(values=cols,guide=FALSE)

  logging::loginfo("Saving RAYT phylogeny plots.")
  ggsave(paste0(data_dir,"/raytTree.png"))
  
}
drawRAYTphylogeny(data_dir)

for(i in 0:5){
	logging::loginfo("Plotting REPINS [i=%d]", i)
    repins_plot=plotREPINs(data_dir,treefile,i,"#40e0d0",2,fontsize)
	logging::logdebug("Plotting done, saving.")
    ggsave(paste0(data_dir, '/', 'repins_',i,'.png'), plot=repins_plot)
	logging::logdebug("Saving done.")

	logging::loginfo("Plotting correlations [i=%d]", i)
    
	correlation_plot = plotCorrelationSingle(data_dir,i,c(0,1),c(0,320),theme,fontsize,"left","bottom")
    ggsave(paste0(data_dir, '/', 'correlations_',i,'.png'), plot=correlation_plot)
}
