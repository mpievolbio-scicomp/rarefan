#!/usr/bin/env Rscript
# Required libraries
library(ggpmisc)
library(ggtree)
library(ape)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)

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

# Define plot routine
plotREPINs=function(folder,treeFile,type,colorBars,bs,fontsize){
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

  t=read.table(paste0(folder,"/presAbs_",type,".txt"),sep="\t")
#popSize=data.frame(
#  )
#


  popSize=data.frame(
                     name=t$strain,
                     rayts=t$numRAYTs,
                     repins=t$numREPINs,
                     prop=(t$mastersequenceFreq/t$numREPINs),
                     propAll=(t$numREPINs/t$"allREP\\REPINFreq"),
                     numClus=t$numberOfRepinClusters,
                     diffRAYTCluster=(t$numberOfRepinClusters-t$numRAYTs)
  )


  tree=read.tree(paste0(folder,"/",treeFile))
  tips=tree$tip.label

  p=ggtree(tree)+
      scale_x_continuous(breaks=scales::pretty_breaks(n=3))+
      xlim_tree(0.045)

  p2=facet_plot(p,
                panel='RAYTs',
                data=popSize,
                geom=geom_segment,
                aes(x=0,
                    xend=rayts,
                    y=y,
                    yend=y),
                size=bs,
                color=colorBars
  )

  p3=facet_plot(p2,
                panel='Number of \nREPIN clusters',
                data=popSize,
                geom=geom_segment,
                aes(x=0,
                    xend=numClus,
                    y=y,
                    yend=y),
                size=bs,
                color=colorBars
  )

  p4=facet_plot(p3,
                panel='Difference \nRAYTs-REPIN \nclusters',
                data=popSize,
                geom=geom_segment,
                aes(x=0,
                    xend=diffRAYTCluster,
                    y=y,
                    yend=y),
                size=bs,
                color=colorBars
  )

  p5=facet_plot(p4,
                panel='REPIN\npopulation\nsize',
                data=popSize,
                geom=geom_segment,
                aes(x=0,
                    xend=repins,
                    y=y,
                    yend=y),
                size=bs,
                color=colorBars
  )

  p6=p5+theme_tree2()

  p7=facet_labeller(p6,
                    c(Tree="",
                      RAYTs="Number of \nRAYTs")
  )

  p8 = p7+
        theme(strip.text.x=element_text(hjust=0),
              strip.background = element_rect(color="black",
                                              fill="white",
                                              linetype="blank"
              )
        )

  p8 = p8 +
          theme(text=element_text(size=fontsize)) +
          themeCurr

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
    t=read.table(paste0(folder,"/presAbs_",type,".txt"),sep="\t")
    t$propMaster=(t$mastersequenceFreq/t$numREPINs)
    t$numRepin=t$numREPINs
    p=ggplot(t,
             aes(x=propMaster,
                 y=numRepin))+
            geom_point()

    p=p+geom_smooth(method=lm,
                    se=FALSE,
                    fullrange=TRUE,
                    formula=y~x,
                    show.legend=F
    )

    p=p+stat_fit_glance(method = 'lm',
                        method.args=list(formula=y~x),
                        aes(label = paste("P-value = ",
                                          signif(..p.value..,
                                                 digits = 2),
                                          sep = "")),
                        label.x=pvLabelX,
                        label.y=pvLabelY,
                        size = fontsize/3)

    p=p+xlim(xlim)+
        ylim(ylim)+
        theme+
        xlab("Proportion master sequence (~Replication rate)")+
        ylab("REPIN population size")

    p=p+theme(axis.text=element_text(size=fontsize),text=element_text(size=fontsize))

    return(p)
}


repins_plot=plotREPINs(data_dir,treefile,0,"#40e0d0",2,fontsize)
ggsave(paste0(data_dir, '/', 'repins.png'), plot=repins_plot)

correlation_plot = plotCorrelationSingle(data_dir,0,c(0,1),c(0,320),theme,fontsize,"left","bottom")
ggsave(paste0(data_dir, '/', 'correlations.png'), plot=correlation_plot)

