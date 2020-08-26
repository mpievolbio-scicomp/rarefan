library(ggpmisc)
library(ggtree)
library(ape)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
theme=theme(axis.line.x = element_line(colour = "black"),legend.key = element_rect(fill = "white"),axis.line.y = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),legend.justification = c(0, 1), legend.position = c(0.5, 1),legend.title=element_blank(),legend.text = element_text(hjust=0),panel.spacing=unit(2,"lines"))
folderNei="/Users/bertels/Documents/Arne/mutationFrequencies/Neisseria/"
fs=10
#folder: folder that contains presAbs_* files
#treeFile: name of newick tree file 
#type: REPIN type that is supposed to be viewed (_*, * is the type)
#colorbars: color for REPIN population size bars in hex e.g. "40e0d0"
#bs: bar size e.g. 2
#fs: font size e.g. 10
plotREPINs=function(folder,treeFile,type,colorBars,bs,fs){
  #rayts=read.table(paste0(folder,"REPINs21/rayts.txt"),colClasses=c(rep("character",3)))
  themeCurr=theme(axis.line.x = element_line(colour = "black"),legend.key = element_rect(fill = "white"),axis.line.y = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),legend.position="none",panel.spacing=unit(2,"lines"),plot.margin = unit(c(5.5,12,5.5,10.5), "pt"))
  t=read.table(paste0(folder,"presAbs_",type,".txt"),sep="\t")
  popSize=data.frame(name=t[,1],rayts=t[,2],repins=t[,3],prop=(t[,5]/t[,3]),propAll=(t[,3]/t[,6]),numClus=t[,7],diffRAYTCluster=t[,7]-t[,2])


 tree=read.tree(paste0(folder,treeFile))
  #tree$tip.label=gsub("\\.","_",tree$tip.label)
  tips=tree$tip.label
  #gons=tips[grep(rootName,tips)]
  #rest=tips[!((1:length(tips))%in%grep(rootName,tips))]
  #gonClade=getMRCA(tree,gons)
  #tree=root(tree,node=gonClade,resolve.root=T)

  #printCorrelations(popSize)
  #printCorrelations(popSize[popSize$name%in%rest,])
  #printCorrelations(popSize[!(popSize$name%in%rest),])

  #gonClade=getMRCA(tree,gons)
  #restClade=getMRCA(tree,rest)

  p=ggtree(tree)+scale_x_continuous(breaks=scales::pretty_breaks(n=3))+xlim_tree(0.045)
  p2=facet_plot(p,panel='RAYTs',data=popSize,geom=geom_segment,aes(x=0,xend=rayts,y=y,yend=y),size=bs,color=colorBars)
  p3=facet_plot(p2,panel='Number of \nREPIN clusters',data=popSize,geom=geom_segment,aes(x=0,xend=numClus,y=y,yend=y),size=bs,color=colorBars)
  p4=facet_plot(p3,panel='Difference \nRAYTs-REPIN \nclusters',data=popSize,geom=geom_segment,aes(x=0,xend=diffRAYTCluster,y=y,yend=y),size=bs,color=colorBars)
  p5=facet_plot(p4,panel='REPIN\npopulation\nsize',data=popSize,geom=geom_segment,aes(x=0,xend=repins,y=y,yend=y),size=bs,color=colorBars)
  p6=p5+theme_tree2()
  p7=facet_labeller(p6,c(Tree="",RAYTs="Number of \nRAYTs"))
  p8=p7+theme(strip.text.x=element_text(hjust=0),
   strip.background = element_rect(
     color="black", fill="white", linetype="blank"
     )
   )#+geom_tiplab()
  #p7=p6+theme_tree2()+geom_cladelabel(node=29,label="",offset=0.1,color="red")+geom_tiplab()
  p8=p8+theme(text=element_text(size=fs))+themeCurr
  return(p8)
}

#folder: folder that contains presAbs_* files
#type: REPIN type that is supposed to be viewed (_*, * is the type)


plotCorrelationSingle=function(folder,type,xlim,ylim,theme,fs,pvLabelX,pvLabelY,subsetSmooth=F,from=F,to=F,repinThreshold=0,name=F,labelOdd){
    t=read.table(paste0(folder,"presAbs_",type,".txt"),sep="\t")
    t$propMaster=t[,5]/t[,3]
    t$numRepin=t[,3]
    p=ggplot(t,aes(x=propMaster,y=numRepin))+geom_point()
    p=p+geom_smooth(method=lm, se=FALSE, fullrange=TRUE,formula=y~x,show.legend=F)
    p=p+stat_fit_glance(method = 'lm',method.args=list(formula=y~x),aes(label = paste("P-value = ", signif(..p.value.., digits = 2), sep = "")),label.x=pvLabelX,label.y=pvLabelY, size = fs/3)

    p=p+xlim(xlim)+ylim(ylim)+theme+xlab("Proportion master sequence (~Replication rate)")+ylab("REPIN population size")
    p=p+theme(axis.text=element_text(size=fs),text=element_text(size=fs))
    return(p)
}


pN=plotREPINs(folderNei,"neisseria.nwk",0,"#40e0d0",2,fs)
pN
pNeisseria=plotCorrelationSingle(folderNei,0,c(0,1),c(0,320),theme,fs,"left","bottom")
pNeisseria


