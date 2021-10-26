#!/usr/bin/python3
import sys
outFolder=sys.argv[1]
id=sys.argv[2]
id=id.split("\\.")[0]
#check for RAYT presence
file=open(outFolder+"/repin_rayt_association.txt.fas","r")


def countLines(fname):
  try:
    file=open(fname)
    lines=file.readlines()
    return(len(lines))
  except IOError as e:
    return (0)
rayts=(countLines(outFolder+"/repin_rayt_association.txt.fas")/2)
overrep=countLines(outFolder+"/"+id+".overrep")
print("RAYTs\t"+str(rayts));
print("NumberOverRep\t"+str(overrep))
for i in range(0,5):
  lc=countLines(outFolder+"/"+id+"_"+str(i)+"/"+id+"_"+str(i)+"_largestCluster.nodes")
  print("largestCluster_"+str(i)+"\t"+str(lc))
