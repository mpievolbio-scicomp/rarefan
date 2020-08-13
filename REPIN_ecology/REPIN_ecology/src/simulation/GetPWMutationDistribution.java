package simulation;

import java.io.*;
import java.util.*;

import util.*;

public class GetPWMutationDistribution {
	
	class NNIndexPairs{
		int index1;
		int index2;
		public NNIndexPairs(int i1,int i2){
			index1=i1;
			index2=i2;
		}
	}
	File Rscript=new File("/usr/local/bin/Rscript");
	public static void main(String args[]){
		File alg=new File(args[0]);
		File[] files=alg.listFiles();
		for(int i=0;i<files.length;i++){
		    if(files[i].getName().endsWith("fas")){
		    	GetPWMutationDistribution gpwmd=new GetPWMutationDistribution(files[i]);
		    }
		}
	}
	
	ArrayList<Fasta> alg;
	ArrayList<NNIndexPairs> indexPairs=new ArrayList<NNIndexPairs>();
	ArrayList<Integer> mutationPositions=new ArrayList<Integer>();
	public GetPWMutationDistribution(File alg){
		this.alg=Fasta.readFasta(alg);
		String split[]=alg.getName().split("\\.|\\/");
		String name=split[split.length-2];
		calculateNearestNeighborIndexPairs();
		calculateMutationPositions();
		//plotHistogram(alg.getParentFile(),name);
	}
	
	private void calculateMutationPositions(){
		for(int i=0;i<indexPairs.size();i++){
			int i1=indexPairs.get(i).index1;
			int i2=indexPairs.get(i).index2;
			mutationPositions.addAll(getMutationPositions(alg.get(i1).getSequence(),alg.get(i2).getSequence()));
		}
	}
	
	private ArrayList<Integer> getMutationPositions(String s1,String s2){
		ArrayList<Integer> pos=new ArrayList<Integer>();
		for(int i=0;i<s1.length();i++){
			char c1=s1.charAt(i);
			char c2=s2.charAt(i);
			if(c1!=c2)pos.add(i);
			
		}
		return pos;
	}
	
//	private void plotHistogram(File outFolder,String name){
//		File out=new File(outFolder+"/"+name+".pdf");
//		RCode rc=R_functions.plot_InitPdf(out, 10, 10);
//		R_functions.makeHistogram(rc, mutationPositions, name);
//		R_functions.runRCode(rc, Rscript);
//		R_functions.writeRCode(rc, new File(out+".R"));
//	}
	
	private void calculateNearestNeighborIndexPairs(){
		Double[][] pairwiseDists=CalculateTraits.calculatePairwiseDistances(alg);
		
		for(int i=0;i<pairwiseDists.length;i++){
			double max=Double.MIN_VALUE;
			int maxIndex=-1;
			for(int j=0;j<pairwiseDists[i].length;j++){
				if(j==i)continue;
				double pw=pairwiseDists[i][j];
				if(pw>max){
					max=pw;
					maxIndex=j;
				}
			}
			indexPairs.add(new NNIndexPairs(i,maxIndex));
			
		}
	}
	
}
