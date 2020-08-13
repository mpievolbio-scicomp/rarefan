package simulation;

import java.io.File;
import java.util.*;

import pairwiseAlignment.*;
import util.*;

public class CalculateTraits {
	
	ArrayList<Fasta> alignment;
	Double pwDist=null;
	Double nnDist=null;
	Double[][] pairwiseDists;
	
	public static void main(String args[]){
		File in=new File(args[0]);
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		CalculateTraits ct=new CalculateTraits(fas);
		System.out.println("NaN\tNaN\t"+ct.getAvgPairwiseDist()+"\t"+ct.getAvgNearestNeighborDist()+"\t"+in+"\t"+fas.size()+"\n");
	}
	
	public CalculateTraits(ArrayList<Fasta> alignment){
		this.alignment=alignment;
		pairwiseDists=calculatePairwiseDistances(alignment);
		calculateAverageNearestNeighborDistance();
		calculateAveragePairwiseDistance();
	}
	
	public static  Double[][] calculatePairwiseDistances(ArrayList<Fasta> alignment){
		Double[][] pairwiseDists=new Double[alignment.size()][alignment.size()];
		for(int i=0;i<alignment.size();i++){
			for(int j=i+1;j<alignment.size();j++){
				String a=alignment.get(i).getSequence();
				String b=alignment.get(j).getSequence();
			    double pw=NeedlemanWunsch.getPairwiseIdentity(a, b);
			    pairwiseDists[i][j]=pw;
			    pairwiseDists[j][i]=pw;
			}
			

		}
		return pairwiseDists;
	}
	
	public Double getAvgNearestNeighborDist(){
		return nnDist;
	}
	
	private void calculateAverageNearestNeighborDistance(){
		if(pairwiseDists==null)pairwiseDists=calculatePairwiseDistances(alignment);
		double sum=0;
		int count=0;
		for(int i=0;i<pairwiseDists.length;i++){
			double max=Double.MIN_VALUE;
			for(int j=0;j<pairwiseDists[i].length;j++){
				if(j==i)continue;
				double pw=pairwiseDists[i][j];
				if(pw>max){
					max=pw;
				}
			}
			sum+=max;
			count++;
		}
		nnDist=(sum)/count;
	}

	public double getAvgPairwiseDist(){
		return pwDist;
	}
	
	private void calculateAveragePairwiseDistance(){
		if(pairwiseDists==null)pairwiseDists=calculatePairwiseDistances(alignment);
		double sum=0;
		int count=0;
		for(int i=0;i<pairwiseDists.length;i++){
			for(int j=i+1;j<pairwiseDists[i].length;j++){
				sum+=pairwiseDists[i][j];
				count++;
			}
		}
		pwDist=(sum)/count;
	}
	
	
}
