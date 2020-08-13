package statistics;

import java.io.*;
import java.util.*;


public class ProduceAndCompareDistribution {
	public static void main(String args[]){
		int rep=Integer.parseInt(args[0]);
		String replicon=args[3];
		ArrayList<Double> original=new ArrayList<Double>();
		if(replicon.equalsIgnoreCase("word")){
			original=AnalyseWordFreqFlank.readWF(new File(args[1]));
		}else{
			original=RepliconStats.readRepliconData(new File(args[1]));
		}
		double min=Double.parseDouble(args[2]);
		System.out.println(produceAndCompare(original, rep, min));
	}
	
	//given a distribution <rep> means are produced by randomly choosing n (length of distribution) members
	//returns how frequently (as a proportion of rep) the mean exceeds the minimum (double min) of the distribution that is to be compared
	public static double produceAndCompare(ArrayList<Double> original,int rep,double min){
		int exceed=0;
		for(int i=0;i<rep;i++){
			double sum=0;
			for(int j=0;j<original.size();j++){
				int rand=(int)(Math.random()*original.size());
				sum+=original.get(rand);
			}
			if((sum/original.size())<min){
				exceed++;
			}
		}
		return (exceed*1.0)/rep;
	}
	
	public static double getSignificance(int reps,File f1,File f2,boolean word){
		ArrayList<Double> original=new ArrayList<Double>();
		double min;
		if(word){
			original=AnalyseWordFreqFlank.readWF(f1);
			min=(new Stats(AnalyseWordFreqFlank.readWF(f2))).getAverage();
		}else{
			original=RepliconStats.readRepliconData(f1);
			min=(new Stats(RepliconStats.readRepliconData(f2))).getAverage();

		}
		return produceAndCompare(original, reps, min);

	}
	
}
