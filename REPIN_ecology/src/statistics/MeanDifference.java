package statistics;

import java.io.*;
import java.util.*;

public class MeanDifference {
	public static void main(String args[]){
		File sample1=new File(args[0]);
		File sample2=new File(args[1]);
		int repeats=Integer.parseInt(args[2]);
		ArrayList<Double> s1=readSample(sample1);
		ArrayList<Double> s2=readSample(sample2);
		ArrayList<Double> mix=new ArrayList<Double>();
		mix.addAll(s2);
		mix.addAll(s1);
		Stats stats=new Stats(s1);
		double mean1=stats.getAverage();
		System.out.println(mean1+" "+stats.getStandardDeviation());
		stats=new Stats(s2);
		double mean2=stats.getAverage();
		System.out.println(mean2+" "+stats.getStandardDeviation());
		double meanDiff=Math.abs(mean1-mean2);
		int count=0;
		for(int i=0;i<repeats;i++){
			ArrayList<Double> subset1=chooseSubset(mix,s1.size());
			ArrayList<Double> subset2=chooseSubset(mix,s2.size());
			stats=new Stats(subset1);
			double meanss1=stats.getAverage();
			stats=new Stats(subset2);
			double meanss2=stats.getAverage();
			double meanDiffss=Math.abs(meanss1-meanss2);
			if(meanDiffss>meanDiff){
				count++;
			}
			if(i%100==0){
				System.out.println(i);
			}
		}
		System.out.println("P-value of hypothesis (two sets originating from the same distribution) being true: "+(count/(repeats*1.0)));
	}
	
	private static ArrayList<Double> chooseSubset(ArrayList<Double> set,int number){
		ArrayList<Double> subset=new ArrayList<Double>();
		
		for(int i=0;i<number;i++){
			subset.add(set.get((int)(set.size()*Math.random())));
		}
		
		return subset;
		
	}
	
	public static ArrayList<Double> readSample(File in){
		ArrayList<Double> samples=new ArrayList<Double>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				samples.add(Double.parseDouble(split[0]));
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return samples;
	}
	
}
