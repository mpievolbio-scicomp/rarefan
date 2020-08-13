package sequenceDistribution;

import java.io.*;
import java.util.*;


import util.*;

public class DistributionRelativeToGenes {
	public static void main(String args[]){
		File positions=new File(args[0]);
		int steps=Integer.parseInt(args[1]);
		String genome=Fasta.readFasta(new File(args[2])).get(0).getSequence();
		int genomeSize=genome.length();
		File outFolder=new File(args[3]);
		String searchInput=args[4];
		ReadSearchOutput rs;
		if(searchInput.equalsIgnoreCase("soap")){
			rs=new ReadSoap(positions);
		}else{
			ArrayList<Integer> pos=Find.getPositions(args[0].toLowerCase(),genome.toLowerCase());
			pos.addAll(Find.getPositions(DNAmanipulations.reverse(args[0]).toUpperCase(), genome.toUpperCase()));
			rs=new ReadBitSetSearch(pos,args[0].length());
		}
		print(getClosestRegion(rs,steps,genomeSize),outFolder);
	}
	
	public static void print(ComplexType ct, File outFolder){
		try{
			
			for(int i=0;i<ct.avg.size();i++){
				System.out.println(ct.query.get(i)+" "+ct.avg.get(i)+" "+ct.region.get(i));
				BufferedWriter bw=new BufferedWriter(new FileWriter(outFolder+"/"+ct.query.get(i)+".out"));
				ArrayList<Integer> avgs=ct.all.get(ct.query.get(i));
				for(int j=0;j<avgs.size();j++){
					bw.write(avgs.get(j)+"\n");
				}
				bw.close();
			}
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static ComplexType getClosestRegion(ReadSearchOutput rs,int steps,int genomeSize){
		ArrayList<Integer> minAvg=new ArrayList<Integer>();
		ArrayList<String> queries=new ArrayList<String>();
		HashMap<String,ArrayList<Integer>> positions=new HashMap<String,ArrayList<Integer>>();
		HashMap<String,ArrayList<Integer>> avgs=new HashMap<String, ArrayList<Integer>>();
		ArrayList<Integer> region=new ArrayList<Integer>();
		for(int i=0;i<rs.getStart().size();i++){
			if(positions.containsKey(rs.getQuery().get(i))){
				positions.get(rs.getQuery().get(i)).add(rs.getStart().get(i));
			}else{
				queries.add(rs.getQuery().get(i));
				ArrayList<Integer> temp=new ArrayList<Integer>();
				temp.add(rs.getStart().get(i));
				positions.put(rs.getQuery().get(i),temp);
			}
		}
		for(int i=0;i<queries.size();i++){
			minAvg.add(Integer.MAX_VALUE);
			region.add(-1);
			avgs.put(queries.get(i), new ArrayList<Integer>());
		}
		for(int j=0;j<genomeSize;j+=steps){
			for(int k=0;k<queries.size();k++){
				ArrayList<Integer> pos=positions.get(queries.get(k));
				int sum=0;
				for(int i=0;i<pos.size();i++){
					sum+=Math.min(Math.abs(pos.get(i)-j), Math.min(genomeSize-pos.get(i),genomeSize-j)+ Math.min(pos.get(i),j));
				}
				avgs.get(queries.get(k)).add(sum/pos.size());
				if(sum/pos.size()<minAvg.get(k)){
					minAvg.set(k,sum/pos.size());
					region.set(k,j);
				}
				

			}
		}
		return new ComplexType(minAvg,region,queries,avgs);
	}
	public static class ComplexType{
		public ComplexType(ArrayList<Integer> Avg,ArrayList<Integer> Region,ArrayList<String> Query,HashMap<String, ArrayList<Integer>> All){
			all=All;
			avg=Avg;
			region=Region;
			query=Query;
		}
		HashMap<String,ArrayList<Integer>> all=new HashMap<String, ArrayList<Integer>>();
		ArrayList<Integer> avg=new ArrayList<Integer>();
		ArrayList<Integer> region=new ArrayList<Integer>();
		ArrayList<String> query=new ArrayList<String>();
		
	}
}
