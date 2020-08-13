package frequencies.util;

import java.io.*;
import java.util.*;


public class CalculateNetworkProperties {
	NetworkProperties nw;
	int numSequences;
	HashMap<String,HashMap<String,Integer>> simNetwork;
	public ArrayList<String> hubs=new ArrayList<String>();
	public ArrayList<String> singlets=new ArrayList<String>();
	double hubThreshold=0.05;
		
	 HashMap<String,Integer> nodes=new HashMap<String,Integer>();
		File Rscript=new File("Rscript");

	public CalculateNetworkProperties(HashMap<String,HashMap<String,Integer>> simNetwork,HashMap<String,Integer> nodes,String genomeID,File outFolder,String regime){
		 this.simNetwork=simNetwork;
		 this.nodes=nodes;
		 nw=new NetworkProperties(regime,genomeID,outFolder);
		 nw.setNumNodes(nodes.size());
		 calculateDegreeDistribution();
		 calculateNumSequences();
		 calculateHubs();
		 calculateSinglets();
		 calculateAvgHubDistance();


	}

	private void calculateDegreeDistribution(){
		int[] degreeDistribution=init(simNetwork.size());
		String[] keys=simNetwork.keySet().toArray(new String[0]);
		for(int i=0;i<keys.length;i++){
			int degree=0;
			String keys2[]=simNetwork.get(keys[i]).keySet().toArray(new String[0]);
			for(int j=0;j<keys2.length;j++){
				if(simNetwork.get(keys[i]).get(keys2[j])==1){
					degree++;
				}
			}
			degreeDistribution[degree]++;
		}
		nw.setDegreeDistribution(degreeDistribution);
	}
	private int[] init(int size){
		int[] array=new int[size];
		for(int i=0;i<size;i++){
			array[i]=0;
		}
		return array;
	}



	//HARDER
	private void calculateAvgHubDistance(){
		double sum=0;
		int count=0;
		for(int i=0;i<hubs.size();i++){
			for(int j=i+1;j<hubs.size();j++){
				sum+=simNetwork.get(hubs.get(i)).get(hubs.get(j));
			    count++;
			}
		}
		nw.setavgHubDist((sum*1.0)/count);
	}
    
	public NetworkProperties getNetworkProperties(){
		return nw;
	}
	
	private void calculateNumSequences(){
		int sum=0;
		String keys[]=nodes.keySet().toArray(new String[0]);
		for(int i=0;i<keys.length;i++){
			sum+=nodes.get(keys[i]);
		}
		numSequences=sum;
	}
	
	private void calculateSinglets(){
		String keys[]=nodes.keySet().toArray(new String[0]);
		for(int i=0;i<keys.length;i++){
			int sum=0;
			for(int j=0;j<keys.length;j++){
				if(j==i)continue;
				if(simNetwork.get(keys[i]).get(keys[j])==1){
					sum++;
				}
			}
			if(sum<1){
				singlets.add(keys[i]);
			}
		}
		nw.setNumSinglets(singlets.size());
	}
	
	private void calculateHubs(){
		String keys[]=nodes.keySet().toArray(new String[0]);
		for(int i=0;i<keys.length;i++){
			double proportion=(nodes.get(keys[i])*1.0)/numSequences;
			if(proportion>hubThreshold){
				hubs.add(keys[i]);
			}
		}
		nw.setNumHubs(hubs.size());
	}
	

	

	



}
