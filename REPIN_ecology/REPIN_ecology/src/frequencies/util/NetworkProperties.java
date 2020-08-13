package frequencies.util;

import java.io.*;
import java.util.*;


public class NetworkProperties {
	static File Rscript=new File("Rscript");
	HashMap<String,Double> propertyHash=new HashMap<String,Double>();
	static String[] properties=new String[]{"numSinglets","avgHubDist","numNodes","numHubs"};
	int[] degreeDistribution;
	File outFolder;
	String genomeID;
	String regime;
	public NetworkProperties(String regime,String genomeID,File outFolder){
		this.regime=regime;
		this.outFolder=outFolder;
		this.genomeID=genomeID;
	}
	public String getRegime(){
		return regime;
	}
	
	public File getOutFolder(){
		return outFolder;
	}
	
	public String getGenomeID(){
		return genomeID;
	}
	
	public void setDegreeDistribution(int[] degreeDist){
		degreeDistribution=degreeDist;
	}

	public void setNumSinglets(int numSinglets){
		propertyHash.put(properties[0], numSinglets*1.0);
	}
	
	public void setavgHubDist(double avgHubDist){
		propertyHash.put(properties[1], avgHubDist*1.0);
	}


	public void setNumNodes(int numNodes){
		propertyHash.put(properties[2], numNodes*1.0);
	}

	public void setNumHubs(int numHubs){
		propertyHash.put(properties[3], numHubs*1.0);
	}
	public ArrayList<Double> getAll(){
		ArrayList<Double> all=new ArrayList<Double>();
		for(int i=0;i<properties.length;i++){
			all.add(propertyHash.get(properties[i]));
		}
		return all;
	}
	
	public void setAll(ArrayList<Double> setproperties){
		for(int i=0;i<setproperties.size();i++){
			propertyHash.put(properties[i],setproperties.get(i));
		}
	}
	public int[] getDegreeDistribution(){
		return degreeDistribution;
	}
	public void writeDegreeDist(File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<degreeDistribution.length;i++){
				bw.write(i+"\t"+degreeDistribution[i]+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}

	}

//	public void plotDegreeDistribution(File out){
//		RCode rc=R_functions.plot_InitPdf(out, 12, 12);
//
//		rc.addIntArray("degreeDist", degreeDistribution);
//		rc.addRCode("hist(degreeDist)");
//		R_functions.runRCode(rc, Rscript);
//		R_functions.writeRCode(rc, new File(out+".R"));
//	}

	private Double get(String property){
		return propertyHash.get(property);
	}
	
	public static void printProperties(ArrayList<NetworkProperties> list,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<list.size();i++){
				String regime=list.get(i).getRegime();
				String genomeID=list.get(i).getGenomeID();
				String heading=regime==null?genomeID:regime;
				bw.write("\t"+heading);
			}
			bw.write("\n");
			for(int i=0;i<properties.length;i++){
				bw.write(properties[i]);
				for(int j=0;j<list.size();j++){
					bw.write("\t"+list.get(j).get(properties[i]));
				}
				bw.write("\n");
			}
			bw.write("goodnessOfFit");
			ArrayList<Double> scores=calculateScore(list);
			for(int i=0;i<scores.size();i++){
				if(scores.get(i)==Double.NaN){
					bw.write("\t-");
				}else{
					bw.write("\t"+scores.get(i));
				}
			}
			bw.write("\n");
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private static ArrayList<Double> calculateScore(ArrayList<NetworkProperties> list){
		ArrayList<Double> scores=new ArrayList<Double>();
		NetworkProperties real=null;
		NetworkProperties avg=null;
		NetworkProperties se=null;
		int index=0;
		for(int k=0;k<list.size();k++){
			NetworkProperties nw=list.get(k);
			if(nw.regime==null){
				real=list.get(k);
				scores.add(Double.NaN);
				index=k;
			}else{
				if((k-index)%2==1){
					avg=list.get(k);

					scores.add(Double.NaN);
				}else if((k-index)%2==0){
					se=list.get(k);

					if(avg!=null){
						scores.add(getScore(real,avg,se));
					}else{
						scores.add(Double.NaN);
					}
				}
			}

		}
		return scores;
	}
	
	private static Double getScore(NetworkProperties real,NetworkProperties sim,NetworkProperties se){				
		int count=0;
		int fit=0;
		for(int i=0;i<properties.length;i++){

			if(real.get(properties[i]).equals(Double.NaN))continue;
		
			if(sim.get(properties[i]).equals(Double.NaN))continue;
			count++;
			fit+=fits(real.get(properties[i]),sim.get(properties[i]),se.get(properties[i]))?1:0;
		}
		return fit*1.0/count;
	}

	
		
	private static boolean fits(double real,double sim,double se){
		if(real<=sim+2*se&&real>=sim-2*se){
			return true;
		}else{
			return false;
		}
	}
	
}
