package frequencies.util;

import java.io.*;
import java.util.*;

import frequencies.RAYTfrequencyGraph;
import util.*;

public class SimilarityNetwork {
	public static void main(String args[]){
		File fasta=new File(args[0]);
		File out=new File(fasta+".nw");
		SimilarityNetwork csn=new SimilarityNetwork(fasta);
		csn.writeCytoscapeInput(out);
	}
	 HashMap<String,Integer> nodes=new HashMap<String,Integer>();
	File Rscript=new File("Rscript");
	 HashMap<String,HashMap<String,Integer>> simNetwork=new HashMap<String,HashMap<String,Integer>>();
	ArrayList<Fasta> fas;
	public SimilarityNetwork(File fasta){
		fas=Fasta.readFasta(fasta);
		 calculateSimilarityNetwork();
	}
	

	
	public HashMap<String,HashMap<String,Integer>> getNetwork(){
		return simNetwork;
	}
	
	
	
	public HashMap<String,Integer> getNodes(){
		return nodes;
	}
	
	public void writeCytoscapeInput(File out){

        System.out.println("Writing Cytoscape input to " + out + ".");
		try{
			if(simNetwork.size()>1){
				BufferedWriter bw=new BufferedWriter(new FileWriter(out));
				String keys1[]=simNetwork.keySet().toArray(new String[0]);
				for(int i=0;i<keys1.length;i++){
					String[] keys2=simNetwork.get(keys1[i]).keySet().toArray(new String[0]);
					for(int j=0;j<keys2.length;j++){
						int diff=simNetwork.get(keys1[i]).get(keys2[j]);
						if(diff==1){
							bw.write(keys1[i]+"\t"+keys2[j]+"\t"+diff+"\n");
						}
					}
				}
				bw.close();
			}else{
				out.delete();
			}
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public void writeNodes(File out){
		try{
			if(nodes.size()>1){
				BufferedWriter bw=new BufferedWriter(new FileWriter(out));
				String[] keys=nodes.keySet().toArray(new String[0]);
				for(int i=0;i<keys.length;i++){
					int freq=nodes.get(keys[i]);
					bw.write(keys[i]+"\t"+freq+"\n");
				}
				bw.close();
			}else{
				out.delete();
			}
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private void calculateSimilarityNetwork(){
		nodes=new HashMap<String,Integer>();
		for(int i=0;i<fas.size();i++){
			String a=fas.get(i).getSequence();
			String namea=a+"_"+RAYTfrequencyGraph.getFreq(fas.get(i).getIdent());
			simNetwork.put(namea, new HashMap<String,Integer>());
			if(!nodes.containsKey(namea))nodes.put(namea,Integer.parseInt(namea.split("_")[1]));
			for(int j=0;j<fas.size();j++){
				if(j==i)continue;
				String b=fas.get(j).getSequence();
				int diff=getNumDifferences(a, b);
				String nameb=b+"_"+RAYTfrequencyGraph.getFreq(fas.get(j).getIdent());
				simNetwork.get(namea).put(nameb, diff);
				if(!nodes.containsKey(nameb))nodes.put(nameb,Integer.parseInt(nameb.split("_")[1]));
			}
		}
	}
	
	public static int getNumDifferences(String a,String b){
		int differences=0;
		for(int i=0;i<a.length();i++){
			if(a.charAt(i)!=b.charAt(i)){
				differences++;
			}
		}
		return differences;
	}
	
}
