package util.phylogenetics;

import java.io.File;
import java.util.*;
import java.util.AbstractMap.SimpleEntry;

import util.*;
/**
 * Creates a UPGMA tree
 * @author bertels
 *
 */
public class UPGMA {
	public class MergeResult{
		public MergeResult(String k1k2,double distK1,double distK2){
			this.k1k2=k1k2;
			this.distK1=distK1;
			this.distK2=distK2;
		}
		String k1k2;
		double distK1;
		double distK2;
	}
	public class UPGMAEntry extends SimpleEntry<String,String> implements Comparable<UPGMAEntry>{

		
		private static final long serialVersionUID = 1L;

		Double distance;
		

		
		public UPGMAEntry(String key,String value,double distance) {
			super(key,value);
			this.distance=distance;
			
		}


		
		public int compareTo(UPGMAEntry o) {

			return this.distance.compareTo(o.distance);
		}
		
		
	}
	
	public static void main(String args[]){
		File alg=new File(args[0]);
		Alignment align=new Alignment(Fasta.readFasta(alg));
		UPGMA upgma=new UPGMA(align);
		System.out.println(upgma.runUPGMA());
	}
	
	HashMap<String,HashMap<String,Double>> distMat=new HashMap<String, HashMap<String,Double>>();
	Alignment alg;
	HashMap<String,Double> lengthMap=new HashMap<String, Double>();
	
	public UPGMA(Alignment alg){
		this.alg=alg;
		calculateDistanceMatrix();
		printDistMat();
	}
	
	public String runUPGMA(){
		while(distMat.size()>1){
			UPGMAEntry ua=getLowest();
			//System.out.println(ua.getKey()+" "+ua.getValue()+" "+ua.distance);
			mergeEntries(ua.getKey(), ua.getValue());
			System.out.println("\n");
			printDistMat();
		}
		return distMat.keySet().toArray(new String[0])[0];
	}
	
	public double getDistance(String k1,String k2){
		return distMat.get(k1).get(k2);
	}

	public int getHammingDistance(String target,String sequence){
		int differences=0;
		for(int i=0;i<target.length();i++){
			if(target.charAt(i)!=sequence.charAt(i)){
				differences++;
			}
		}
		return differences;
	}
	
	private UPGMAEntry getLowest(){
		String[] klist=distMat.keySet().toArray(new String[0]);
		int minI=0;
		int minJ=0;
		double min=Double.POSITIVE_INFINITY;
		for(int i=0;i<klist.length;i++){
			for(int j=i+1;j<klist.length;j++){
				double dist=distMat.get(klist[i]).get(klist[j]);
				if(min>dist){
					min=dist;
					minI=i;
					minJ=j;
				}
			}
		}
		return new UPGMAEntry(klist[minI], klist[minJ], min);
	}
	
	private void delete(String key){
		String[] klist=distMat.keySet().toArray(new String[0]);
		distMat.remove(key);
		for(int i=0;i<klist.length;i++){
			if(!key.equals(klist[i])){
				distMat.get(klist[i]).remove(key);
			}
		}
	}
	
	public MergeResult /*merged name*/ mergeEntries(String key1,String key2){
		String[] klist=distMat.keySet().toArray(new String[0]);
		double k1k2dist=distMat.get(key1).get(key2)/2;
		double k1l=lengthMap.containsKey(key1)?lengthMap.get(key1):0;
		double k2l=lengthMap.containsKey(key2)?lengthMap.get(key2):0;

		String k1k2="("+key1+":"+(k1k2dist-k1l)+","+key2+":"+(k1k2dist-k2l)+")";
		lengthMap.put(k1k2, k1k2dist);
		distMat.put(k1k2, new HashMap<String, Double>());
		for(int i=0;i<klist.length;i++){
			if(!klist[i].equals(key1)&&!klist[i].equals(key2)){

				double newDist=(distMat.get(key1).get(klist[i])+distMat.get(key2).get(klist[i]))/2;

				distMat.get(k1k2).put(klist[i], newDist);
				distMat.get(klist[i]).put(k1k2, newDist);
				
			}
		}
		
		delete(key1);
		delete(key2);
		return new MergeResult(k1k2,k1l,k2l);
	}
	
	public void printDistMat(){
		String[] keys=distMat.keySet().toArray(new String[0]);
		for(int i=0;i<keys.length;i++){
			System.out.print("\t"+keys[i]);
		}
		System.out.println();
		for(int i=0;i<keys.length;i++){
			System.out.print(keys[i]);
			for(int j=0;j<keys.length;j++){
				if(i!=j)System.out.print("\t"+distMat.get(keys[i]).get(keys[j]));
				else System.out.print("\t0");
			}
			System.out.println();
		}
	}
	
	private void calculateDistanceMatrix(){
		ArrayList<String> ids=alg.getIdents();
		
		for(int i=0;i<ids.size();i++){
			distMat.put(ids.get(i).replace(":", "_"),new HashMap<String, Double>());
			for(int j=0;j<ids.size();j++){
				double dist=getHammingDistance(alg.getSequence(i), alg.getSequence(j));
				distMat.get(ids.get(i).replace(":", "_")).put(ids.get(j).replace(":", "_"), dist);
			}
		}
	}
}
