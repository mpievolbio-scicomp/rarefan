package extragenicSequenceSimulation;
//uses an ArrayList of String as input, to create an array list of bitsets as random extragenic regions
//these regions have the same length and GC content as the original regions
//import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import util.EntryExpanded;

public class SimulateExtragenicSequences {
	
	ArrayList<Integer> lengths;
	ArrayList<ArrayList<EntryExpanded<Double,Character>>> ATGCProbs;
	
//	public static void main(String args[]){
//		HashMap<String,StringBuilder> rf=ReadFasta.readFasta(new File(args[0]));
//		String genome=rf.values().toArray(new StringBuilder[0])[0].toString();
//		File artemis=new File(args[1]);
//		GenerateExtragenicSequences ge=new GenerateExtragenicSequences(genome,artemis);
//		System.out.println("Extragenic sequences done.");
//		SimulateExtragenicSequences ses=new SimulateExtragenicSequences(ge.getSequences());
//	}
	
	public SimulateExtragenicSequences(ArrayList<String> original){
		calcParameters(original);
			
	}
	
	private void calcParameters(ArrayList<String> original){
		ATGCProbs=new ArrayList<ArrayList<EntryExpanded<Double,Character>>>();
		lengths=new ArrayList<Integer>();
		for(int i=0;i<original.size();i++){
			String seq=original.get(i);
			lengths.add(seq.length());
			ATGCProbs.add(convertToProbabilitySet(getATGCContent(seq)));
		}
	}
	
	public ArrayList<String> simulate(){
		ArrayList<String> simulation=new ArrayList<String>();
		for(int j=0;j<lengths.size();j++){
			int seqLength=lengths.get(j);
			simulation.add(shuffle(ATGCProbs.get(j),seqLength));
		}
		return simulation;
	}
	
	private HashMap<Character,Double> getATGCContent(String sequence){
		HashMap<Character,Double> ATGC=new HashMap<Character, Double>();
		sequence=sequence.toUpperCase();
		double part=1.0/sequence.length();
		for(int i=0;i<sequence.length();i++){
			Character c=sequence.charAt(i);
			if(!ATGC.containsKey(c)){
				ATGC.put(c, part);
			}else{
				ATGC.put(c, ATGC.get(c)+part);
				
			}
		}
		return ATGC;
		
	}
	
	
	private ArrayList<EntryExpanded<Double,Character>> convertToProbabilitySet(HashMap<Character,Double> ATGC){
		Iterator<Entry<Character,Double>> it=ATGC.entrySet().iterator();
		ArrayList<EntryExpanded<Double, Character>> BaseProbabilities=new ArrayList<EntryExpanded<Double,Character>>();
		double sum=0;
		while(it.hasNext()){
			Entry<Character,Double> e=it.next();
			sum+=e.getValue();
			EntryExpanded<Double, Character> ee=new EntryExpanded<Double, Character>(sum,e.getKey());
			BaseProbabilities.add(ee);
		}
		return BaseProbabilities;
	}
	
	private String shuffle(ArrayList<EntryExpanded<Double,Character>> ATGC,int length){
		StringBuilder sequence=new StringBuilder();
		for(int i=0;i<length;i++){
			double rand=Math.random();
			int j=0;
			while(rand>ATGC.get(j).getKey() && j<ATGC.size()){
				j++;
			}
			sequence.append(ATGC.get(j).getValue());
			
		}
		return sequence.toString();
	}
	

	
	
}
