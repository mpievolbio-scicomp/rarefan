package overrepresented;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import util.DNAmanipulations;
import util.ReadFasta;

public class PrintWordFrequencyBothStrands {
	public static void main(String args[]){
		File fasta=new File(args[0]);
		int wl=Integer.parseInt(args[2]);
		int start=Integer.parseInt(args[1]);
		File out=new File(args[3]);
		
		writeWords(start,wl,fasta,out);
			
		
	}
	


	public static void writeWords(int start,int wl,File fasta,File out){
		out.delete();
		HashMap<String,StringBuilder> hm=ReadFasta.readFasta(fasta);
		Iterator<Entry<String, StringBuilder>> it=hm.entrySet().iterator();
		StringBuilder genome=new StringBuilder();
		while(it.hasNext()){
			Entry<String,StringBuilder> e=it.next();
			genome.append("|"+e.getValue());
		}
			HashMap<BitSet,Integer> wordsBitSet=new HashMap<BitSet, Integer>();

			for(int i=start;i<=wl;i++){
				for(int j=0;j<genome.length()-1-i;j++){
					String key=genome.substring(j, j+i);
					BitSet key2BitSet=new BitSet();
					key2BitSet=DNAmanipulations.codeDNA(key.toUpperCase());
					if(key2BitSet==null){
						continue;
					}
					BitSet complement=DNAmanipulations.reverse(key2BitSet);
					if(wordsBitSet.containsKey(key2BitSet) ){
						wordsBitSet.put(key2BitSet,wordsBitSet.get(key2BitSet)+1);
					}else if( wordsBitSet.containsKey(complement)){
						wordsBitSet.put(complement,wordsBitSet.get(complement)+1);

					}else{
						wordsBitSet.put(key2BitSet,1);
					}
				}
				write(wordsBitSet,out,i);
				wordsBitSet.clear();
			}
	}


	private static void write(HashMap<BitSet,Integer> words,File out,int wl){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out,true));
			Iterator<Entry<BitSet,Integer>> it=words.entrySet().iterator();
			while(it.hasNext()){
				Entry<BitSet,Integer> e=it.next();
				bw.write(DNAmanipulations.decodeDNA(e.getKey())+" "+e.getValue()+"\n");
			}
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}


	
}
