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

public class PrintWordFrequency {
	public static void main(String args[]){
		File fasta=new File(args[0]);
		int wl=Integer.parseInt(args[2]);
		int start=Integer.parseInt(args[1]);
		File out=new File(args[3]);
		out.delete();
		HashMap<String,StringBuilder> hm=ReadFasta.readFasta(fasta);
		Iterator<Entry<String, StringBuilder>> it=hm.entrySet().iterator();
		while(it.hasNext()){
			Entry<String,StringBuilder> e=it.next();

				writeWords(start,wl,e.getValue().toString(),out);
			
		}
		
	}

	private static void writeWords(int start,int wl,String genome,File out){

			HashMap<BitSet,Integer> wordsBitSet=new HashMap<BitSet, Integer>();
			HashMap<Long,Integer> wordsLong=new HashMap<Long, Integer>();

			for(int i=start;i<=wl;i++){
				for(int j=0;j<genome.length()-1-i;j++){
					String key=genome.substring(j, j+i);
					BitSet key2BitSet=new BitSet();
					long key2Long=0;
					if(i>31){
						key2BitSet=DNAmanipulations.codeDNA(key.toUpperCase());
						if(key2BitSet==null){
							continue;
						}
						if(wordsBitSet.containsKey(key2BitSet)){
							wordsBitSet.put(key2BitSet,wordsBitSet.get(key2BitSet)+1);
						}else{
							wordsBitSet.put(key2BitSet,1);
						}
					}else{
						key2Long=DNAmanipulations.codeDNALong(key.toUpperCase());
						if(key2Long==-1){
							continue;
						}
						if(wordsLong.containsKey(key2Long)){
							wordsLong.put(key2Long,wordsLong.get(key2Long)+1);
						}else{
							wordsLong.put(key2Long,1);
						}
					}
					

				}
				if(i>31)write(wordsBitSet,out,i);
				else writeQuick(wordsLong,out,i);
				wordsBitSet.clear();
				wordsLong.clear();
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
	private static void writeQuick(HashMap<Long,Integer> words,File out,int wl){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out,true));
			Iterator<Entry<Long,Integer>> it=words.entrySet().iterator();
			while(it.hasNext()){
				Entry<Long,Integer> e=it.next();
				bw.write(DNAmanipulations.decodeDNA(e.getKey(),wl)+" "+e.getValue()+"\n");
			}
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}

	
}
