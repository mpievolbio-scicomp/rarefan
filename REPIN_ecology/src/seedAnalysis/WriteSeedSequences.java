package seedAnalysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;


import util.BitSetIndexHash;
import util.DNAmanipulations;
import util.Fasta;
import util.GenerateMutatedSequences;
import util.ReadBitSetSearch;
import util.ReadFasta;
import util.SortPositions;
import extragenicSequenceSimulation.GenerateExtragenicSequences;

public class WriteSeedSequences {
	public static void main(String[] args){
		File genome=new File(args[0]);
		int mutations=Integer.parseInt(args[2]);
		File out=new File(args[3]);
		File artemis=new File(args[1]);
		//int minOcc=Integer.parseInt(args[4]);
		ArrayList<String> words=new ArrayList<String>();
		for(int i=4;i<args.length;i++){
			words.add(args[i]);
		}
		writeSeedSequencesOld(words,genome,artemis,mutations,out);
		
	}
	
	private static ArrayList<BitSet> transformWords(ArrayList<String> words){
		ArrayList<BitSet> bs=new ArrayList<BitSet>();
		for(int i=0;i<words.size();i++){
			bs.add(DNAmanipulations.codeDNA(words.get(i)));
		}
		return bs;
	}
	
	public static void writeSeedSequences(ArrayList<String> wordsRaw,File genome,int mutations,File out){
		ArrayList<BitSet> words=transformWords(wordsRaw);
		int size=words.get(0).length();
		ArrayList<String> fasta=toString(Fasta.readFasta(genome));
		BitSetIndexHash bsih=new BitSetIndexHash(DNAmanipulations.toBitSet(fasta),size,true);
		System.out.println("Initialised index...");
		GenerateMutatedSequences gm=new GenerateMutatedSequences(words,mutations);
		System.out.println("Mutated sequences...");
		writeOccurrences(bsih,gm.getList(mutations),out);
		System.out.println("Done.");
	}
	
	public static void writeSeedSequencesConnected(ArrayList<String> wordsRaw,File genome,int mutations,File out){
		ArrayList<BitSet> words=transformWords(wordsRaw);
		int size=words.get(0).length();
		ArrayList<String> fasta=toString(Fasta.readFasta(genome));
		System.out.println("Initialize index..."+genome );
		BitSetIndexHash bsih=new BitSetIndexHash(DNAmanipulations.toBitSet(fasta),size,true);
		System.out.println("Get all connected sequences...");
		ArrayList<BitSet> connected=getAllConnectedSequences(bsih,words,mutations);
		connected=minimizeDistance(connected,words.get(0));
		System.out.println("Write occurrences...");
		writeOccurrences(bsih,connected,out);
		System.out.println("Done.");
	}
	
	private static ArrayList<BitSet> minimizeDistance(ArrayList<BitSet> all,BitSet ref){
		HashSet<BitSet> newSet=new HashSet<BitSet>();
		for(int i=0;i<all.size();i++) {
			BitSet forward=all.get(i);
			BitSet reverse=DNAmanipulations.reverse(forward);
			int distF=getDistance(ref,forward);
			int distR=getDistance(ref,reverse);
			BitSet add=distF<distR?forward:reverse;
			if(!newSet.contains(add)) {
				newSet.add(add);
			}
		}
		return new ArrayList<BitSet>(newSet);
	}
	
	private  static int getDistance(BitSet one,BitSet two) {
		String oneS=DNAmanipulations.decodeDNA(one);
		String twoS=DNAmanipulations.decodeDNA(two);
		return getDifferences(oneS,twoS);
	}
	
	private static int getDifferences(String one,String two) {
		int diff=0;
		for(int i=0;i<one.length();i++) {
			if(one.charAt(i)!=two.charAt(i)) {
				diff++;
			}
		}
		return diff;
	}
	
	private static ArrayList<BitSet> getAllConnectedSequences(BitSetIndexHash bsih,ArrayList<BitSet> words,int mutations){
		HashSet<BitSet> all=new HashSet<BitSet>();
		for(int i=0;i<words.size();i++){
			getConnectedMutations(words.get(i),bsih,all,mutations);
		}
		
		return toArrayList(all);
	}
	
	private static ArrayList<BitSet> toArrayList(HashSet<BitSet> all){
		ArrayList<BitSet> list=new ArrayList<BitSet>();
		BitSet[] array=all.toArray(new BitSet[0]);
		for(int i=0;i<array.length;i++){
			list.add(array[i]);
		}
		
		return list;
		
	}
	
	private static void getConnectedMutations(BitSet word,BitSetIndexHash bsih,HashSet<BitSet> all,int mutations){
		ArrayList<BitSet> words=new ArrayList<BitSet>();
		words.add(word);
		GenerateMutatedSequences gm=new GenerateMutatedSequences(words,mutations);
		ArrayList<BitSet> results=gm.getList(mutations);
		for(int i=0;i<results.size();i++){
			if(bsih.getNumber(results.get(i))>0){
				/*System.out.println(DNAmanipulations.decodeDNA(word)+" "+DNAmanipulations.decodeDNA(results.get(i))+" "+bsih.getNumber(results.get(i))+" "+mutations+" "+all.size());
				if(all.size()>10000) {
					System.exit(0);
				}*/
				if(!all.contains(results.get(i))&&!all.contains(DNAmanipulations.reverse(results.get(i)))){
					all.add(results.get(i));
					getConnectedMutations(results.get(i),bsih,all,mutations);
				}
			}
		}
	}
	
	public static void writeOccurrences(BitSetIndexHash bsih,ArrayList<BitSet> words,File out){
		words=GenerateMutatedSequences.eliminateDoubles(words);
		ArrayList<Fasta> wordFasta=new ArrayList<Fasta>();
		for(int i=0;i<words.size();i++){
			int count=bsih.getNumber(words.get(i));
			if(count>0){
				wordFasta.add(new Fasta("wordOccurs"+count,DNAmanipulations.decodeDNA(words.get(i))));
			}
		}
		Fasta.write(wordFasta, out);
	}
	
	private static ArrayList<String> toString(ArrayList<Fasta> fas){
		ArrayList<String> list=new ArrayList<String>();
		for(int i=0;i<fas.size();i++){
			String seq=fas.get(i).getSequence();
			String newseq=seq.replaceAll("[^ATGC]", "");
			newseq=newseq.replace("N", "");
			list.add(newseq);
			
		}
		return list;
	}
	
	
	public static void writeSeedSequencesOld(ArrayList<String> wordsRaw,File genome,File artemis,int mutations,File out){
		ArrayList<BitSet> words=transformWords(wordsRaw);
		int minOcc=1;
		int size=words.get(0).length();
		String fasta=ReadFasta.readFasta(genome).values().toArray(new StringBuilder[0])[0].toString();
		GenerateExtragenicSequences ge=new GenerateExtragenicSequences(fasta,artemis);
		System.out.println("Generated extragenic sequences...");
		BitSetIndexHash bsih=new BitSetIndexHash(DNAmanipulations.toBitSet(ge.getSequences()),size,ge.getMap(),true);
		System.out.println("Initialised index...");
		GenerateMutatedSequences gm=new GenerateMutatedSequences(words,mutations);
		System.out.println("Mutated sequences...");
		ArrayList<Integer> seqPositions=subtractOverLaps(bsih.getPosMap(gm.getList(mutations)),size/2);
		ReadBitSetSearch rbs=new ReadBitSetSearch(seqPositions,size/2);
		writeOld(rbs.getStart(), rbs.getEnd(), rbs.getQuery(),fasta, out,minOcc,bsih);
		System.out.println("Done.");
	}
	
	public static ArrayList<Integer> subtractOverLaps(ArrayList<Integer> pos,int size){
		SortPositions sa=new SortPositions(pos);
		ArrayList<Integer> temp=sa.getList();
		ArrayList<Integer> positions=new ArrayList<Integer>();
		for(int i=0;i<temp.size()-1;i++){
		if(temp.get(i)+size<temp.get(i+1)){
		positions.add(temp.get(i));
		}
		}
		if(temp.size()>0)positions.add(temp.get(temp.size()-1));
		return positions;
	}
	

	
	public static void writeOld(ArrayList<Integer> start,ArrayList<Integer> end,ArrayList<String> id,String genome,File out,int minOcc,BitSetIndexHash bsih){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			HashMap<String,Integer> exists=new HashMap<String, Integer>();
			int sum=0;
			for(int i=0;i<start.size();i++){
				String word="";
				if(start.get(i)<end.get(i))
					word=genome.substring(start.get(i),end.get(i)).toUpperCase();
				else
					word=DNAmanipulations.reverse(genome.substring(end.get(i),start.get(i)).toUpperCase());
				if(!exists.containsKey(word)){
					int occ=bsih.getNumber(DNAmanipulations.codeDNA(word));
					if(minOcc<=occ){
						sum+=occ;
						exists.put(word,occ);
						bw.write(">"+word+occ+"\n"+word+"\n");
					}
				}
			}
			bw.close();
			
			System.out.println(sum + "total occurrences.");
			System.out.println("Bp frequencies:");
			System.out.println(writeBPFreq(exists));
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
	public static String writeBPFreq(HashMap<String,Integer> wordfreq){
		String[] words=wordfreq.keySet().toArray(new String[0]);
		ArrayList<HashMap<Character,Integer>> position=new ArrayList<HashMap<Character,Integer>>();
		for(int i=0;i<words[0].length();i++){
			HashMap<Character,Integer> temp=new HashMap<Character, Integer>();
			for(int j=0;j<words.length;j++){
				char x=words[j].toUpperCase().charAt(i);
				if(temp.containsKey(x)){
					temp.put(x,temp.get(x)+wordfreq.get(words[j]));
				}else{
					temp.put(x, wordfreq.get(words[j]));

				}
				
			}
			position.add(temp);
		}
		return print(position);
	}
	
	public static String print(ArrayList<HashMap<Character,Integer>> pos){
		StringBuilder sb=new StringBuilder();
		for(int i=0;i<pos.size();i++){
			Iterator<Entry<Character,Integer>> it=pos.get(i).entrySet().iterator();
			sb.append("Position "+i+":\n");
			while(it.hasNext()){
				Entry<Character,Integer> e=it.next();
				sb.append(e.getKey()+" "+e.getValue()+"\n");
			}
			
		}
		return sb.toString();
	}

}
