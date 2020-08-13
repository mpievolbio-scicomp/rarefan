package extragenicSequenceSimulation;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map.Entry;

import util.BitSetIndexHash;
import util.DNAmanipulations;
import util.GenerateMutatedSequences;
import util.ReadFasta;
import util.SequencePositions;

public class MutatedSequenceOccurrences {

	ArrayList<ArrayList<Integer>> statistics;
	ArrayList<ArrayList<Integer>> exSpace;
	private int maxMut;
	private int number=0;

	public static void main(String args[]){
		HashMap<String,StringBuilder> rf=ReadFasta.readFasta(new File(args[0]));
		String genome=rf.values().toArray(new StringBuilder[0])[0].toString();
		File artemis=new File(args[1]);
		int maxmutations=Integer.parseInt(args[2]);
		File summary=new File(args[3]);
		int maxsimulations=Integer.parseInt(args[4]);
		File wordFasta=new File(args[5]);
		
		try{
			wordFasta.createNewFile();
			summary.createNewFile();	
		}catch(IOException e){
			System.err.println(e.toString());
			System.exit(1);
		}
		int printMutations=Integer.parseInt(args[6]);
		
		ArrayList<BitSet> words=new ArrayList<BitSet>();
		
		
		boolean full=args[7].equals("full");
		int size;
		size=args[8].length();

		for(int i=8;i<args.length;i++){
			if(args[i].length()!=size){
				System.err.println("Words have to have the same size!");
				System.exit(1);
			}
			words.add(DNAmanipulations.codeDNA(args[i]));
		}

		size=words.get(0).length();
		System.out.println("Start...");
		GenerateExtragenicSequences ge;
		if(artemis.exists()){
			 ge=new GenerateExtragenicSequences(genome,artemis);
		}else{
			ge=new GenerateExtragenicSequences(genome);
		}
		System.out.println("Extragenic sequence generation done.");
		GenerateMutatedSequences gm=new GenerateMutatedSequences(words,maxmutations);
		System.out.println("Original data...");
		ArrayList<String> seqs=ge.getSequences();
		MutatedSequenceOccurrences mso=new MutatedSequenceOccurrences(maxmutations);
		BitSetIndexHash bsih=new BitSetIndexHash(DNAmanipulations.toBitSet(seqs),size,false);
		mso.makeStatistics(bsih,size,gm,full);
		mso.write(summary,false);
		mso.writeWordsWithOccurrence(wordFasta,gm,bsih,printMutations,ge.getMap());
		System.out.println("Random data...");
		SimulateExtragenicSequences ses=new SimulateExtragenicSequences(seqs);
		mso=new MutatedSequenceOccurrences(maxmutations);
		for(int i=0;i<maxsimulations;i++){
			ArrayList<String> randomSeqs=ses.simulate();
			bsih=new BitSetIndexHash(DNAmanipulations.toBitSet(randomSeqs),size,false);
			mso.makeStatistics(bsih,size,gm,full);
		}
		mso.write(summary,true);
	}

	public void writeWordsWithOccurrence(File out,GenerateMutatedSequences gm,BitSetIndexHash bsih,int mut,HashMap<Integer,Integer> map){
		
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			
			ArrayList<BitSet> words=gm.getList(mut);
			
			for(int i=0;i<words.size();i++){
				ArrayList<SequencePositions> sq;
				if((sq=bsih.getPos(words.get(i)))!=null){
					bw.write(">"+DNAmanipulations.decodeDNA(words.get(i))+sq.size()+" Positions:");
					for(int j=0;j<sq.size();j++)
						bw.write(map.get(sq.get(j).sequence)+sq.get(j).position+";");
					bw.write("\n"+DNAmanipulations.decodeDNA(words.get(i))+"\n");
				}
			}
			
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
			
		}
		
	}
	
	

	
	
	/*public void writeSequences(File out){
		
	}*/
	
	public MutatedSequenceOccurrences(int maxmutations){
		maxMut=maxmutations;
		statistics=new ArrayList<ArrayList<Integer>>();
		exSpace=new ArrayList<ArrayList<Integer>>();
		for(int i=0;i<=maxMut;i++){
			statistics.add(new ArrayList<Integer>());
			statistics.get(i).add(0);
			exSpace.add(new ArrayList<Integer>());
			exSpace.get(i).add(0);

		}
	}

	public int getOccurrences(int mut){
		return statistics.get(mut).get(0);
	}
	
	public int getOccupiedExSpaces(int mut){
		return exSpace.get(mut).get(0);
	}
	
	public void makeStatistics(BitSetIndexHash bsih,int size,GenerateMutatedSequences gm,boolean full){
		number++;
		ArrayList<BitSet> newWords=gm.getList(0);
		System.out.println("Trial "+number);
		for(int i=0;i<=maxMut;i++){
			System.out.println("\t"+i+" mutations. Number of sequences: "+newWords.size());
			ArrayList<Integer> stats=checkOverLap(bsih.getPos(newWords),size,full);
			//int occurrences=bsih.getNumber(newWords);
			int occurrences=stats.get(0);
			statistics.get(i).add(occurrences);
			statistics.get(i).set(0,statistics.get(i).get(0)+occurrences);
			int spaces=stats.get(1);
			exSpace.get(i).add(spaces);
			exSpace.get(i).set(0,exSpace.get(i).get(0)+spaces);
			if(i<maxMut){
				newWords=gm.getList(i+1);
			}
		}
	}

	private static ArrayList<Integer> checkOverLap(ArrayList<SequencePositions> pos,int size,boolean full){
		HashMap<Integer,ArrayList<Integer>> hm=new HashMap<Integer,ArrayList<Integer>>();
		ArrayList<Integer> stats=new ArrayList<Integer>();
		HashMap<Integer,Boolean> seqHM=new HashMap<Integer,Boolean>();
		int sum=0;
		for(int i=0;i<pos.size();i++){
			int seq=pos.get(i).sequence;
			seqHM.put(seq,true);
			int posi=pos.get(i).position;
			if(hm.containsKey(seq)){
				hm.get(seq).add( posi);
			}else{
				ArrayList<Integer> al=new ArrayList<Integer>();
				al.add(posi);
				hm.put(seq, al);
			}
		}
		Iterator<Entry<Integer,ArrayList<Integer>>> it=hm.entrySet().iterator();
		while(it.hasNext()){
			Entry<Integer,ArrayList<Integer>> e=it.next();
			TreeMap<Integer,Boolean> tm=new TreeMap<Integer,Boolean>();
			for(int i=0;i<e.getValue().size();i++){
				tm.put(e.getValue().get(i),true);
				
			}
			sum+=checkOverLap(tm,size,full);
		}
		stats.add(sum);
		stats.add(seqHM.size());
		return stats;
		
	}

	private static int checkOverLap(TreeMap<Integer,Boolean> tm,int size,boolean full){
		Integer[] pos=tm.keySet().toArray(new Integer[0]);
		int number=0;
		int fullO=full?1:2;
		for(int i=0;i<pos.length-1;i++){

			if(pos[i]+(size/fullO)<pos[i+1]){
				number++;
			}
		}
		number++;
		return number;
	}
	
	public ArrayList<Integer> getStatistics(int mutation){
		return statistics.get(mutation);
	}

	
	public void write(File Summary,boolean append){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(Summary,append));
			bw.write("\tSummary\t");
			for(int j=1;j<statistics.get(0).size();j++){
				bw.write("Trial "+j+"\t");
			}
			bw.write("\n");
			for(int i=0;i<statistics.size();i++){
				bw.write(i+" mutations:\t");
				bw.write((int)((statistics.get(i).get(0)*1.0)/number)+"\t");
				for(int j=1;j<statistics.get(i).size();j++){
					bw.write(statistics.get(i).get(j)+"\t");
				}
				bw.write("\n");
			}
			bw.close();

		}catch(IOException e){
			System.err.println(e.toString());
			System.exit(1);
		}
	}




}
