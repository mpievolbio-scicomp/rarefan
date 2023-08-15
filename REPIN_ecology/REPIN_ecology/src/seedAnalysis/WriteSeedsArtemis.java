package seedAnalysis;
import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;

import util.BitSetIndexHash;
import util.DNAmanipulations;
import util.GenerateMutatedSequences;
import util.ReadBitSetSearch;
import util.ReadFasta;
import util.SortPositions;
import util.WriteGenomeAnnotation;
import extragenicSequenceSimulation.GenerateExtragenicSequences;


public class WriteSeedsArtemis {
	public static void main(String args[]){
		File genome=new File(args[0]);
		String fasta=ReadFasta.readFasta(genome).values().toArray(new StringBuilder[0])[0].toString();
		File artemis=new File(args[1]);
		ArrayList<BitSet> words=new ArrayList<BitSet>();
		int mutations=Integer.parseInt(args[2]);
		File outArt=new File (args[3]);
		for(int i=4;i<args.length;i++){
			words.add(DNAmanipulations.codeDNA(args[i]));
		}
		int size=words.get(0).length();
		GenerateExtragenicSequences ge=new GenerateExtragenicSequences(fasta,artemis);
		System.out.println("Generated extragenic sequences...");
		BitSetIndexHash bsih=new BitSetIndexHash(DNAmanipulations.toBitSet(ge.getSequences()),size,ge.getMap(),false);
		System.out.println("Initialised index...");
		GenerateMutatedSequences gm=new GenerateMutatedSequences(words,mutations);
		System.out.println("Mutated sequences...");
		ArrayList<Integer> seqPositions=subtractOverLaps(bsih.getPosMap(gm.getList(mutations)),size/2);
		ReadBitSetSearch rbs=new ReadBitSetSearch(seqPositions,size/2);
		WriteGenomeAnnotation.writeTab(seqPositions, rbs.getEnd(), rbs.getQuery(), outArt);
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
	

}
