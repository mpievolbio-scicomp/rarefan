package extragenicSequenceSimulation;

import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

import util.Info;
import util.ReadArtemis;
import util.ReadGenbank;

//Generates extragenic sequences given either an artemis file or genbank file and a genome fasta file (one sequence)
//might want to unify it...too lazy atm

public class GenerateExtragenicSequences {
	ArrayList<String> sequence;
	HashMap<Integer,Integer> map=new HashMap<Integer, Integer>();
//	public static void main(String args[]){
//		HashMap<String,StringBuilder> rf=ReadFasta.readFasta(new File(args[0]));
//		String genome=rf.values().toArray(new StringBuilder[0])[0].toString();
//		File artemis=new File(args[1]);
//		GenerateExtragenicSequences ge=new GenerateExtragenicSequences(genome,artemis);
//	} 
	boolean Gbk;
	ArrayList<Info> extraIntervals=new ArrayList<Info>();
	
	public GenerateExtragenicSequences(String genome){
		sequence=new ArrayList<String>();
		sequence.add(genome);
		map.put(0, genome.length());
		Gbk=false;
		extraIntervals.add(new Info(0,genome.length(),"genome"));
	}
	
	public GenerateExtragenicSequences(String Genome,File Artemis){
		ReadArtemis ra=new ReadArtemis(Artemis);
		BitSet intra=ra.getBoolArray("CDS");
		sequence=new ArrayList<String>();
		map = new HashMap<Integer, Integer>();
		setSequence(intra,Genome);
		Gbk=false;
	}
	public GenerateExtragenicSequences(String Genome,File Genbank,boolean gbk){
		ReadGenbank rb=new ReadGenbank(Genbank);
		ArrayList<Info> intra=rb.getIntervals("CDS");
		sequence=new ArrayList<String>();
		map = new HashMap<Integer, Integer>();
		setSequence(intra,Genome);
		Gbk=true;
	}
	public GenerateExtragenicSequences(String Genome,File Genbank,boolean gbk,int distance){
		ReadGenbank rb=new ReadGenbank(Genbank);
		ArrayList<Info> intra=rb.getIntervals("CDS",distance);
		sequence=new ArrayList<String>();
		map = new HashMap<Integer, Integer>();
		setSequence(intra,Genome);
		Gbk=true;
	}
	public HashMap<Integer,Integer> getMap(){
		return map;
	}
	//not tested!
	private void setSequence(ArrayList<Info> intra,String Genome){
		int lastEnd=0;
		for(int i=0;i<intra.size();i++){
			
			int start=intra.get(i).getStart();
			//search for end!
			while(i+1<intra.size() && intra.get(i).getEnd()>intra.get(i+1).getStart()){
				i++;
			}

			sequence.add(Genome.substring(lastEnd,start));
			extraIntervals.add(new Info(lastEnd+1,start-1,"extragenic"));
			map.put(i, start);
			lastEnd=intra.get(i).getEnd();
		}

	}
	private void setSequence(BitSet intra,String Genome){
		StringBuilder sb=new StringBuilder();
		boolean swit=true;
		int seq=0;
		for(int i=0;i<Genome.length();i++){
			if(intra.size()>i){
				if(intra.get(i)&&swit==false){
					sequence.add(sb.toString());
					sb=new StringBuilder();
					swit=true;
				}else if(false==intra.get(i)){
					if(swit==true){
						map.put(seq,i);
						seq++;
					}
					swit=false;
					sb.append(Genome.charAt(i));
				}
			}else{
				System.err.println("Boolean Array of ReadArtemis not long enough! Wrong genome?");
			}
		}
		sequence.add(sb.toString());

	}	
	
	public ArrayList<Info> getIntervals(){
		if(Gbk){
			return extraIntervals;
		}
		System.err.println("No extragenic intervals set.");
		return null;
	}
	
	public ArrayList<String> getSequences(){
		return sequence;
	}
}
