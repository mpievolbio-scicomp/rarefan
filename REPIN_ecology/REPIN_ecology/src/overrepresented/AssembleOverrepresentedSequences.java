package overrepresented;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import util.DNAmanipulations;
import util.ReadFasta;


public class AssembleOverrepresentedSequences {
	ArrayList<String> sequences=new ArrayList<String>();
	ArrayList<Integer> occurrences=new ArrayList<Integer>();
	int minOverlap=0;
	boolean check=false;
	public static void main(String args[]){
		File over=new File(args[0]);
		File out=new File(args[1]);
		int minOverlap=Integer.parseInt(args[2]);
		AssembleOverrepresentedSequences ao;
		if(args.length>3){
			File genome=new File(args[3]);
			StringBuilder[] Genome=ReadFasta.readFasta(genome).values().toArray(new StringBuilder[0]);
			ao=new AssembleOverrepresentedSequences(over,minOverlap,Genome[0].toString().toUpperCase());
		}else{
			ao=new AssembleOverrepresentedSequences(over,minOverlap);
		}
		ao.write(out);
		
	}
	
	public AssembleOverrepresentedSequences(File in,int MinOverlap){
		minOverlap=MinOverlap;
		read(in);
		assembleAndMerge();
	}
	public AssembleOverrepresentedSequences(File in,int MinOverlap,String Genome){
		minOverlap=MinOverlap;
		read(in);
		assembleAndMerge();
		checkInGenome(Genome);
		check=true;
	}
	
	private void checkInGenome(String Genome){
		for(int i=0;i<sequences.size();i++){
			int j=0;
			int from=0;
			while(-1!=(from=Genome.indexOf(sequences.get(i),from))){
				j++;
				from++;
			}
			from=0;
			while(-1!=(from=Genome.indexOf(DNAmanipulations.reverse(sequences.get(i)),from))){
				j++;
				from++;
			}
			occurrences.add(j);
		}
	}
	
	private void assembleAndMerge(){
		ArrayList<String> assembly=new ArrayList<String>();
		while((assembly=assemble()).size()>0){
			
			sequences.addAll(assembly);
			
			merge();

		}
	
	}
	
	private  ArrayList<String> assemble(){
		ArrayList<String> al=new ArrayList<String>();
		
		for(int i=0;i<sequences.size();i++){
			for(int j=i+1;j<sequences.size();j++){
				String help="";
				if((help=assemble(sequences.get(i),sequences.get(j))).length()>0){
					al.add(help);
				}else{
					if((help=assemble(DNAmanipulations.reverse(sequences.get(i)),sequences.get(j))).length()>0){
						al.add(help);
					}	
				}
			}
		}
		
		return al;
	}
	
	private void merge(){
		ArrayList<String> newAl=new ArrayList<String>();
		while(sequences.size()>0){
			for(int j=1;j<sequences.size();j++){
				if(sequences.get(0).contains(sequences.get(j))){
					sequences.remove(j);
					j=0;
				}else if(sequences.get(j).contains(sequences.get(0))){
					sequences.remove(0);
					j=0;
				}else if(sequences.get(0).contains(DNAmanipulations.reverse(sequences.get(j)))){
					sequences.remove(j);
					j=0;
				}else if(sequences.get(j).contains(DNAmanipulations.reverse(sequences.get(0)))){
					sequences.remove(0);
					j=0;
				}
			}
			newAl.add(new String(sequences.get(0)));
			sequences.remove(0);
			
		}
		
		sequences=newAl;
	}
	
	private  void write(File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<sequences.size();i++){
				
				if(!check){
					bw.write(">seq"+i+"\n");
					bw.write(sequences.get(i)+"\n");
				}
				else bw.write(sequences.get(i)+"\t"+occurrences.get(i)+"\n");
			}
			
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
	private void read(File in){
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				sequences.add(split[0]);
			}
			br.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	

	
	private String assemble(String seq1,String seq2){
				if(seq1.equals(seq2))return seq1;
				if(seq1.length()<seq2.length()){
					String help=new String(seq1);
					seq1=seq2;
					seq2=help;
				}
				if(seq1.contains(seq2))return seq1;
				int score =0;
				String scoreS="";
				for(int l=0;l<seq1.length()+seq2.length();l++){
					if(l<=seq1.length() && l<=seq2.length()){
						
						String c1=seq1.substring(0, l);
						String c2=seq2.substring(seq2.length()-l,seq2.length());
						if(c1.equals(c2) && c1.length()>minOverlap& c1.length()>score){
							scoreS=seq2+seq1.substring(l,seq1.length());
							score=c1.length();
						}
						
					}else if(l>seq2.length() && l<=seq1.length()){
						//this case does not happen cause of the contains condition further up
						
					}else if(l>seq1.length()){
						int start1=l-seq2.length();
						int end1=seq1.length();
						int start2=0;
						int end2=seq1.length()-l+seq2.length();
						String c1=seq1.substring(start1,end1);
						String c2=seq2.substring(start2,end2);
						if(c1.length()!=c2.length()){
							System.err.println(c1+"!="+c2);
							System.exit(-1);
						}
						if(c1.equals(c2)&& c1.length()>minOverlap& c1.length()>score){
							scoreS=seq1+seq2.substring(end2,seq2.length());
							score=c1.length();
						}
					}
				}

		return scoreS;
	}
	
}
