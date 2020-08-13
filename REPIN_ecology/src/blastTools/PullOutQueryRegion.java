package blastTools;

import java.io.*;
import java.util.*;
import util.*;
//pulls out genes of the original query length usually from a big genome database given a blast input file (m8) and a database fasta file
public class PullOutQueryRegion {
	
	public static void main(String[] args){
		File blast=new File(args[0]);
		File db=new File(args[1]);
		int queryLength=Integer.parseInt(args[2]);
		File out=new File (args[3]);
		ReadBlast rb=new ReadBlast(blast);
		Fasta.write(getFastaSequences(db, getGeneIntervals(rb, queryLength, true),true),out);
	}
	

	
	public  static ArrayList<Fasta> getFastaSequences(File database,HashMap<String,ArrayList<Info>> sequenceIntervals,boolean translate){
		ArrayList<Fasta> sequences=new ArrayList<Fasta>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(database));
			String line="";
			boolean found=false;
			StringBuilder sequence=new StringBuilder();
			ArrayList<Info> currentIntervals=null;
			String ident=null;
			while((line=br.readLine())!=null){
				if(line.startsWith(">")){
					if(sequence.length()>0 && currentIntervals!=null){
						String seq=sequence.toString();
						for(int i=0;i<currentIntervals.size();i++){
							boolean inverted=currentIntervals.get(i).info.contains("inverted");
							int start=currentIntervals.get(i).getStart();
							int end=currentIntervals.get(i).getEnd();
							String subseq=seq.substring(start-1,end );

							subseq=!inverted?subseq:DNAmanipulations.reverse(subseq);
							subseq=translate?DNAmanipulations.translate(subseq, DNAmanipulations.code()):subseq;
							String newIdent=inverted?ident+"_"+start+"_"+end+"_complement":ident+"_"+start+"_"+end;
							sequences.add(new Fasta(newIdent,subseq));
						}
						
					}
					sequence=new StringBuilder();
					String split[]=line.split("\\s+|>");
					
					if(sequenceIntervals.containsKey(split[1])){
						found=true;
						ident=split[1];
						currentIntervals=getUniqueIntervals(sequenceIntervals.get(split[1]));
					}else{
						currentIntervals=null;
						found = false;
						ident=null;
					}
				}else if(found){
					sequence.append(line.replaceAll("\\s+", ""));
				}
			}
			if(sequence.length()>0 && currentIntervals!=null){
				String seq=sequence.toString();
				for(int i=0;i<currentIntervals.size();i++){
					boolean inverted=currentIntervals.get(i).info.contains("inverted");
					int start=currentIntervals.get(i).getStart();
					int end=currentIntervals.get(i).getEnd();
					String subseq=seq.substring(start-1,end-1 );
					subseq=!inverted?subseq:DNAmanipulations.reverse(subseq);
					subseq=translate?DNAmanipulations.translate(subseq, DNAmanipulations.code()):subseq;
					String newIdent=inverted?ident+"_"+start+"_"+end+"_complement":ident+"_"+start+"_"+end;
					sequences.add(new Fasta(newIdent,subseq));
				}
				
			}
		}catch(IOException e){
			e.printStackTrace();
		}
		
		
		return sequences;
	}
	
	//only for proteinQuery and nucleotide DB or nuc nuc or prot prot
	public  static HashMap<String,ArrayList<Info>> getGeneIntervals(ReadBlast rb,int queryLength,boolean protQueryNucDB){
		HashMap<String,ArrayList<Info>> entries=new HashMap<String, ArrayList<Info>>();
		for (int i=0;i<rb.getStartDB().size();i++){
			int startQ=rb.getStartQuery().get(i);
			int endQ=rb.getEndQuery().get(i);
			int startDB=rb.getStartDB().get(i);
			int endDB=rb.getEndDB().get(i);
			boolean inverted=startDB>endDB;
			int startDist=inverted?queryLength-endQ:startQ-1;
			int endDist=!inverted?queryLength-endQ:startQ-1;
			if(protQueryNucDB){
				startDist=3*startDist;
				endDist=3*endDist;
			}
			int newStart=inverted?endDB-startDist:startDB-startDist;
			int newEnd=inverted?startDB+endDist:endDB+endDist;
			Info entry=new Info(newStart,newEnd,inverted?"inverted":"");
			if(entries.containsKey(rb.getDatabase().get(i))){
				entries.get(rb.getDatabase().get(i)).add(entry);
			}else{
				ArrayList<Info> genes=new ArrayList<Info>();
				genes.add(entry);
				entries.put(rb.getDatabase().get(i), genes);
			}
//			int expectedLength=queryLength*3;
//			int length=newEnd-newStart;
//			if(length!=expectedLength){
//				System.err.println(rb.get(i));
//				System.err.println("Error, length incorrect: "+length+" as opposed to the correct: "+expectedLength);
//			}
			
		}
		
		
		return entries;
	}
	
	public static ArrayList<Info> getUniqueIntervals(ArrayList<Info> intervals){
		ArrayList<Info> unique=new ArrayList<Info>();
		InfoTree itree=new InfoTree();
		for(int i=0;i<intervals.size();i++){
			ArrayList<Info> al=new ArrayList<Info>();
			itree.search(intervals.get(i), al);
			if(al.size()==0){
				itree.insert(intervals.get(i));
				unique.add(intervals.get(i));
			}
		}
		return unique;
	}
	
	
}
