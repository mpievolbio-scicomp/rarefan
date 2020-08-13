package overrepresented;

import java.io.*;
import java.util.*;

import util.Info;
import util.ReadFasta;
import yafMSignificance.GetFlankingGenes;

//determines the number of occurrences of a certain group of sequences (output of GroupSeedSequences.java) by counting the number of clusters that are formed
public class GroupOccurrences {

	public static void main(String args[]){
		File group=new File(args[0]);
		String genome=ReadFasta.readFasta(new File(args[1])).values().toArray(new StringBuilder[0])[0].toString();
		ArrayList<String> words=readGroup(group);
		GroupOccurrences go=new GroupOccurrences(genome);
		for(int i=0;i<words.size();i++){
			ArrayList<Info> intervals=GetFlankingGenes.find(genome, words.get(i));	
			System.out.print(words.get(i)+"\t");
			go.setOccupied(intervals);
			System.out.println(go.countOccurrences());
		}
	}
	String genome;
	BitSet occupied;
	public  GroupOccurrences(String Genome){
		genome=Genome;
		occupied=new BitSet(Genome.length());
	}
	
	public int countOccurrences(){
		int count=0;
		boolean prev=false;
		for(int i=0;i<occupied.length();i++){
			if(occupied.get(i)!=prev&& occupied.get(i)){
				count++;
				prev=occupied.get(i);
			//System.out.println(i);
			}else if(!occupied.get(i)){
				prev=occupied.get(i);
			}
		}
		return count;
	}
	
	public void setOccupied(ArrayList<Info> intervals){
		for(int i=0;i<intervals.size();i++){
			setOccupied(intervals.get(i));
		}
	}
	public void setOccupied(Info interval){
			occupied.set(interval.getStart(), interval.getEnd());
	}	
	public static ArrayList<String> readGroup(File group){
		ArrayList<String> al=new ArrayList<String>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(group));
			String line="";
			while((line=br.readLine())!=null){
				al.add(line.split("\\s+")[0]);
			}
				
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return al;
	}
	
}
