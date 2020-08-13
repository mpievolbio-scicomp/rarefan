package overrepresented;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

import util.DNAmanipulations;

//this program computes the p-Value and word occurrence of all sequences in a fasta file that are of a specified length (atm 16bp)
//it uses the program GetFrequencyBelowPvalue.java for calculating the pValue

//as input are required: folder for word frequency files
//a fasta file where the first word in the description has to be the same as the filename before the extension in the wordFrequency folder

//words that are longer than 20bp are chopped into 16bp long words and are subsequently analysed

//input: fastaFile wordFrequencyFileFolder OutputFilename
//output: FastaDescription\tword\toccurrences\tp-Value
public class PValueWords {
	public static void main(String args[]){
		File fasta=new File(args[0]);
		File wordFreqs=new File(args[1]);
		int wordlength=Integer.parseInt(args[2]);
		File out=new File(args[3]);
		//HashMap<BitSet,Integer> Occurrences=null;
		ArrayList<String> ids=new ArrayList<String>();
		ArrayList<String> words=new ArrayList<String>();
		readFasta(fasta,ids,words,wordlength);
		GetFrequencyBelowPvalue gp = null;
		File frequencyE=new File("");
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<ids.size();i++){
				File frequency=new File(wordFreqs+"/"+ids.get(i).split("\\s+")[0]+".out");
				
				if(gp==null || !frequency.equals(frequencyE)){
					//Occurrences=null;
					gp=new GetFrequencyBelowPvalue(frequency,wordlength);
					//Occurrences=readMap(frequency,wordlength);
				}
				int occ=gp.getFreq(words.get(i));
				double pValue=gp.getpValue(occ, words.get(i).length());
				bw.write(ids.get(i)+"\t"+words.get(i)+"\t"+occ+"\t"+pValue+"\n");
				System.out.println(ids.get(i)+"\t"+words.get(i)+"\t"+occ+"\t"+pValue);
				frequencyE=frequency;
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		
	}
	
	public static void readFasta(File fas,ArrayList<String> ids,ArrayList<String> words,int wl){
		try{
			BufferedReader br = new BufferedReader(new FileReader(fas));
			String line="";
			String id="";
			BufferedWriter bw=new BufferedWriter(new FileWriter(new File("/home/frederic/auckland/fastaSeqs/palindromesMorethan20bp.fas")));
			while((line=br.readLine())!=null){
				if(line.startsWith(">")){
					id=line.substring(1);
				}else {
					line=line.replace("\n","");
					line=line.replace("\r","");
					if(line.length()==wl){
						ids.add(id);
						words.add(line);
					}else{
						bw.write(">"+id+"\n"+line+"\n");
						chopWord(line,id,ids,words,wl);
					}
				}
			}
			br.close();
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		
	}
	
	public static void chopWord(String word,String id,ArrayList<String> ids,ArrayList<String> words,int wl){
		for(int i=1;i<=word.length()-wl-1;i++){
			ids.add(id+"."+i);
			words.add(word.substring(i-1,i+wl-1));
			//System.out.println(word.substring(i-1,i+15)+" "+word);
		}
	}
	
	public static int getOccurrences(HashMap<BitSet,Integer> Occ,String word){
		BitSet code=DNAmanipulations.codeDNA(word);
		BitSet rev=DNAmanipulations.reverse(code);
		if(Occ.containsKey(code)){
			return Occ.get(code);
		}else if(Occ.containsKey(rev)){
			return Occ.get(rev);
		}else return 0;
	}
	public static HashMap<BitSet,Integer> readMap(File in,int wl){
		HashMap<BitSet,Integer> hm=new HashMap<BitSet, Integer>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				if(split[0].length()!=wl)continue;
				BitSet coded=DNAmanipulations.codeDNA(split[0]);
				hm.put(coded, Integer.parseInt(split[1]));
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return hm;
	}
}
