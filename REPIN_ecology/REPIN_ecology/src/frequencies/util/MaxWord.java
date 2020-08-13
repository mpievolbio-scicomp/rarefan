package frequencies.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import overrepresented.PrintWordFrequencyBothStrands;
import util.Fasta;
import util.Find;

public class MaxWord {
	File wfr;
	File maxWordFile;
	String genomeID;
	File outFolder;
	String wfrExt=".wfr";
	String maxWordExt=".mw";

	String maxWord;
	int frequency;

	public MaxWord(File fas,int wordlength,File outFolder,String genomeID){
		this.genomeID=genomeID;
		this.outFolder=outFolder;
		wfr=new File(outFolder+"/"+genomeID+wfrExt);
		System.out.println(wfr);
		if(!wfr.exists())PrintWordFrequencyBothStrands.writeWords(wordlength, wordlength, fas, wfr);
		maxWordFile=new File(outFolder+"/"+genomeID+maxWordExt);
		calculateMaxWords(maxWordFile);
	}
	public MaxWord(File fas,String word,int wordlength,File outFolder,String genomeID){
		this.genomeID=genomeID;
		this.outFolder=outFolder;
		maxWordFile=new File(outFolder+"/"+genomeID+maxWordExt);
		maxWord=word;
		frequency = calculateFrequency(fas,word);
		
		calculateMaxWords(maxWordFile);
	}
	
	private int calculateFrequency(File fasFile,String word){
		int f=0;
		ArrayList<Fasta> fas=Fasta.readFasta(fasFile);
		for(int i=0;i<fas.size();i++){
			
			f+=Find.getOccurrences(word,fas.get(i).getSequence());
		}
		return f;
	}
	
	public int getFrequency(){
		return frequency;
	}
	
	
	public String getMaxWord(){
		return maxWord;
	}
	private void calculateMaxWords(File out){
		if(!out.exists()){
			writeMaxWord(out);
		}else{
			readMaxWord(out);
		}
	}

	private void readMaxWord(File in){
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line;
			while((line=br.readLine())!=null){
				String split[]=line.split("\t");
				frequency=Integer.parseInt(split[1]);
				maxWord=split[0];
				
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private void writeMaxWord(File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			if(wfr!=null&&wfr.exists())calculateMaxWord(wfr);
			bw.write(maxWord+"\t"+frequency+"\n");
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}

	}
	private void calculateMaxWord(File in){
		int max=0;
		String maxword="";
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String split[]=line.split("\\s+");
				int wc=Integer.parseInt(split[1]);
				String word=split[0];
				if(wc>max){
					max=wc;
					maxword=word;
				}
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		maxWord= maxword;
		frequency=max;
	}
}
