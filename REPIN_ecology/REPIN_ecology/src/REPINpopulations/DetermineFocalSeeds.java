package REPINpopulations;

import java.io.*;
import overrepresented.GroupSeedSequences;
import overrepresented.PrintWordFrequencyBothStrands;
import overrepresented.SelectOverRepresentedWords;

public class DetermineFocalSeeds {
	public static void main(String args[]) {
		File genome=new File(args[0]);
		File outFolder=new  File(args[1]);
		int minFreq=Integer.parseInt(args[2]);
		int wordLength=Integer.parseInt(args[3]);
			DetermineFocalSeeds dfs=new DetermineFocalSeeds(genome,outFolder,minFreq,wordLength,30);
	}
	File genome;
	int wl;
	File genomeFolder;
	String genomeID;
	int minRepFreq;
	int flanking =30;
	String[] seedSequences;
	File outFolder;
	public DetermineFocalSeeds(File genome,File outFolder,int minRepFreq,int wl, int flanking) {
		this.genome=genome;
		this.flanking=flanking;
		this.outFolder=outFolder;
		this.wl=wl;
		genomeFolder=genome.getParentFile();
		genomeID=genome.getName().split("\\.")[0];
		this.minRepFreq=minRepFreq;
		File wordFreq=getWordFrequencies();
		File overRepresentedWords=getOverrepresentedWords(wordFreq);
		seedSequences=getSeedSequences(overRepresentedWords);

	}
	
	public String[] getFocalSeeds() {
		return seedSequences;
	}
	
	private String[] getSeedSequences(File overRepresentedWords) {
		File outFolder=new File(this.outFolder+"/groupSeedSequences/");
		outFolder.mkdir();
		return GroupSeedSequences.groupSeedSequences(genome, overRepresentedWords, outFolder, flanking);
	}
	
	private File getOverrepresentedWords(File wordFreq) {
		File out=new File(outFolder+"/"+genomeID+".overrep");
		SelectOverRepresentedWords.readAndWrite(wordFreq, out, minRepFreq, wl);
		return out;
	}
	
	private File getWordFrequencies() {
		File wordFreqOut=new File(outFolder+"/"+genomeID+".wfr");
		PrintWordFrequencyBothStrands.writeWords(wl, wl, genome, wordFreqOut);
		return wordFreqOut;
	}
	
	
}
