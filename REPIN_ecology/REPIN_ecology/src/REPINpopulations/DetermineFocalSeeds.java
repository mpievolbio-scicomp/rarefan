package REPINpopulations;

import java.io.*;
import java.util.*;

import identifyRAYTs.BlastRAYTs;
import overrepresented.GroupSeedSequences;
import overrepresented.PrintWordFrequencyBothStrands;
import overrepresented.SelectOverRepresentedWords;
import util.*;

public class DetermineFocalSeeds {
	public static void main(String args[]) {
		File genome=new File(args[0]);
		int minFreq=Integer.parseInt(args[1]);
		int wordLength=Integer.parseInt(args[2]);
		File queryRAYT=new File(args[3]);
		if(args.length>4) {
			String legacyBlastPerlLocation=args[4];
			DetermineFocalSeeds dfs=new DetermineFocalSeeds(genome,minFreq,wordLength);
		}else {
			DetermineFocalSeeds dfs=new DetermineFocalSeeds(genome,minFreq,wordLength);

		}
	}
	File genome;
	int wl;
	File genomeFolder;
	String genomeID;
	int minRepFreq;
	int flanking =30;
	String[] seedSequences;
	public DetermineFocalSeeds(File genome,int minRepFreq,int wl) {
		this.genome=genome;
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
		File outFolder=new File(genomeFolder+"/groupSeedSequences/");
		outFolder.mkdir();
		return GroupSeedSequences.groupSeedSequences(genome, overRepresentedWords, outFolder, flanking);
	}
	
	private File getOverrepresentedWords(File wordFreq) {
		File out=new File(genomeFolder+"/"+genomeID+".overrep");
		SelectOverRepresentedWords.readAndWrite(wordFreq, out, minRepFreq, wl);
		return out;
	}
	
	private File getWordFrequencies() {
		File wordFreqOut=new File(genomeFolder+"/"+genomeID+".wfr");
		PrintWordFrequencyBothStrands.writeWords(wl, wl, genome, wordFreqOut);
		return wordFreqOut;
	}
	
	
}
