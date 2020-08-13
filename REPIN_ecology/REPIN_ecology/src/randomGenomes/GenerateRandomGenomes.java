package randomGenomes;

import java.io.*;
import java.util.*;

import overrepresented.PrintWordFrequencyBothStrands;
import util.*;

public class GenerateRandomGenomes {
	public static void main(String args[]) {
		File inFolder=new File(args[0]);
		File outFolder=new File(args[1]);
		generateRandomSequences(inFolder, outFolder);
	}
	static int wordlength=16;
	public static void generateRandomSequences(File inFolder,File outFolder) {
		File[] list=inFolder.listFiles();
		for(int i=0;i<list.length;i++) {
			File in=list[i];
			if(in.getName().endsWith("fna")) {
				ArrayList<Fasta> fas=Fasta.readFasta(in);
				int seqLength=getSequenceLength(fas);
				double GC=getGCContent(fas);
				String randGenome=DNAmanipulations.generateRandomSequence(seqLength, GC);
				File out=new File(outFolder+"/"+in.getName());
				writeGenome(randGenome,out);
				File wfr=new File(outFolder+"/"+out.getName().split("\\.")[0]+".wfr");
				PrintWordFrequencyBothStrands.writeWords(wordlength,wordlength,out,wfr);
			}
		}
		
		
	}
	
	private static int getSequenceLength(ArrayList<Fasta> fas) {
		int seqLength=0;
		for(int i=0;i<fas.size();i++) {
			seqLength+=fas.get(i).getSequence().length();
			
		}
		return seqLength;
		
	}
	
	private static double getGCContent(ArrayList<Fasta> fas) {
		int gc=0;
		int seqLength=0;
		for(int i=0;i<fas.size();i++) {
			String seq=fas.get(i).getSequence();
			seqLength+=seq.length();
			for(int j=0;j<seq.length();j++) {
				if(seq.charAt(j)=='c'||seq.charAt(j)=='g'||seq.charAt(j)=='C'||seq.charAt(j)=='G') {
					gc+=1;
				}
			}
			fas.get(i).getSequence().length();
			
		}
		return gc/(1.0*seqLength);
		
	}
	
	private static void writeGenome(String seq,File out) {
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
		fas.add(new Fasta("random",seq));
		Fasta.write(fas, out);
	}
	
	
}
