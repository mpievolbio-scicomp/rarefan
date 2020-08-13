package ecoliREP;

import java.io.File;

public class BlastAllExtragenicSpaceFlankingSequences {
	public static void main(String args[]){
		File inputFile=new File(args[0]); //same input as for REP comparison: name\tgenomeFasta\tgenbank\n
		File out=new File(args[1]); //output directory
		int distance=Integer.parseInt(args[2]); //clustering of extragenic spaces...not really applicable
		boolean IsGenome=args[3].equalsIgnoreCase("genome"); //flanking against flanking (false) against genome (true)
		BlastREPs bR=new BlastREPs(out,IsGenome);
		bR.createExtragenicSequenceFiles(inputFile,distance);
		bR.writeFrequency();
		bR.writeFlankingSequences(distance,true,true,0);
		bR.blastFlanking(1e-35); 
	
	}
}
