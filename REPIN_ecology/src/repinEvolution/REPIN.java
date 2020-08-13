package repinEvolution;

import util.*;

public class REPIN {
	int start;
	int end;
	int REPINtype;
	Fasta seq;
	String genomeID;
	int isREPIN;
	int isLargestCluster;
	public REPIN(int start,int end,int REPINtype,Fasta seq,String genomeID,int isREPIN,int isLargestCluster) {
		this.start=start;
		this.end=end;
		this.REPINtype=REPINtype;
		this.seq=seq;
		this.genomeID=genomeID;
		this.isREPIN=isREPIN;
		this.isLargestCluster=isLargestCluster;
	}
	
	public REPIN() {
		this.start=-1;
		this.end=-1;
		this.REPINtype=-1;
		this.seq=null;
		this.genomeID="none";
		this.isREPIN=-1;
		this.isLargestCluster=-1;

	}
	
	public REPIN(String genomeID) {
		this.start=-1;
		this.end=-1;
		this.REPINtype=-1;
		this.seq=null;
		this.genomeID=genomeID;
		this.isREPIN=-1;
		this.isLargestCluster=-1;

	}
	
	public String toString() {
		return genomeID+"\t"+start+"\t"+end+"\t"+REPINtype+"\t"+isREPIN+"\t"+seq.getIdent()+"\t"+seq.getSequence()+"\t"+isLargestCluster+"\n";
	}
	public static String getHeading() {
		return "genomeID\tstart\tend\tREPINtype\tREPIN\tseqIdentifier\tsequence\tisLargestCluster\n";

	}
	
	public Fasta toFasta() {
		String ident=genomeID+","+start+","+end+","+REPINtype+","+isLargestCluster;
		Fasta fas=new Fasta(ident,seq.getSequence());
		return fas;
	}
	
}
