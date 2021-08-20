package repinEvolution;

import java.io.*;
import java.util.*;

import frequencies.*;
import util.Fasta;

//this class is supposed to read and align all fasta (*.fas) files in the inputfolder
//it then identifies all REPINs in each of the genomes defined by REP sequences provided as input?
//it then uses the alignment files to group REPINs from the same position in the genome into clusters
public class MakeREPINClusters {
	
	int numMuts=1;
	double minFrac=0.01;
	ClusterREPINs cr;
	GenomeAlignments ga;
	String reference;
	String[] reps;
	boolean analyseREPIN;
	HashMap<String,ArrayList<Fasta>> repinAlignments;
	public static void main(String[] args) {
		File inFolder=new File(args[0]);
		String[] focalSeeds=Arrays.copyOfRange(args, 1, args.length);
		File out=new File(inFolder+"/repinClusters.txt");
		MakeREPINClusters mrc=new MakeREPINClusters(inFolder, focalSeeds,true);
		mrc.writeClusters(out);
		File fasAlignments=new File(out.getParentFile()+"/fas/");
		fasAlignments.mkdir();
		mrc.writeREPINAlignments(fasAlignments);
	}
	
	public MakeREPINClusters(File inFolder,String[] reps,boolean analyseREPIN) {
		this.analyseREPIN=analyseREPIN;
		File[] genomes=getGenomes(inFolder);
		reference=AlignTwoGenomes.getID(genomes[0]);
		ArrayList<Fasta> refGenome=Fasta.readFasta(genomes[0]);
		ga=alignGenomesToReference(genomes);
		File outFolder=new File(inFolder+"/repinProp/");
		outFolder.mkdir();
		this.reps=reps;
		HashMap<String,REPINpopulations> rp=getREPINProperties(genomes,outFolder);
		clusterREPINsReference(ga,rp,refGenome.get(0).getSequence());
		analyseREPINs(inFolder);
		repinAlignments=cr.getREPINAlignments();
		AnalyseREPINAlignments ara=new AnalyseREPINAlignments(repinAlignments);
		File REPINDivOut=new File(outFolder+"/repinDiversity.txt");
		ara.writePropNotEqual(REPINDivOut);
	}
	
	private void analyseREPINs(File inFolder) {
		AnalyseREPINClusters arc=new AnalyseREPINClusters(cr.rc,reference, getQueries());
		File REPINTransitions=new File(inFolder+"/REPINTransitions.txt");
		arc.writeTransitionStats(REPINTransitions);
	}
	public HashSet<String> getQueries(){
		String[] genomeIDs=ga.getAllKeys();
		HashSet<String> queries=new HashSet<String>();
		for(int i=0;i<genomeIDs.length;i++) {
			if(!genomeIDs[i].equals(reference))queries.add(genomeIDs[i]);
		}
		return queries;
	}
	

	public void writeClusters(File out) {
		cr.writeREPINClusters(out);
	}

	
	public ClusterREPINs getREPINClusters() {
		return cr;
	}
	
	public void writeREPINAlignments(File fasFolder) {
		fasFolder.mkdir();

		String[] names=repinAlignments.keySet().toArray(new String[0]);
		for(int i=0;i<names.length;i++) {
			Fasta.write(repinAlignments.get(names[i]), new File(fasFolder+"/"+names[i]+".fas"));
		}
	}
	
	private void clusterREPINsReference(GenomeAlignments ga,HashMap<String,REPINpopulations> rp,String refGenomeSeq) {
		String[] genomeIDs=ga.getAllKeys();
		cr=new ClusterREPINs(rp.get(reference),reference,refGenomeSeq,reps);
		for(int j=0;j<genomeIDs.length;j++) {
			String g2=genomeIDs[j];
			if(!g2.equals(reference))
				cr.addREPINsToCluster(rp.get(g2), ga.get(reference, g2),g2);
		}
	}

	
	private HashMap<String/*genome ID*/,REPINpopulations> getREPINProperties(File[] genomes,File outFolder){
		HashMap<String,REPINpopulations> rp=new HashMap<String, REPINpopulations>();
		for(int i=0;i<genomes.length;i++) {
			String id=AlignTwoGenomes.getID(genomes[i]);
			rp.put(id, new REPINpopulations());
			for(int j=0;j<reps.length;j++) {
				rp.get(id).addREPINPopulation(reps[j], new REPINProperties(outFolder,id+"_"+j,genomes[i],reps[j].length(),numMuts,minFrac,null,reps[j],false,analyseREPIN, 1));
			}
		}
		
		return rp;
	}
	
	private File[] getGenomes(File inFolder) {
		ArrayList<File> genomes=new ArrayList<File>();
		File[] files=inFolder.listFiles();
		for(int i=0;i<files.length;i++) {
			if(files[i].getName().endsWith(".fas")) {
				genomes.add(files[i]);
			}
		}
		return genomes.toArray(new File[0]);
	}
	
	private GenomeAlignments alignGenomesToReference(File[] genomes) {
		GenomeAlignments alignments=new GenomeAlignments();
		String reference=AlignTwoGenomes.getID(genomes[0]);
		for(int j=1;j<genomes.length;j++) {
			String id2=AlignTwoGenomes.getID(genomes[j]);
			alignments.put(reference,id2,new AlignTwoGenomes(genomes[0],genomes[j]));
		}
		return alignments;
	}
	
}
