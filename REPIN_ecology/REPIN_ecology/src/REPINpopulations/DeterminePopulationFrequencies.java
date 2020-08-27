package REPINpopulations;
import java.io.*;
import java.util.*;

import frequencies.REPINProperties;
import identifyRAYTs.BlastRAYTs;
import util.*;
//the idea is to determine the REPIN populations that are present in a number of focal strains
//(runnning word freqs, identify groups etc) at first I will supply the sequence seeds (focal sequences)
//for each sequence group we determine the frequency in each of the given strains
//this data will then be displayed on a tree using R
import util.phylogenetics.RunTreePrograms;

public class DeterminePopulationFrequencies {
	//requires mcl, andi, clustDist and BLAST+
	String focalSeeds[];
	ArrayList<File> genomes;
	File inFolder;
	int numMuts=1;
	double minFrac=0.01;
	//distance from repin to rayt, if within vicinity then repin cluster is associated with that rayt
	String legacyBlastPerlLocation;
	File queryRAYT;
	File genomeFolder;
	String e;
	boolean analyseREPIN;
	HashMap<String/*genomes*/,HashMap<String/*focal seed*/,Integer/*pop size*/>> results=new HashMap<String,HashMap<String,Integer>>();
	public static void main(String args[]) {
		File inFolder=new File(args[0]);
		String focalSeedGenome=args[1];
		int minRepFreq=Integer.parseInt(args[2]);
		int wordlength=Integer.parseInt(args[3]);
		File queryRAYT=new File(args[4]);
		File treeFile=new File(args[5]);
		String evalue=args[6];
		boolean analyseREPIN=args[7].equalsIgnoreCase("true");
		File out=new File(inFolder+"/results.txt");
		DeterminePopulationFrequencies dpf;
		String program="tblastn";
		if(args.length>8) {
			String legacyBlastPerlLocation=args[8];
			dpf=new DeterminePopulationFrequencies(inFolder, focalSeedGenome,minRepFreq,wordlength,queryRAYT,program,treeFile,legacyBlastPerlLocation,evalue,analyseREPIN);

		}else {
			dpf=new DeterminePopulationFrequencies(inFolder, focalSeedGenome,minRepFreq,wordlength,queryRAYT,program,treeFile,"",evalue,analyseREPIN);
		}


		dpf.print(out);
	}
	
	
	
	public DeterminePopulationFrequencies(File inFolder,String focalSeedGenome,int minRepFreq,int wordlength,File queryRAYT,String program,File treeFile,String legacyBlastPerlLocation,String evalue,boolean analyseREPIN){
		this.inFolder=inFolder;
		genomes=getFiles(inFolder);
		this.legacyBlastPerlLocation=legacyBlastPerlLocation;
		this.queryRAYT=queryRAYT;
		this.focalSeeds=getFocalSeeds(focalSeedGenome,minRepFreq,wordlength);
		this.genomeFolder=inFolder;
		this.analyseREPIN=analyseREPIN;
		e=evalue;

		calculateResults();

		BlastRAYTs.runProgram(inFolder, queryRAYT, inFolder, e, program, getREPtype(), "yafM_relatives.fna",analyseREPIN);
		treeFile=new File(inFolder+"/"+treeFile);
		if(!treeFile.exists()) {
			generateTree(treeFile);
		}
	}

	private void generateTree(File treeFile) {
		String filenames=generateFileNameString();
		String treeID=treeFile.getName().split("\\.")[0];
		File distFile=new File(inFolder+"/"+treeID+".dist");
		RunTreePrograms.runProgram("andi "+filenames, "", inFolder,distFile);
		RunTreePrograms.runProgram("clustDist "+distFile, "", inFolder, treeFile);
	}

	private String generateFileNameString() {
		StringBuffer sb=new StringBuffer();
		File[] files=inFolder.listFiles();
		for(int i=0;i<files.length;i++) {
			if(files[i].getName().endsWith("fas")||files[i].getName().endsWith("fna")) {
				sb.append(" "+files[i]);
			}
		}
		return sb.toString();
	}

	private String[] getREPtype() {
		ArrayList<String> list=new ArrayList<String>();
		for(int i=0;i<focalSeeds.length;i++) {
			list.add(i+"");
		}
		return list.toArray(new String[0]);
	}

	private String[] getFocalSeeds(String genome,int minRepFreq,int wl) {
		File fsg=new File(inFolder+"/"+genome);
		DetermineFocalSeeds dfs=new DetermineFocalSeeds(fsg,minRepFreq,wl);
		return dfs.getFocalSeeds();
	}

	public void print(File out) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			String[] genomes=results.keySet().toArray(new String[0]);
			for(int i=0;i<genomes.length;i++) {
				String[] seeds=results.get(genomes[i]).keySet().toArray(new String[0]);
				for(int j=0;j<seeds.length;j++) {
					bw.write(genomes[i].replace("_", "\t")+"\t"+seeds[j]+"\t"+results.get(genomes[i]).get(seeds[j])+"\n");
				}
			}
			bw.close();

		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private void calculateResults() {
		REPIN_RAYT_prox rrp=new REPIN_RAYT_prox();

		for(int i=0;i<genomes.size();i++) {
			String[] split=genomes.get(i).getAbsolutePath().split("\\/|\\.");
			for(int j=0;j<focalSeeds.length;j++) {
				String genomeID=split[split.length-2]+"_"+j;

				results.put(genomeID, new HashMap<String,Integer>());

				File outFolder=new File(inFolder+"/"+genomeID+"/");
				outFolder.mkdir();
				int wl=focalSeeds[j].length();
				
				REPINProperties rp=new REPINProperties(outFolder,genomeID,genomes.get(i),wl,numMuts,minFrac,null,focalSeeds[j],false,analyseREPIN);
				System.out.println("Write REPINs as artemis files for "+genomeID+"...");

				writeREPINArtemis(new File(outFolder+"/"+genomeID+"_largestCluster.ss"),j);
				writeREPINArtemis(new File(outFolder+"/"+genomeID+".ss"),j);
				File cluster;
				int k=0;
				while((cluster=new File(outFolder+"/"+genomeID+"_"+k+".ss")).exists()){
					writeREPINArtemis(cluster,k);
					k++;
				}
				System.out.println("Write RAYT locations "+genomeID+"...");

				ArrayList<Info> raytPos=writeRAYTLocation(genomes.get(i));
				int popsize=rp.getPopSize();
				results.get(genomeID).put(focalSeeds[j],popsize);
				System.out.println("REPIN RAYT proximity calculation for "+genomeID+"...");
				rrp.addRAYTREPINProximity(j, genomeID, outFolder, rp, raytPos);;
			}
		}
		rrp.writeStats(new File(inFolder+"/prox.stats"));

	}


	private ArrayList<Info> writeRAYTLocation(File genome) {
		String genomeID=genome.getName().split("\\.")[0];
		ArrayList<Info> RAYTLocations;
		if(legacyBlastPerlLocation!="") {
			RAYTLocations=BlastRAYTs.blastQuery(genome, queryRAYT, genomeFolder, e, "tblastn",legacyBlastPerlLocation);
		}else {
			RAYTLocations=BlastRAYTs.blastQuery(genome, queryRAYT, genomeFolder, e, "tblastn");
		}
		WriteArtemis.write(RAYTLocations, new File(genomeFolder+"/rayt_"+genomeID+".tab"));
		return RAYTLocations;
	}


	private void writeREPINArtemis(File in,int group) {
		if(in.exists()) {
			ArrayList<Fasta> fas=Fasta.readFasta(in);
			ArrayList<Info> pos=new ArrayList<Info>();
			for(int i=0;i<fas.size();i++) {
				String ident=fas.get(i).getIdent();
				pos.addAll(getPos(ident,"Gr_"+group));
			}
			String genomeID=in.getName().split("\\.")[0];
			WriteArtemis.write(pos, new File(in.getParent()+"/"+genomeID+".tab"));

		}
	}

	private ArrayList<Info> getPos(String ident,String inf){
		String split[]=ident.split("\\s+");
		ArrayList<Info> pos=new ArrayList<Info>();
		for(int i=1;i<split.length;i++) {
			String[] split2=split[i].split("_");
			int start=Integer.parseInt(split2[0]);
			int end=Integer.parseInt(split2[1]);
			pos.add(new Info(start,end,inf+" "+split[0]));
		}
		return pos;
	}


	private ArrayList<File> getFiles(File inFolder) {
		ArrayList<File> genomes=new ArrayList<File>();
		File[] all=inFolder.listFiles();
		for(int i=0;i<all.length;i++) {
			if(all[i].getAbsolutePath().endsWith(".fas")||all[i].getAbsolutePath().endsWith(".fna")) {
				genomes.add(all[i].getAbsoluteFile());
			}
		}
		return genomes;
	}

}
