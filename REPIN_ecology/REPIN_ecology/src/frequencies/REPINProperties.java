package frequencies;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.*;

import REPINpopulations.REPINposition;
import frequencies.util.*;
import seedAnalysis.WriteSeedSequences;
import util.Fasta;
import util.phylogenetics.RunTreePrograms;

public class REPINProperties {
	//File mclPath=new File("/home/bertels/Programs/mcl-14-137//src/shmcl/mcl");
    //File mclPath=new File("/Users/bertels/Programs/mcl-12-068/src/shmcl/mcl");
    String mclPath="mcl";
	File wordFrequencies;
	String genomeID;
	File seedSequence;
	HashMap<String,Double> mutFreqs;
	File mutFreqFile;
	File nodes;
	File slopes;
	File fas;
	File mutOccFile;
	int wordlength;
	static String seedExt=".ss";
	String mutFreqExt=".mf";
	String mutOccExt=".mo";
	File wfr;
	MaxWord maxWord;
	File outFolder;
	String simNetExt=".nw";
	String ddExt=".dd";
	String nodesExt=".nodes";
	String nodeDistExt=".nodeDist";
	String mutRateExt=".mr";
	File degreeDist;
	File simNet;
	int numSeqs;
	double propNonHub;
	//only simulation
	String regime;
	NetworkProperties nw;
	File largestCluster;
	File mutHist;
	File mutFreqTimeOut;
	File mutationClassHist;
	boolean analyseREPIN=true;
	File fitnessOut;
	int mutclasses;
	int numDifferencesToCluster=2;
	String word=null;
	int popsize;
	boolean needsToContainWord=false;
	HashMap<String,ArrayList<REPINposition>> repinPositions=new HashMap<String,ArrayList<REPINposition>>();
	HashMap<String,ArrayList<REPINposition>> largestClusterRepinPositions;

	public String getRegime(){
		return regime;
	}
	
	public HashMap<String,ArrayList<REPINposition>> getREPINPositions(){
		return repinPositions;
	}
	
	public HashMap<String,ArrayList<REPINposition>> getLargestClusterREPINPositions(){
		if(largestClusterRepinPositions==null) {
			largestClusterRepinPositions=getSeedSequencePositions();
		}
		return largestClusterRepinPositions;

	}
	
	private HashMap<String,ArrayList<REPINposition>> getSeedSequencePositions(){
		ArrayList<Fasta> fas=Fasta.readFasta(seedSequence);
		HashMap<String,ArrayList<REPINposition>> pos=new HashMap<String,ArrayList<REPINposition>>();
		for(int i=0;i<fas.size();i++) {
			String line=fas.get(i).getIdent();
			String[] split=line.split("\\s+");
			String seq=fas.get(i).getSequence();
			pos.put(seq,new ArrayList<REPINposition>());
			for(int j=1;j<split.length;j++) {
				String posSplit[]=split[j].split("_");
				int id=Integer.parseInt(posSplit[0]);
				int start=Integer.parseInt(posSplit[1]);
				int end=Integer.parseInt(posSplit[2]);
				pos.get(seq).add(new REPINposition(start, end, id));
			}
		}
		return pos;

	}
	
	

	public REPINProperties(File outFolder,String genomeID,File fas,int wordlength,int numMuts,double minFrac,File mutRate,String word,boolean needsToContainWord,boolean analyseREPIN){
		this.word=word;
		this.needsToContainWord=needsToContainWord;
		this.outFolder=outFolder;
		this.genomeID=genomeID;
		this.fas=fas;
		this.analyseREPIN=analyseREPIN;
		this.wordlength=wordlength;
		seedSequence=new File(outFolder+"/"+genomeID+seedExt);
		System.out.println("Processing genome "+genomeID);
		System.out.println("Calculate word frequencies...");
		maxWord=word==null?new MaxWord(fas,wordlength,outFolder,genomeID):new MaxWord(fas,word,wordlength,outFolder,genomeID);
		if(maxWord.getFrequency()>11){
			System.out.println("Write seed sequences...");
			//NEEDS TO BE CONVERTED TO REPINS!!!!
			writeSeedSequence(numMuts);
			//if(regime!=null&&regime.equals("cluster")){
			System.out.println("Determine largest sequence cluster...");
			determineLargestSequenceCluster();
			//}
			if(analyseREPIN)writeAllSequenceClusters();
			//System.out.println("Write mutation frequencies...");
			//writeMutationFrequencies();
			//System.out.println("Calculate mutation rates...");
			numSeqs=getNumSeqs(seedSequence);
			if(numSeqs>0){
			    System.out.println("Determine mutation histogram...");
			    calculateHistogram();
				//mutationRate=new MutationRate(mutFreqs, minFrac, outFolder, genomeID,wordlength,numSeqs);
				//mutationRate.write(mutRate,regime);
				System.out.println("Calculate Similarity Network...");
				calculateSimilarityNetwork();
				calculatePropNonMainHub();
			}
			System.out.println("Done with REPIN properties.");
		}
	}


	private void calculateHistogram(){
		mutHist=new File(outFolder+"/"+genomeID+".hist");
		mutationClassHist=new File(outFolder+"/"+genomeID+".mch");

		//if(!mutHist.exists()){
			CalculateDistanceHist cdh=new CalculateDistanceHist();
			System.out.println(seedSequence);
			cdh.write(seedSequence, mutHist);
			//cdh.writeMCH(seedSequence, mutationClassHist);
		//}
	}


	private void writeAllSequenceClusters(){
		File mclout=new File(outFolder+"/"+genomeID+".mcl");
		File seeds=new File(outFolder+"/"+genomeID+seedExt);
		ArrayList<ArrayList<String>> seqSelection=getAllClusters(mclout);
		for(int i=0;i<seqSelection.size();i++) {
			File cluster=new File(outFolder+"/"+genomeID+"_"+i+".ss");
			replaceSeedSequences(cluster,seqSelection.get(i),seeds);

		}
	}


	//will have to overwrite "File seedSequence"
	public void determineLargestSequenceCluster(){
		simNet=new File(outFolder+"/"+genomeID+"_allSeed.nw");
		largestCluster=new File(outFolder+"/"+genomeID+"_largestCluster.nodes");
		nodes=new File(outFolder+"/"+genomeID+nodesExt);
		File mclout=new File(outFolder+"/"+genomeID+".mcl");
		File newSeedSequences=new File(outFolder+"/"+genomeID+"_largestCluster"+seedExt);
		//if(!simNet.exists()){
			//simNet=new File(outFolder+"/"+genomeID+simNetExt);

			SimilarityNetwork sn=new SimilarityNetwork(seedSequence);
			sn.writeCytoscapeInput(simNet);
			sn.writeNodes(nodes);

			if(analyseREPIN){
                //>>> DEBUG
                System.out.println("Working Directory = " + System.getProperty("user.dir"));
                String mcl_command =mclPath+" "+simNet+" "+" -I 1.2 --abc -o "+mclout;
                System.out.println("mcl_command = " + mcl_command);
                //<<< DEBUG
				RunTreePrograms.runProgram(
                        mcl_command,
                        "",
                        new File(System.getProperty("user.dir")));
				int maxNode=getMaxNode(nodes);
				ArrayList<String> seqSelection=getLargestCluster(mclout,maxNode,true);
				seqSelection=addSequences(seqSelection,numDifferencesToCluster);//genomeID.equals("DC3000")||genomeID.equals("putGB1")?addSequences(seqSelection,numDifferencesToCluster+1):addSequences(seqSelection,numDifferencesToCluster);
				replaceSeedSequences(newSeedSequences,seqSelection,seedSequence);
				writeNodes(largestCluster,seqSelection);
				seedSequence=newSeedSequences;
			}else{
				sn.writeNodes(largestCluster);
				try{
					Files.copy(seedSequence.toPath(), newSeedSequences.toPath(),StandardCopyOption.REPLACE_EXISTING);
				}catch(IOException e){
					e.printStackTrace();
					System.exit(-1);
				}
			}
		//}
			if(analyseREPIN){
				seedSequence=newSeedSequences;

			}
			popsize=getNumREPINs(nodes);
	}

	public int getPopSize() {
		return popsize;
	}

	private String makeAs() {
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<wordlength;i++) {
			sb.append('A');
		}
		return sb.toString();
	}
	private int getMaxNode(File in){
		int maxNode=-1;
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				int nodefreq=Integer.parseInt(split[1]);


				if(maxNode<nodefreq && !split[0].contains(makeAs())){
					maxNode=nodefreq;
				}
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return maxNode;
	}

	private ArrayList<String> addSequences(ArrayList<String> seqSelection,int muts){
		HashSet<String> seqSelectionHash=new HashSet<String>(seqSelection);
		ArrayList<String> newSelection=new ArrayList<String>();
		ArrayList<Fasta> allSeqs=Fasta.readFasta(seedSequence);
		for(int i=0;i<allSeqs.size();i++){
			String current=allSeqs.get(i).getSequence();

			for(int j=0;j<seqSelection.size();j++){
				String ss=seqSelection.get(j).split("_")[0];
				//System.out.println(ss+" "+current);
				int diff=SimilarityNetwork.getNumDifferences(ss, current);
				if(diff<=muts){
					int occ=Integer.parseInt(allSeqs.get(i).getIdent().substring(10).split("\\s+")[0]);
					if(occ<2||seqSelectionHash.contains(current+"_"+occ))newSelection.add(current+"_"+occ);
					break;
				}
			}
		}
		return newSelection;
	}

	private void writeNodes(File out,ArrayList<String> nodes){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<nodes.size();i++){
				bw.write(nodes.get(i)+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}


	private ArrayList<String> getLargestCluster(File in,int maxNode,boolean useMaxNode){
		ArrayList<String> ids=new ArrayList<String>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				if(word!=null && needsToContainWord){
					if(line.contains(word)&& !line.contains(makeAs())){
						String split[]=line.split("\\s+");
						for(int i=0;i<split.length;i++){
							ids.add(split[i]);
						}
						break;

					}
				}else if((line.contains("_"+maxNode+"")||!useMaxNode)&& !line.contains(makeAs())){
					String split[]=line.split("\\s+");
					for(int i=0;i<split.length;i++){
						ids.add(split[i]);
					}
					break;
				}
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return ids;
	}

	private ArrayList<ArrayList<String>> getAllClusters(File in){
		ArrayList<ArrayList<String>> allClusters=new ArrayList<ArrayList<String>>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String split[]=line.split("\\s+");
				ArrayList<String> ids=new ArrayList<String>();

				for(int i=0;i<split.length;i++){
					ids.add(split[i]);
				}
				allClusters.add(ids);

			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return allClusters;
	}

	private HashMap<String,Fasta> getHash(ArrayList<Fasta> fas){
		HashMap<String,Fasta> hm=new HashMap<String,Fasta>();
		for(int i=0;i<fas.size();i++){
			String id=fas.get(i).getIdent();
			id=id.replaceAll("[A-Za-z]", "").split("\\s+")[0];
			int count=Integer.parseInt(id);
			id=fas.get(i).getSequence()+"_"+count;
			hm.put(id,new Fasta(fas.get(i).getIdent(), fas.get(i).getSequence()));
		}
		return hm;
	}

	private void replaceSeedSequences(File out,ArrayList<String> ids,File seeds){
		HashMap<String,Fasta> fas=getHash(Fasta.readFasta(seeds));

		ArrayList<Fasta> largestCluster=new ArrayList<Fasta>();
		for(int i=0;i<ids.size();i++){

			largestCluster.add(fas.get(ids.get(i)));
		}
		Fasta.write(largestCluster, out);
	}



	public REPINProperties(File outFolder,REPINProperties real,String regime,double minFrac,File mutRate,File fitnessOut,int mutclasses){
		this.mutclasses=mutclasses;
		this.fitnessOut=fitnessOut;
		this.outFolder=outFolder;
		this.genomeID=real.genomeID;
		this.wordlength=real.wordlength;
		this.regime=regime;
		seedSequence=new File(outFolder+"/"+genomeID+seedExt);
//		if(real.mutationRate!=null){
//
//			System.out.println("Simulation for genome "+genomeID);
//			simulateSeqEvol(real);
//			System.out.println(real);
//
//			writeMutationFrequencies();
//			mutationRate=new MutationRate(mutFreqs,minFrac,outFolder,genomeID,wordlength,numSeqs);
//			mutationRate.write(mutRate,regime);
			System.out.println("Calculate Similarity Network...");
			calculateSimilarityNetwork();
			System.out.println("Determine mutation histogram...");
			calculateHistogram();

	}


	private void calculatePropNonMainHub(){
		double maxSeed=getNumMaxSeqs(seedSequence);

		propNonHub= (1.0*maxSeed)/numSeqs;
	}

	public NetworkProperties getNetworkProperties(){
		return nw;
	}

//	private void simulateSeqEvol(REPINProperties real){
//			int max=real.numSeqs;
//			mutationRate=real.mutationRate;
//			mutOccFile=new File(outFolder+"/"+genomeID+mutOccExt);
//			mutFreqTimeOut=new File(outFolder+"/"+genomeID+".mutFreqTime");
//			//double mr=mutationRate.getMutationRate();
//			String id=seedSequence.getName().split("\\.")[0];
//			double GC=calculateGC(real.fas);
//			SimulateREPEvolution sse=new SimulateREPEvolution( max, regime,real.propNonHub,real.outFolder,id,GC,fitnessOut,mutclasses);
//			sse.printMutationOccurrenceTimes(mutOccFile);
//			sse.writeFasta(seedSequence);
//			sse.printMutFreqTime(mutFreqTimeOut);
//	}



//	private double calculateGC(File fas){
//		ArrayList<Fasta> seqs=Fasta.readFasta(fas);
//		int count=0;
//		int GC=0;
//		for(int i=0;i<seqs.size();i++){
//			String seq=seqs.get(i).getSequence();
//			seq=seq.toUpperCase();
//			for(int j=0;j<seq.length();j++){
//				count++;
//				if(seq.charAt(j)=='G'||seq.charAt(j)=='C'){
//					GC++;
//				}
//			}
//		}
//		return (GC*1.0)/count;
//	}

	static int getNumMaxSeqs(File in){
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		int max=0;
		for(int i=0;i<fas.size();i++){
			int curr=RAYTfrequencyGraph.getFreq(fas.get(i).getIdent());
			if(curr>max){
				max=curr;
			}
		}
		return max;
	}

	public static String getMaxSequence(File in){
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		int max=0;
		String seq="";
		for(int i=0;i<fas.size();i++){
			int curr=RAYTfrequencyGraph.getFreq(fas.get(i).getIdent());
			if(curr>max){
				max=curr;
				seq=fas.get(i).getSequence();
			}
		}
		return seq;
	}

	static int getNumSeqs(File in){
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		int sum=0;
		for(int i=0;i<fas.size();i++){
			sum+=RAYTfrequencyGraph.getFreq(fas.get(i).getIdent());
		}
		return sum;
	}

	int getNumREPINs(File in){
		int sum=0;

		try {
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			String as=makeAs();
			while((line=br.readLine())!=null) {
				String[] split=line.split("_|\\s+");
				int num=Integer.parseInt(split[1]);
				String seq=split[0];
				if(!seq.contains(as)||!analyseREPIN) {
					sum+=num;
				}

			}
			br.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return sum;
	}

	private void calculateSimilarityNetwork(){
		degreeDist=new File(outFolder+"/"+genomeID+ddExt);
		SimilarityNetwork sn=new SimilarityNetwork(seedSequence);
		if(simNet==null){
			simNet=new File(outFolder+"/"+genomeID+simNetExt);
			nodes=new File(outFolder+"/"+genomeID+nodesExt);
			sn.writeCytoscapeInput(simNet);
		    sn.writeNodes(nodes);
		}
		CalculateNetworkProperties cnp=new CalculateNetworkProperties(sn.getNetwork(), sn.getNodes(),genomeID,outFolder,regime);
		nw=cnp.getNetworkProperties();

		nw.writeDegreeDist(degreeDist);
	}


	private void writeSeedSequence(int numMuts){
		ArrayList<String> rawWords=new ArrayList<String>();
		rawWords.add(maxWord.getMaxWord());
		File seedSequenceREP=new File(seedSequence+".REP");
		if(!seedSequenceREP.exists())WriteSeedSequences.writeSeedSequencesConnected(rawWords, fas,  numMuts, seedSequenceREP);
		//if(!seedSequence.exists()) {//&& !ecoli){
			ConvertToREPIN ctr=new ConvertToREPIN(seedSequenceREP,fas,130,maxWord.getMaxWord());
			ctr.write(seedSequence);
			repinPositions=ctr.getPositions();

		//}else {
		//	if(ecoli) {
		//		util.FileHandler.copy(seedSequenceREP, seedSequence);
		//	}
		//}

	}
//	private void writeMutationFrequencies(){
//		mutFreqs=new HashMap<String,Double>();
//		mutFreqFile=new File(outFolder+"/"+genomeID+mutFreqExt);
//
//		RAYTfrequencyGraph rfg=new RAYTfrequencyGraph(seedSequence);
//		if(!mutFreqFile.exists()){
//			rfg.writeMutationFrequencies(mutFreqFile);
//		}
//		System.out.println(mutFreqFile);
//		mutFreqs=rfg.getMutationFrequencies();
//	}



}
