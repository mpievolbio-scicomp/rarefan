package ecoliREP;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

import blastTools.PerformBlast;
import blastTools.ReadBlast;


import extragenicSequenceSimulation.GenerateExtragenicSequences;

import util.*;

//This program is supposed to print the numbers of REP elements found in each genome provided in a file including 
//the path to the genome, blast databases and genbank file (/home/frederic/auckland/outputFiles/REPAnalysisEcoli/REPComparison.in):
//NAME\tPathDataBase(genome fasta file)\tPathGenbank\thomologue(y\n)
//another input is the e-Value threshold for the results that are to be considered
//input: e-Value inputFile FileContainingNonOverlappingREPSequences outputFolder distance(length of flanking sequence considered 400bp standard) genome/notGenome(comparison to the whole genome or only to flanking sequences)
public class BlastREPs {
	public static void main(String args[]){
		double eValue=Double.parseDouble(args[0]);// eValue thresholds for when a REP sequence is a REP sequence
		File inputFile=new File(args[1]);
		File REPseqs=new File(args[2]);
		File out=new File(args[3]);
		int distance=Integer.parseInt(args[4]);
		boolean IsGenome=checkInput(args[5],"genome","flanking");
		boolean IsGenbank=checkInput(args[6],"genbank","noGenbank");
		//distance from REP to flanking sequence, in order to avoid the inclusion of REPs in the flanking sequence 
		int buffer=Integer.parseInt(args[7]);
		BlastREPs bR=new BlastREPs(out,IsGenome);
		bR.createBlast(inputFile,"blastn",eValue,REPseqs,false,false);
		bR.writeFrequency();
		bR.writeFlankingSequences(distance,false,IsGenbank,buffer);
		bR.blastFlanking(1e-35); // blast against genomes
		//blastFlanking(flanking,flanking,queries,out,names,1e-40,frequency,false);  //blast against flanking sequences
		
	}
	
	public static boolean checkInput(String input,String one,String two){
		
			if(input.equalsIgnoreCase(one))return true;
			else if(input.equalsIgnoreCase(two))return false;
				
			System.err.println("6th argument has to be either \"Genome\" or \"Flanking\".");
			System.exit(-1);
			return false;
	}
	File outputDir;
	boolean isGenome;
	public BlastREPs(File out,boolean IsGenome){
		outputDir=out;
		isGenome=IsGenome;
	}
	HashMap<String,HashMap<String,Integer>> distanceMatrix=new HashMap<String, HashMap<String,Integer>>();
	HashMap<String,Integer> frequency;
	//blasts the flanking sequences against another set of sequences (either genomes or again flanking sequences)
	public  void blastFlanking(double eValue){
		File temp=new File(outputDir+"/temp/");
		temp.deleteOnExit();
		if(!temp.exists()){
			temp.mkdir();
		}
		
		for(int i=0;i<names.size();i++){
			HashMap<String,File> outs=new HashMap<String, File>();
			for(int j=0;j<names.size();j++){
				if(i==j)continue;
				File out=new File(temp+"/"+j+".out");
				outs.put(names.get(j),out);
				out.deleteOnExit();
				File in=flanking.get(names.get(i));
				Input info=inputInfo.get(names.get(j));
				File database=isGenome?info.genome:flanking.get(names.get(j));
				PerformBlast.blast("","","blastn",eValue,out,in,database,false,false,true,false);		
			}
			analyse(outs,names.get(i));
		}
		writeDistanceMatrix();
		
	}
	//writes the number of REPs that are different between the files into a distance matrix file
	public  void writeDistanceMatrix(){
		try{
			String Genome=isGenome?"Genomes":"Flanking";
			BufferedWriter bw=new BufferedWriter(new FileWriter(new File(outputDir+"/REPDistanceMatrix"+Genome+".out")));
			bw.write("\t");
			for(int i=0;i<names.size();i++){
				bw.write(names.get(i)+"\t");
			}
			bw.write("\n");
			for(int i=0;i<names.size();i++){
				bw.write(names.get(i)+"\t");
				for(int j=0;j<names.size();j++){
					if(j==i){
						bw.write("0\t");
					}
					if(!distanceMatrix.get(names.get(i)).containsKey(names.get(j))){
						bw.write("0\t");
					}else{
						double value=distanceMatrix.get(names.get(i)).get(names.get(j))/(1.0*frequency.get(names.get(i)));
						bw.write(value+"\t");
					}
				}
				bw.write("\n");
			}
			
			bw.close();
			
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static TreeMap<Integer,String> sortHits(HashMap<String, String> hits,ArrayList<String> query){
		TreeMap<Integer,String> tm=new TreeMap<Integer, String>();
		for(int i=0;i<query.size();i++){
			String line=hits.containsKey(query.get(i))?query.get(i)+"\t"+hits.get(query.get(i)):query.get(i)+"\tmissing";
			String[] split=line.split("_");
			int pos=Integer.parseInt(split[1]);
			tm.put(pos,line);
		}
		return tm;
	}
	
	public  void writeMatches(HashMap<String,String> hits,File matches,ArrayList<String> query){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(matches));
			TreeMap<Integer,String> sorted=sortHits(hits,query);
			Iterator<Entry<Integer, String>> it=sorted.entrySet().iterator();
			
			while(it.hasNext()){
				Entry<Integer,String> e=it.next();
				String value=e.getValue();
				bw.write(value+"\n");
				
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	//creates the distances between files
	public  void analyse(HashMap<String,File> outs,String name){
		HashMap<String,String> hits=new HashMap<String, String>();
		String Genome=isGenome?"Genomes":"Flanking";
		File matchDir=new File(outputDir+"/matches/"+Genome+"/"+name+"/");
		if(!matchDir.exists()){
			matchDir.mkdirs();
		}
		HashMap<String,ArrayList<String>> notFound=new HashMap<String,ArrayList<String>>();
		distanceMatrix.put(name,new HashMap<String, Integer>());
		ArrayList<String> query=queries.get(name);
		for(int i=0;i<names.size();i++){
			
			if(names.get(i).equals(name))continue;
			ReadBlast rb=new ReadBlast(outs.get(names.get(i)));
			hits=new HashMap<String, String>();
			for(int j=0;j<rb.getQuery().size();j++){
				if(!hits.containsKey(rb.getQuery().get(j))) {
					String value="";
					if(isGenome){
						int start=rb.getStartDB().get(j)>rb.getEndDB().get(j)?rb.getEndDB().get(j):rb.getStartDB().get(j);
						int end=rb.getStartDB().get(j)<rb.getEndDB().get(j)?rb.getEndDB().get(j):rb.getStartDB().get(j);
						value=rb.getDatabase().get(j).replace("_", "")+"_"+start+"_"+end+"_-1\t"+rb.getEvalue().get(j);
					}else{
						value=rb.getDatabase().get(j)+"\t"+rb.getEvalue().get(j);
					}
					hits.put(rb.getQuery().get(j),value);
					
				}
			}
			
			writeMatches(hits,new File(matchDir+"/"+names.get(i)+".out"),query);
			for(int j=0;j<query.size();j++){
				if(!hits.containsKey(query.get(j))){
					if(!notFound.containsKey(query.get(j))){
						notFound.put(query.get(j), new ArrayList<String>());
					}
					if(!distanceMatrix.get(name).containsKey(names.get(i))){
						distanceMatrix.get(name).put(names.get(i),1);
					}else{
						distanceMatrix.get(name).put(names.get(i),distanceMatrix.get(name).get(names.get(i))+1);
					}
					notFound.get(query.get(j)).add(names.get(i));
				}
				
			}
		}
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(new File(outputDir+"/"+name+"_FlankingsNotFound"+Genome+".out")));
			Iterator<Entry<String,ArrayList<String>>> it=notFound.entrySet().iterator();
			while(it.hasNext()){
				Entry<String,ArrayList<String>> e=it.next();
				ArrayList<String> genomeNames=e.getValue();
				String queryName=e.getKey();
				bw.write(queryName+"\t"+genomeNames.size()+"\t");
				for(int i=0;i<genomeNames.size();i++){
					bw.write(genomeNames.get(i)+"\t");
				}
				bw.write("\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static ArrayList<Integer> getStartOrEnd(ArrayList<Integer> x,ArrayList<Integer> y,boolean starts){
		ArrayList<Integer> res=new ArrayList<Integer>();
		for(int i=0;i<x.size();i++){
			if(starts){
				res.add(x.get(i)>y.get(i)?y.get(i):x.get(i));
			}else{
				res.add(x.get(i)<y.get(i)?y.get(i):x.get(i));
			}
		}
		return res;
	}
	HashMap<String,ArrayList<String>> queries=new HashMap<String, ArrayList<String>>();
	HashMap<String,File> flanking=new HashMap<String, File>();
	

	public ArrayList<Info> makeInfoList(ArrayList<Integer> start,ArrayList<Integer> end){
		ArrayList<Info> al=new ArrayList<Info>();
		for(int i=0;i<start.size();i++){
			al.add(new Info(start.get(i),end.get(i),"Region"));
		}
		return al;
	}
	
	 

	
	public ArrayList<Info> getExSpaces(Input in,String genome,ArrayList<Integer> start,ArrayList<Integer> end,int distance,boolean genbank,int buffer){
		ArrayList<Info> REPextragenicSpaces=new ArrayList<Info>();
		if(genbank){
			GenerateExtragenicSequences ge=new GenerateExtragenicSequences(genome,in.genbank,true,distance);
			ArrayList<Info>  extras=ge.getIntervals();//clusterExSpace(SortArrayList.sort(ge.getIntervals()),distance,genome.length());
			IntervalTree<Info> exSpaceTree=new IntervalTree<Info>();
			for(int i=0;i<extras.size();i++){
				exSpaceTree.insert(extras.get(i));
			}
			for(int i=0;i<start.size();i++){
				ArrayList<Info> res=new ArrayList<Info>();
				exSpaceTree.search(new Info(start.get(i),end.get(i),""), res);
				if(res.size()!=1){
					System.err.println("Error: REP sequences overlap with "+res.size()+" extragenic spaces!\nAll are appended!");
					System.err.println(in.name);
					System.err.println("REP position: "+start.get(i)+" "+end.get(i));

				}
				for(int j=0;j<res.size();j++){
					exSpaceTree.remove(res.get(j));
				}
				REPextragenicSpaces.addAll(res);
			}
		}else{
			REPextragenicSpaces=cluster(makeInfoList(start, end),distance,genome.length(),buffer);
		}
		return REPextragenicSpaces;
	}
	
	//writes the REP flanking sequences into a file
	//int distance= clusters all exSpaces together that are not separated by a gene of at least 300bp
	//AllexSpaces -> if all extragenic spaces are to be used this is true other wise, the regions found in the fasta file are matched with extragenic spaces from a genbak file and only the exSpaces that overlap are used
	//genbank specifies only when AllExSpaces=false whether or not extragenic spaces are to be used for the definition of REPs
	//buffer specifies the distance from REP to flanking sequence (applies only to AllExSpaces=false and genbank=false)
	public  void writeFlankingSequences(int distance,boolean AllExSpaces,boolean genbank,int buffer){
		try{
			for(int i=0;i<names.size();i++){
				File flank=new File(outputDir+"/"+names.get(i)+".fas");
				Input in=inputInfo.get(names.get(i));
				//delete existing database files
				PerformBlast.deleteDatabases(flank);
				BufferedWriter bw=new BufferedWriter(new FileWriter(flank));
				String genome=ReadFasta.readFasta(in.genome).values().toArray(new StringBuilder[0])[0].toString();
				ReadBlast rb=new ReadBlast(in.blastoutREP);

				ArrayList<Integer> starts=getStartOrEnd(rb.getStartDB(),rb.getEndDB(),true);
				ArrayList<Integer> ends=getStartOrEnd(rb.getStartDB(),rb.getEndDB(),false);
				starts=SortArrayList.sort(starts);
				ends=SortArrayList.sort(ends);
				ArrayList<Info> infos=AllExSpaces?makeInfoList(starts,ends):getExSpaces(in,genome,starts,ends,distance,genbank,buffer);
				if(genbank)infos=exSpaceToFlanking(infos,distance,genome.length());	
				queries.put(names.get(i),new ArrayList<String>());
				for(int j=0;j<infos.size();j++){
					int start=infos.get(j).getStart();
					int end=infos.get(j).getEnd();
					int number=j/2;
					String seq=genome.substring(start,end);
					String name=infos.get(j).getInfo().replace("_", "")+"_"+start+"_"+end+"_"+number;
					queries.get(names.get(i)).add(name);
					bw.write(">"+name+"\n"+seq+"\n");
				}

				bw.close();
				flanking.put(names.get(i),flank);
			}
		}catch(IOException e ){
			e.printStackTrace();
		}
	}
	
		//writes regions from blast file
		public  void writeFlankingSequences(){
		try{
			for(int i=0;i<names.size();i++){
				File flank=new File(outputDir+"/"+names.get(i)+".fas");
				Input in=inputInfo.get(names.get(i));
				//delete existing database files
				PerformBlast.deleteDatabases(flank);
				BufferedWriter bw=new BufferedWriter(new FileWriter(flank));
				String genome=ReadFasta.readFasta(in.genome).values().toArray(new StringBuilder[0])[0].toString();
				ReadBlast rb=new ReadBlast(in.blastoutREP);

				ArrayList<Integer> starts=getStartOrEnd(rb.getStartDB(),rb.getEndDB(),true);
				ArrayList<Integer> ends=getStartOrEnd(rb.getStartDB(),rb.getEndDB(),false);
				starts=SortArrayList.sort(starts);
				ends=SortArrayList.sort(ends);
				ArrayList<Info> infos=makeInfoList(starts,ends);	
				queries.put(names.get(i),new ArrayList<String>());
				for(int j=0;j<infos.size();j++){
					int start=infos.get(j).getStart();
					int end=infos.get(j).getEnd();
					int number=j;
					String seq=genome.substring(start,end);
					String name=infos.get(j).getInfo().replace("_", "")+"_"+start+"_"+end+"_"+number;
					queries.get(names.get(i)).add(name);
					bw.write(">"+name+"\n"+seq+"\n");
				}

				bw.close();
				flanking.put(names.get(i),flank);
			}
		}catch(IOException e ){
			e.printStackTrace();
		}
	}
	
//	public static ArrayList<Info> clusterExSpace(ArrayList<Info> infos,int distance,int max){
//		ArrayList<Info> al=new ArrayList<Info>();
//		for(int i=0;i<infos.size();i++){
//			if(!(infos.get(i).getStart()-distance<0) ){
//				int currentI=i;
//				while(i<infos.size()-1 && (infos.get(i).getEnd()+distance>infos.get(i+1).getStart())){
//					i++;
//				}
//				if(!(infos.get(i).getEnd()+distance>=max)){
//				
//					al.add(new Info(infos.get(currentI).getStart(),infos.get(i).getEnd(),""));
//				}
//			}
//		}
//		
//		return al;
//	}


	
	public static ArrayList<Info> cluster(ArrayList<Info> infos,int distance,int max,int buffer){
		int start=0;
		ArrayList<Info> temp=new ArrayList<Info>();

		for(int i=0;i<=infos.size();i++){
			if(i==infos.size()){
				if(max-start>distance)temp.add(new Info(start,max,"test"));
				break;
			}
			if(infos.get(i).getStart()-start>distance){
				temp.add(new Info(start,infos.get(i).getStart(),"test"));
				
			}
			start=infos.get(i).getEnd();
		}
		ArrayList<Info> al=new ArrayList<Info>();

		for(int i=0;i<temp.size();i++){	
			if(!(temp.get(i).getEnd()-distance<0)&&!(temp.get(i).getStart()+distance>=max)){
				if(temp.get(i).getStart()!=0)
				al.add(new Info(temp.get(i).getStart()+buffer,temp.get(i).getStart()+distance+buffer,"REPright"));
				if(temp.get(i).getEnd()!=max)
				al.add(new Info(temp.get(i).getEnd()-distance-buffer,temp.get(i).getEnd()-buffer,"REPleft"));


			}
			
		}

		return al;
	}
	
	public static ArrayList<Info> exSpaceToFlanking(ArrayList<Info> infos,int distance,int max){
		ArrayList<Info> al=new ArrayList<Info>();
		for(int i=0;i<infos.size();i++){
			if(!(infos.get(i).getStart()-distance<0) &&!(infos.get(i).getEnd()+distance>=max)){
				
				
					al.add(new Info(infos.get(i).getStart()-distance,infos.get(i).getStart(),"REPleft"));
					al.add(new Info(infos.get(i).getEnd(),infos.get(i).getEnd()+distance,"REPright"));
				
			}
		}
		
		return al;
	}

	
	public  void writeFrequency(){
		frequency=new HashMap<String, Integer>();
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(outputDir+"/frequency.out"));
			for(int i=0;i<names.size();i++){
				Input in=inputInfo.get(names.get(i));
				ReadBlast rb=new ReadBlast(in.blastoutREP);
				bw.write(names.get(i)+"\t"+in.homologue+"\t"+rb.getQuery().size()+"\n");
				frequency.put(names.get(i), rb.getQuery().size());
			}
			bw.close();
		}catch(IOException e ){
			e.printStackTrace();
		}
	}
	

	public static boolean sequenceFile(File fasta,File out,int distance){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			String genome=ReadFasta.readFasta(fasta).values().toArray(new StringBuilder[0])[0].toString();
			for(int i=0;i<genome.length()-distance;i+=distance){
				int start=i;
				int end=i+distance;
				int length=distance;
	
				bw.write("extragenicSequence\t"+fasta.toString()+"\t-1\t"+length+"\t-1\t-1\t-1\t-1\t"+start+"\t"+end+"\t0\t10000\n");
			}
			
			bw.close();
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
		
		
	}
	public static boolean extragenicFile(File fasta,File genbank,File out,int distance){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			String genome=ReadFasta.readFasta(fasta).values().toArray(new StringBuilder[0])[0].toString();
			GenerateExtragenicSequences ge=new GenerateExtragenicSequences(genome,genbank,true,distance);
			ArrayList<Info> extras=ge.getIntervals();
			for(int i=0;i<extras.size();i++){
				int length=extras.get(i).getEnd()-extras.get(i).getStart();
	
				bw.write("extragenicSequence\t"+fasta.toString()+"\t-1\t"+length+"\t-1\t-1\t-1\t-1\t"+extras.get(i).getStart()+"\t"+extras.get(i).getEnd()+"\t0\t10000\n");
			}
			
			bw.close();
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
		
		
	}
	public  void createSequenceFiles(File input,int distance){
		File folder=new File(outputDir+"/temp/");
		folder.deleteOnExit();
		if(folder.mkdir()){
			try{
				BufferedReader br=new BufferedReader(new FileReader(input));
				String line="";
				while((line=br.readLine())!=null){
					if(line.startsWith("//"))continue;
					String[] split=line.split("\\s+");
					File genome=new File(split[1]);
					File genbank=new File(split[2]);
					File out=new File(folder+"/"+split[0]+".out");
					out.deleteOnExit();
					System.out.println("Calculate extragenic space for:"+split[0]);
					if(sequenceFile(genome,out,distance)){					
						File blastoutREP=out;
						String name=split[0];
						names.add(name);
						String homologue=split[4];
						inputInfo.put(name,new Input(blastoutREP,genome,genbank, homologue,name));
					}
				}
				br.close();
			}catch(IOException e){
				e.printStackTrace();
			}
			
		}else{
			System.err.print("Cannot create temp dir: "+outputDir+"/temp/");
			System.exit(1);
			
		}
	}
	public  void createExtragenicSequenceFiles(File input,int distance){
		File folder=new File(outputDir+"/temp/");
		folder.deleteOnExit();
		if(folder.mkdir()){
			try{
				BufferedReader br=new BufferedReader(new FileReader(input));
				String line="";
				while((line=br.readLine())!=null){
					if(line.startsWith("//"))continue;
					String[] split=line.split("\\s+");
					File genome=new File(split[1]);
					File genbank=new File(split[2]);
					File out=new File(folder+"/"+split[0]+".out");
					out.deleteOnExit();
					System.out.println("Calculate extragenic space for:"+split[0]);
					if(extragenicFile(genome,genbank,out,distance)){					
						File blastoutREP=out;
						String name=split[0];
						names.add(name);
						String homologue=split[4];
						inputInfo.put(name,new Input(blastoutREP,genome,genbank, homologue,name));
					}
				}
				br.close();
			}catch(IOException e){
				e.printStackTrace();
			}
			
		}else{
			System.err.print("Cannot create temp dir: "+outputDir+"/temp/");
			System.exit(1);
			
		}
	}
	
	public  void createBlast(File input,String program,double eValue,File REPsequences,boolean redo,boolean verbose){
		File folder=new File(outputDir+"/temp/");
		//folder.deleteOnExit();
		folder.mkdir();//if(folder.mkdir()){
			
			try{
				BufferedReader br=new BufferedReader(new FileReader(input));
				String line="";
				while((line=br.readLine())!=null){
					
					if(line.startsWith("#"))continue;
					String[] split=line.split("\\s+");
					File out=new File(folder+"/"+split[0]+".out");
					//out.deleteOnExit();
					File genome=new File(split[1]);
					boolean DNA=program.endsWith("blastn");
					if((!redo &&out.exists())||PerformBlast.blast("","",program,eValue,out,REPsequences,genome,verbose,false,DNA,false)){					
						File blastoutREP=out;
						String name=split[0];
						String homologue=split[4];
						File genbank=new File(split[2]);
						names.add(name);
						inputInfo.put(name,new Input(blastoutREP,genome,genbank,homologue,name));
					}
				}
				br.close();
			}catch(IOException e){
				e.printStackTrace();
			}
			
		//}else{
		//	System.err.print("Cannot create temp dir: "+outputDir+"/temp/");
		//	System.exit(1);
			
		//}
	}
	HashMap<String,InfoTree> repTree=null;
	
	public ArrayList<String> getNames(){
		return names;
	}
	
	public int countREPs(String name,Info interval){
		ArrayList<Info> al=new ArrayList<Info>();
		if(repTree==null){
			repTree=new HashMap<String, InfoTree>();
			for(int i=0;i<names.size();i++){
				Input in=inputInfo.get(names.get(i));
				ReadBlast rb=new ReadBlast(in.blastoutREP);
				ArrayList<Integer> start=getStartOrEnd(rb.getStartDB(), rb.getEndDB(), true);
				ArrayList<Integer> end=getStartOrEnd(rb.getStartDB(), rb.getEndDB(), false);
				ArrayList<Info> repInters=makeInfoList(start, end);
				InfoTree temp=new InfoTree();
				for(int j=0;j<repInters.size();j++){
					temp.insert(repInters.get(j));
				}
				repTree.put(names.get(i), temp);
			}
		}
		repTree.get(name).search(interval, al);
		return al.size();
	}
	
	ArrayList<String> names=new ArrayList<String>();
	public HashMap<String,Input> inputInfo=new HashMap<String, Input>();

	 
	
}
