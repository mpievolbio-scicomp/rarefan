package ecoliREP;

import java.io.*;
import java.util.*;

import util.Info;
import util.ReadGenbank;



//this program analyses match files created by BlastREPs.java
//the input is a tab separated file that includes information about the location of name, fasta file, genbank file and match file folder for the strains that are to be compared
//name\tfastaFile\tgenbankfile\tmatchFileFolder\tlist of strains that are to be compared with that strain
//second input is the output folder
//the output folder is going to contain a list of all recombinations including gene names, a matrix with the number of recombinations, a matrix with the numbers of REP and one that contains the quotient
//compares only pairs, like 1 with 2, 3 with 4...
//recombination event file:
//REP number\tstrain1\tstrain2\tstartPos\tendPos\tGeneInfo
public class GetREPdifference {
	public static void main(String args[]){
		System.out.println("Start...");
		File infos=new File(args[0]);
		File out=new File(args[1]);
		GetREPdifference grd=new GetREPdifference(infos);
		grd.write(out);
	}
	
	HashMap<String, HashMap<String,Integer>> numREPs=new HashMap<String, HashMap<String,Integer>>();
	HashMap<String,HashMap<String,Integer >> numRecs=new HashMap<String, HashMap<String,Integer>>();
	HashMap<String,HashMap<String,Double >> quo=new HashMap<String, HashMap<String,Double>>();

	ArrayList<PosAndMore> results=new ArrayList<PosAndMore>();
	ArrayList<Infos> infos=new ArrayList<Infos>();
	
	
	public GetREPdifference(File infoFile){
		infos=Infos.readInfos(infoFile);
		findMatches();
	}
	
	
	public  void findMatches(){
		results=new ArrayList<PosAndMore>();
		HashMap<String,ReadGenbank> genbanks=new HashMap<String, ReadGenbank>();

		for(int i=0;i<infos.size();i++){
			String name1=infos.get(i).name;
			for(int j=0;j<infos.size();j++){
				if(i==j)continue;
				String name2=infos.get(j).name;
				File match1with2=new File(infos.get(i).matchFolder+"/"+name2+".out");
				analyseMatchFile(match1with2);
				System.out.println(name1);
				if(!genbanks.containsKey(name1)){
					System.out.println("Reading "+infos.get(i).genbank);
					ReadGenbank rgb=new ReadGenbank(infos.get(i).genbank);
					genbanks.put(infos.get(i).name, rgb);
				}
					
				getGeneOverlaps(genbanks.get(name1));
				if(numREPs.containsKey(name1)){
					numREPs.get(name1).put(name2, reps);
				}else{
					HashMap<String,Integer> hm=new HashMap<String, Integer>();
					hm.put(name2,reps);
					numREPs.put(name1, hm);
				}
				if(numRecs.containsKey(name1)){
					numRecs.get(name1).put(name2, recs);
				}else{
					HashMap<String,Integer> hm=new HashMap<String, Integer>();
					hm.put(name2,recs);
					numRecs.put(name1, hm);
				}
				if(quo.containsKey(name1)){
					quo.get(name1).put(name2, recs/(reps*1.0));
				}else{
					HashMap<String,Double> hm=new HashMap<String, Double>();
					hm.put(name2,recs/(reps*1.0));
					quo.put(name1, hm);
				}
				results.add(new PosAndMore(annotation,name1,name2));
				
			}
		}
	}
	
	static void summarize(ArrayList<Info> inf,Info overlap){
		for(int i=0;i<inf.size();i++){
			overlap.append(inf.get(i));
		}
	}
	
	public  void getGeneOverlaps(ReadGenbank rgb){
		for(int i=0;i<annotation.size();i++){
			summarize(rgb.getFeatures(annotation.get(i)),annotation.get(i)); 
			
		}
	}
	
	public  void writeMatrix(File out,HashMap Matrix){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			
			bw.write("\t");
			for(int i=0;i<infos.size();i++){
				bw.write(infos.get(i).name+"\t");
			}
			bw.write("\n");
			for(int i=0;i<infos.size();i++){
				bw.write(infos.get(i).name+"\t");
				for(int j=0;j<infos.size();j++){
					if(j==i){
						bw.write("0\t");
						continue;
					}
					if(!((HashMap<String,HashMap<String,Object>>)Matrix).get(infos.get(i).name).containsKey(infos.get(j).name)){
						bw.write("0\t");
					}else{
						Object value=((HashMap<String,HashMap<String,Object>>)Matrix).get(infos.get(i).name).get(infos.get(j).name);
						
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
	public  void write(File out){
		try{
			File recMatrix=new File(out+"/recMatrix.out");
			File repMatrix=new File(out+"/repMatrix.out");
			File quoMatrix=new File(out+"/quoMatrix.out");
			File resList=new File(out+"/resList.out");
			BufferedWriter bw=new BufferedWriter(new FileWriter(resList));
			bw.write("Present in\tnot present in\tREP number\tstart\tend\n");
			for(int i=0;i<results.size();i++){
				for(int j=0;j<results.get(i).pos.size();j++){
					bw.write(results.get(i).name1+"\t"+results.get(i).name2+"\t"+results.get(i).pos.get(j).repNumber+"\t"+results.get(i).pos.get(j).getStart()+"\t"+results.get(i).pos.get(j).getEnd()+"\t"+results.get(i).pos.get(j).getInfo()+"\n");
				}
			}
			bw.close();
			writeMatrix(recMatrix, numRecs);
			writeMatrix(repMatrix, numREPs);
			writeMatrix(quoMatrix, quo);

		}catch(IOException e){
			e.printStackTrace();
		}
	}
	ArrayList<Info> annotation=new ArrayList<Info>();
	int reps=0;
	int recs=0;
	public  void analyseMatchFile(File match){
		try{
			annotation=new ArrayList<Info>();
			BufferedReader br=new BufferedReader(new FileReader(match));
			String line="";
			String prev[]=null;
			reps=0;
			recs=0;
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+|_");
				if(prev!=null){
					if(prev[0].contains("left")&&prev[4].contains("missing")&&!split[4].contains("missing")){
						annotation.add(new Info(Integer.parseInt(prev[1]),Integer.parseInt(prev[2]),"",reps/2));
						recs++;
					}else if(split[0].contains("right")&&split[4].contains("missing")&&!prev[4].contains("missing")){
						annotation.add(new Info(Integer.parseInt(split[1]),Integer.parseInt(split[2]),"",reps/2));
						recs++;
					}
				}
				prev=split;
				reps++;
			}
		}catch(IOException e){
			e.printStackTrace();
		}
		
		
	}
	

	

	

	
	public static class PosAndMore{
		public PosAndMore(ArrayList<Info> Pos,String n1,String n2){
			pos=Pos;
			name1=n1;
			name2=n2;
		}
		ArrayList<Info> pos;
		String name1;
		String name2;
		
	}
	
}
