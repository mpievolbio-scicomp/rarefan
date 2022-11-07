package ecoliREP;

import java.io.*;
import java.util.*;

import util.*;

//input is a REP output folder like /home/frederic/auckland/outputFiles/REPAnalysisEcoli/noGenbank/matches/Flanking and a 
//genome comparison output folder like /home/frederic/auckland/outputFiles/EcoliExtragenicSpace/Chop300bp/matches/Flanking
//and a file like REPComparison.in that specifies the species names...
//all comparison files are searched for irregularities (inversions, deletions, duplications) 
//each REP irregularity is then compared to the genome comparison the distance to the end of the irregularity in the genome comparison file is used as a quality value, 
//meaning the smaller the distance the more likely is a REP involvement for the recombination event

//!!!!!!!!!!!!!!
//if there is a potentially long deletion in the REP file that overlaps with two deleted regions in the genome file this deletion may only be due to two separate deletions of two or more REP regions
//At this point I will just add a quality value the Recombination class called intervalProportion, which is the proportion of genome that is deleted in comparison to the REP interval that is deleted   
//!!!!!!!!!!!!!!

public class AnalyseForRecombinations {
	
	public static void main(String args[]){
		File inFolderREP=new File(args[0]);
		File inFolderGenome=new File(args[1]);
		File config=new File(args[2]);
		File outFolder=new File(args[3]);
		//a subset of recombinations is recorded in the output folder under recombinationsBelow"threshold".out, only recombinations are recorded where the distance of the REP recombination event from the genome recombination event is less than threshold
		int threshold=Integer.parseInt(args[4]);
		
		ArrayList<Infos> infos=Infos.readInfos(config);
		writeRecombinations(inFolderGenome,inFolderREP,infos,outFolder,threshold);
	}
	
	public static void writeRecombinations(File inFolderGenome,File inFolderREP,ArrayList<Infos> infos,File outFolder,int threshold){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(outFolder+"/recombinationsBelow"+threshold+".out"));
			for(int i=0;i<infos.size();i++){
				String name1=infos.get(i).name;
				File currentFolder=new File(outFolder+"/"+name1+"/");
				if(!currentFolder.exists())currentFolder.mkdirs();
				ReadGenbank rgbk=new ReadGenbank(infos.get(i).genbank);
				for(int j=0;j<infos.size();j++){
					if(i==j)continue;
					String name2=infos.get(j).name;
					//System.out.println(name1+"\t"+name2);
					File inGenome=new File(inFolderGenome+"/"+name1+"/"+name2+".out");
					File inREP=new File(inFolderREP+"/"+name1+"/"+name2+".out");
					File out=new File(currentFolder+"/"+name2+".out");
					Matches genomeMatches=new Matches(inGenome);
					Matches REPMatches=new Matches(inREP);
					ArrayList<Recombination> recombs=getRecombinations(genomeMatches,REPMatches,rgbk.getInfoTree("CDS"));
					Recombination.write(recombs,out,threshold);
					bw.write(name1+"\t"+name2+"\n");
					bw.write(checkThreshold(recombs,threshold));

				}
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static String checkThreshold(ArrayList<Recombination> recs,int t){
		StringBuilder sb=new StringBuilder();
		for(int i=0;i<recs.size();i++){
			if(Math.abs(recs.get(i).qualityStart)<t && Math.abs(recs.get(i).qualityEnd)<t){
				sb.append(recs.get(i).toString()+"\n");
			}
		}
		return sb.toString();
	}
	
	public static ArrayList<Recombination> getRecombinations(Matches genome,Matches REP,InfoTree gbk){
		ArrayList<Recombination> recs=new ArrayList<Recombination>();
		System.out.println("Deletions:");
		recs.addAll(generateRecQuality(genome.findDeletions(),REP.findDeletions()));
		//System.out.println("Inversions");
		//recs.addAll(generateRecQuality(genome.findInversions(),REP.findInversions()));
		//System.out.println("Duplications");
		//recs.addAll(generateRecQuality(genome.findDuplications(),REP.findDuplications()));
		recs=addGeneInfo(recs,gbk);
		return recs;
	}
	
	public static ArrayList<Recombination> addGeneInfo(ArrayList<Recombination> recs,InfoTree gbkGenes){
		for(int i=0;i<recs.size();i++){
			ArrayList<Info> overlaps=new ArrayList<Info>();
			gbkGenes.search(recs.get(i).infoGenome, overlaps);
			recs.get(i).infoGenome.info+="_genes:";
			for(int j=0;j<overlaps.size();j++){
				recs.get(i).infoGenome.info+=overlaps.get(j).getGeneOrLocusInfo()+";";
				
			}
		}
		return recs;
	}
	
	public static ArrayList<Recombination> generateRecQuality(ArrayList<Info> recGenome,ArrayList<Info> recREP){
		ArrayList<Recombination> recs=new ArrayList<Recombination>();
		InfoTree itree=new InfoTree();
		for(int i=0;i<recGenome.size();i++){
			itree.insert(recGenome.get(i));
		}
		for(int i=0;i<recREP.size();i++){
			ArrayList<Info> overlaps=new ArrayList<Info>();
			itree.search(recREP.get(i), overlaps);
			if(overlaps.size()>=1){
				recs.add(mergeAndCheckRecombinations(overlaps,recREP.get(i)));
				
			}else {
				System.err.println("No overlapping recombinations found.");
			}
		}
		
		
		return recs;
	}
	
	private static Recombination mergeAndCheckRecombinations(ArrayList<Info> overlaps,Info recREP){
		int start=Integer.MAX_VALUE;
		int end=-1;
		int intervalSize=0;
		for(int i=0;i<overlaps.size();i++){
			intervalSize+=overlaps.get(i).getEnd()-overlaps.get(i).getStart();
			start=min(overlaps.get(i).getStart(),start);
			end=max(overlaps.get(i).getEnd(),end);
		}
		int qS=recREP.getStart()-start;
		int qE=end-recREP.getEnd();
		double intervalProportion=intervalSize/(1.0*recREP.getEnd()-recREP.getStart());
		return new Recombination(recREP,new Info(start,end,overlaps.get(0).info),qS,qE,intervalProportion);
		
	}
	
	static int max(int n1,int n2){
		return n1>n2?n1:n2;
	}
	static int min(int n1,int n2){
		return n1>n2?n2:n1;
	}
	public static class Recombination{
		
		public static void write(ArrayList<Recombination> recs,File out,int threshold){
			try{
				BufferedWriter bw=new BufferedWriter(new FileWriter(out));
				
//				for(int i=0;i<recs.size();i++){
//					bw.write(recs.get(i).toString()+"\n");
//				}
				bw.write(checkThreshold(recs, threshold));
				bw.close();
			}catch(IOException e){
				e.printStackTrace();
			}
		}
		
		Info infoREP;
		Info infoGenome;
		int qualityStart;
		int qualityEnd;
		double intervalProportion;
		public Recombination(Info InfoREP,Info InfoGenome,int QualityStart,int QualityEnd,double IntervalProportion){
			infoREP=InfoREP;
			infoGenome=InfoGenome;
			qualityEnd=QualityEnd;
			qualityStart=QualityStart;
			intervalProportion=IntervalProportion;
		}
		public String toString(){
			return infoREP+"\t"+infoGenome+"\t"+qualityStart+"\t"+qualityEnd+"\t"+intervalProportion;
		}
	}
	
}
