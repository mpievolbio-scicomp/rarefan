package yafMSignificance;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import pairwiseAlignment.NeedlemanWunsch;

import ecoliREP.*;
import util.*;
import util.IntervalTree.Flank;

public class REPInfluenceOnEvolution {
	public static void main(String args[]){
		File config=new File(args[0]);//like /home/frederic/auckland/outputFiles/mutationRate/flankingCandidates.in
		int minWordNumber=Integer.parseInt(args[1]);
		File outFolder=new File(args[2]);
		File genbank=new File(args[3]);
		String genome=Fasta.readFasta(new File(args[4])).get(0).getSequence();
		File wordFreqFolder=new File(args[5]);
		String name=args[6];
		int length=Integer.parseInt(args[7]);
		File mat=new File(args[8]);
		HashMap<Character,HashMap<Character,Integer>> matrix=NeedlemanWunsch.readSimilarityMatrix(mat);
		ArrayList<Info> intervals = GetFlankingGenes.getExSpacesAboveWordNumber(minWordNumber,
				genbank, genome, new File(wordFreqFolder + "/" + name + ".out"),length);
		InfoTree geneTree=GetFlankingGenes.makeGeneTree(genbank);
		ArrayList<Info> genes=changeNames(geneTree.parseTree(),name);
		geneTree=new InfoTree(genes);
		HashMap<String,HashMap<String, Info>> homologues=getGenes(outFolder, config,genes, genome);
		ArrayList<Info> genesFlanked=getFlankingGenes(geneTree, intervals, name);
		System.out.println("flanked genes:");
		produceResults(outFolder, config, genesFlanked, homologues, genome, name, "flanked",matrix,15,6.6);
		ArrayList<Info> genesNotFlanked=geneTree.delete(genesFlanked);
		System.out.println("unflanked genes:");
		produceResults(outFolder, config, genesNotFlanked, homologues, genome, name, "notFlanked",matrix,15,6.6);

	}
	public static ArrayList<Info> getFlankingGenes(InfoTree geneTree,
			ArrayList<Info> intervals, String name) {
		InfoTree flankingTree = new InfoTree();
		for (int i = 0; i < intervals.size(); i++) {
			Info interval = intervals.get(i);
			if (!GetFlankingGenes.overlaps(geneTree, interval, false)) {
				Flank<Info> f = geneTree.findFlank(interval);
				if (!GetFlankingGenes.overlaps(flankingTree, f.prev, false)) {

					flankingTree.insert(f.prev);
				}
				if (!GetFlankingGenes.overlaps(flankingTree, f.next, false)) {
					
					flankingTree.insert(f.next);
				}

			}
		}
		return flankingTree.parseTree();
	}
	public static ArrayList<Info> changeNames(ArrayList<Info> genes,String name){
		ArrayList<Info> newgenes=new ArrayList<Info>();
		for(int i=0;i<genes.size();i++){
			String complement=genes.get(i).info.contains("complement")?"_complement":"";
			Info inf=new Info(genes.get(i).getStart(),genes.get(i).getEnd(),name+"_"+genes.get(i).getGeneOrLocusInfo()+complement);
			newgenes.add(inf);
			
			
		}
		return newgenes;
	}
	
	public static void produceResults(File outFolder,File config,ArrayList<Info> genes,HashMap<String,HashMap<String,Info>> homologues,String genome,String name,String extension,HashMap<Character,HashMap<Character,Integer>> matrix,double dO,double dE){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(outFolder+"/results"+extension+".out"));
			NamesAndInput nai=new NamesAndInput(config);
			HashMap<String,Integer> hist=new HashMap<String, Integer>();
			ArrayList<String> names=nai.names;;
			HashMap<String,Input> input=nai.input;
			ArrayList<Fasta> genomes=getGenomes(input,names);
			for(int i=0;i<genes.size();i++){
				String key=genes.get(i).info+"_"+genes.get(i).getStart()+"_"+genes.get(i).getEnd();
				if(i==0)System.out.println(key);
				if(!present(homologues, key))continue;
				System.out.println(genes.get(i).info);
				ArrayList<Fasta> sequences=new ArrayList<Fasta>();
				sequences.add(Fasta.makeFasta(genes.get(i), genome,false ));
				bw.write(genes.get(i)+"\n");
				for(int j=0;j<names.size();j++){
					Info gene=homologues.get(names.get(j)).get(key);
					gene.info=gene.info.split("\\|")[0];
					sequences.add(Fasta.makeFasta(gene,genomes.get(j).getSequence(),false));
					bw.write(gene+"\n");
				}
				
				ArrayList<Double> identities=getIdentities(sequences,matrix,dO,dE);
				ArrayList<String> allNames=concatenate(name,names);
				writeIdentities(bw,identities,allNames);
				if(allNames.size()==3){
					String tt=treeTopo(allNames, identities);
					bw.write(tt+"\n");
					if(hist.containsKey(tt)){
						hist.put(tt, hist.get(tt)+1);
					}else{
						hist.put(tt,1);
					}
				}
				bw.write("\n");
			}
			bw.close();
			writeHist(hist,new File(outFolder+"/topo"+extension+".out"));
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static void writeHist(HashMap<String ,Integer> h,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			Iterator<Entry<String,Integer>> it=h.entrySet().iterator();
			while(it.hasNext()){
				Entry<String,Integer> e=it.next();
				bw.write(e.getKey()+" "+e.getValue()+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static void writeIdentities(BufferedWriter bw,ArrayList<Double> ident,ArrayList<String> names){
		int k=0;
		for(int i=0;i<names.size();i++){
			for(int j=i+1;j<names.size();j++){
				try{
					bw.write(names.get(i)+"<->"+names.get(j)+" "+ident.get(k)+"\n");	
				}catch(IOException e){
					e.printStackTrace();
				}
				k++;
			}
			
		}
	}
	
	public static ArrayList<Double> getIdentities(ArrayList<Fasta> sequences,HashMap<Character,HashMap<Character,Integer>> matrix,double dO,double dE){
		ArrayList<Double> identities=new ArrayList<Double>();
		for(int i=0;i<sequences.size();i++){
			for(int j=i+1;j<sequences.size();j++){
				NeedlemanWunsch nw=new NeedlemanWunsch(sequences.get(i).getSequence(),sequences.get(j).getSequence(),matrix,dO,dE);
				identities.add(nw.getPairwiseIdentity());
			}
		}
		return identities;
	}
	
	private static ArrayList<String> concatenate(String name,ArrayList<String> names){
		ArrayList<String> newNames=new ArrayList<String>();
		newNames.add(name);
		for(int i=0;i<names.size();i++){
			newNames.add(names.get(i));
		}
		return newNames;
	}
	
	public static ArrayList<Fasta> getGenomes(HashMap<String,Input> input,ArrayList<String> names){
		ArrayList<Fasta> genomes=new ArrayList<Fasta>();
		for(int i=0;i<names.size();i++){
			genomes.add(Fasta.readFasta(input.get(names.get(i)).genome).get(0));
		}
		return genomes;
	}
	//works only for three organisms
	public static String treeTopo(ArrayList<String> names,ArrayList<Double> identities){
		if(names.size()!=3||identities.size()!=3){
			System.err.println("Tree topo currently only works for three branches.");
			return "unknown";
		}
		return max(identities)==0?"(("+names.get(0)+","+names.get(1)+")"+names.get(2)+")":max(identities)==1?"(("+names.get(0)+","+names.get(2)+")"+names.get(1)+")":"(("+names.get(1)+","+names.get(2)+")"+names.get(0)+")";
	}
	
	public static int max(ArrayList<Double> items){
		double max=Double.MIN_VALUE;
		int index=-1;
		for(int i=0;i<items.size();i++){
			if(max<items.get(i)){
				max=items.get(i);
				index=i;
			}
		}
		return index;
	}
	
	public static boolean present(HashMap<String,HashMap<String,Info>> hm,String gene){
		Iterator<Entry<String,HashMap<String, Info>>> it=hm.entrySet().iterator();
		
		while(it.hasNext()){
			Entry<String,HashMap<String,Info>> e=it.next();
			if(!e.getValue().containsKey(gene)){
				return false;
			}
		}
		return true;
	}
	
	public static HashMap<String,HashMap<String,Info>> getGenes(File outFolder,File inputFile,ArrayList<Info> genes,String refGenome){
		
		HashMap<String,HashMap<String,Info>> blastResults=new HashMap<String,HashMap<String, Info>>();
		
		BlastREPs br = new BlastREPs(outFolder, true);
		File out=new File(outFolder+"/temp.out");
		Fasta.write(Fasta.makeFasta(genes, refGenome, false), out);
		br.createBlast(inputFile, "blastn", 1e-20, out,false,true);
		ArrayList<String> names = br.getNames();
		HashMap<String, Input> inputInfo = br.inputInfo;
		for (int i = 0; i < names.size(); i++) {
			Input in = inputInfo.get(names.get(i));
			InfoTree geneTree = GetFlankingGenes.makeGeneTree(in.genbank);
			ArrayList<Info> blastIntervals = GetFlankingGenes.blastToIntervals(in.blastoutREP);
			ArrayList<Info> geneOverlaps = GetFlankingGenes.getOverlappingGenes(geneTree,
					blastIntervals,false);
			blastResults.put(names.get(i),new HashMap<String, Info>());
			for(int j=0;j<geneOverlaps.size();j++){
				String key=geneOverlaps.get(j).info.split("\\|")[1];
				if(j==0)System.out.println(key);
				blastResults.get(names.get(i)).put(key, geneOverlaps.get(j));
			}
		}
		return blastResults;
	}
	
}
