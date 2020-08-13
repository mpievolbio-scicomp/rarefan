package yafMSignificance;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import pairwiseAlignment.NeedlemanWunsch;

import util.Fasta;
import util.Info;
import util.InfoTree;
import util.ReadFasta;
import ecoliREP.BlastREPs;
import ecoliREP.Input;

public class GetFlankingGenesBothSides {
	public static void main(String args[]) {
		File genbank = new File(args[0]);
		HashMap<String, StringBuilder> genomeFasta = ReadFasta
				.readFasta(new File(args[1]));
		String genome = genomeFasta.values().toArray(new StringBuilder[0])[0]
				.toString();
		String name = args[2];
		File out = new File(args[3]);
		File inputFile = new File(args[4]);
		int length = Integer.parseInt(args[5]);
		File wordFreqFolder = new File(args[6]);
		int minWordNumber = Integer.parseInt(args[7]);
		File mat=new File(args[8]);
		HashMap<Character,HashMap<Character,Integer>> matrix=NeedlemanWunsch.readSimilarityMatrix(mat);
		
		ArrayList<Info> intervals = GetFlankingGenes.getExSpacesAboveWordNumber(minWordNumber,
				genbank, genome, new File(wordFreqFolder + "/" + name + ".out"),length);
		InfoTree geneTree=GetFlankingGenes.makeGeneTree(genbank);
		ArrayList<Info> genesBoth = GetFlankingGenes.getFlankingGenesBothSides( geneTree,intervals, name);
		File flankingGenesBoth = new File(out + "/" + name + "_flankedGenesBoth.fas");
		ArrayList<Fasta> bothSidesFas=Fasta.makeFasta(genesBoth, genome, false);
		Fasta.write(bothSidesFas, flankingGenesBoth);
		
		ArrayList<Info> genesLeft = GetFlankingGenes.getFlankingGenesLeftSide(geneTree, genesBoth, name);
		File flankingGenesLeft = new File(out + "/" + name + "_flankedGenesLeft.fas");
		ArrayList<Fasta> LeftSideFas=Fasta.makeFasta(genesLeft, genome, false);
		Fasta.write(LeftSideFas, flankingGenesLeft);
		
		ArrayList<Info> genesRight = GetFlankingGenes.getFlankingGenesRightSide(geneTree, genesBoth, name);
		File flankingGenesRight = new File(out + "/" + name + "_flankedGenesRight.fas");
		ArrayList<Fasta> RightSideFas=Fasta.makeFasta(genesRight, genome, false);
		Fasta.write(RightSideFas, flankingGenesRight);
		
		for(int i=0;i<bothSidesFas.size();i++){
			File tempFile=new File(out+"/temp.fas");
			tempFile.deleteOnExit();
			writeFasta(LeftSideFas.get(i), tempFile);
			ArrayList<Fasta> lefties=getGenes(out,inputFile,tempFile,genome);
			for(int j=0;j<lefties.size();j++){
				String a=LeftSideFas.get(i).getSequence();
				String b=lefties.get(j).getSequence();
				System.out.println("Left: "+LeftSideFas.get(i).getIdent()+"<->"+lefties.get(j).getIdent());
				NeedlemanWunsch nw=new NeedlemanWunsch(a,b,matrix,15,6.6);
				System.out.println(nw.getPairwiseIdentity());
				for(int k=j+1;k<lefties.size();k++){
					a=lefties.get(k).getSequence();
					b=lefties.get(j).getSequence();
					System.out.println("Left: "+lefties.get(k).getIdent()+"<->"+lefties.get(j).getIdent());
					nw=new NeedlemanWunsch(a,b,matrix,15,6.6);
					System.out.println(nw.getPairwiseIdentity());
				
				}

			}
			System.out.println("_______________________");
			writeFasta(bothSidesFas.get(i), tempFile);
			ArrayList<Fasta> both=getGenes(out,inputFile,tempFile,genome);
			for(int j=0;j<both.size();j++){
				String a=bothSidesFas.get(i).getSequence();
				String b=both.get(j).getSequence();
				System.out.println("Both: "+bothSidesFas.get(i).getIdent()+"<->"+both.get(j).getIdent());
				NeedlemanWunsch nw=new NeedlemanWunsch(a,b,matrix,15,6.6);
				System.out.println(nw.getPairwiseIdentity());
				for(int k=j+1;k<both.size();k++){
					a=both.get(k).getSequence();
					b=both.get(j).getSequence();
					System.out.println("Both: "+both.get(k).getIdent()+"<->"+both.get(j).getIdent());
					nw=new NeedlemanWunsch(a,b,matrix,15,6.6);
					System.out.println(nw.getPairwiseIdentity());				
				}
			}
			System.out.println("_______________________");

			writeFasta(RightSideFas.get(i), tempFile);
			ArrayList<Fasta> right=getGenes(out,inputFile,tempFile,genome);
			for(int j=0;j<right.size();j++){
				String a=RightSideFas.get(i).getSequence();
				String b=right.get(j).getSequence();
				System.out.println("Right: "+RightSideFas.get(i).getIdent()+"<->"+right.get(j).getIdent());
				NeedlemanWunsch nw=new NeedlemanWunsch(a,b,matrix,15,6.6);
				System.out.println(nw.getPairwiseIdentity());				
				for(int k=j+1;k<right.size();k++){
					a=right.get(k).getSequence();
					b=right.get(j).getSequence();
					System.out.println("Right: "+right.get(k).getIdent()+"<->"+right.get(j).getIdent());
					nw=new NeedlemanWunsch(a,b,matrix,15,6.6);
					System.out.println(nw.getPairwiseIdentity());				
				}
			}
			System.out.println("\n");
		}
		
		
//		writeFastaGenes(out,inputFile,flankingGenesBoth,"Both",genome);
//		ArrayList<Info> genesLeft = GetFlankingGenes.getFlankingGenesLeftSide(geneTree, genesBoth, name);
//		File flankingGenesLeft = new File(out + "/" + name + "_flankedGenesLeft.fas");
//		Fasta.write(Fasta.makeFasta(genesLeft, genome, true), flankingGenesLeft);
//		writeFastaGenes(out,inputFile,flankingGenesLeft,"Left",genome);		
//		ArrayList<Info> genesRight = GetFlankingGenes.getFlankingGenesRightSide(geneTree, genesBoth, name);
//		File flankingGenesRight = new File(out + "/" + name + "_flankedGenesRight.fas");
//		Fasta.write(Fasta.makeFasta(genesRight, genome, true), flankingGenesRight);
//		writeFastaGenes(out,inputFile,flankingGenesRight,"Right",genome);	

	}
	
	public static void writeFasta(Fasta seq,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write(seq.getIdent()+"\n"+seq.getSequence()+"\n");
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			
		}
	}
	public static ArrayList<Fasta> getGenes(File outFolder,File inputFile,File gene,String refGenome){
		BlastREPs br = new BlastREPs(outFolder, true);
		br.createBlast(inputFile, "blastn", 1e-20, gene,true,true);
		ArrayList<String> names = br.getNames();
		HashMap<String, Input> inputInfo = br.inputInfo;
		ArrayList<Fasta> genes=new ArrayList<Fasta>();
		for (int i = 0; i < names.size(); i++) {
			//System.out.println(names.get(i));
			// find all homologous genes, intervals
			Input in = inputInfo.get(names.get(i));
			HashMap<String, StringBuilder> genomeFastaLocal = ReadFasta
					.readFasta(in.genome);
			String genomeLocal = genomeFastaLocal.values().toArray(
					new StringBuilder[0])[0].toString();
			InfoTree geneTree = GetFlankingGenes.makeGeneTree(in.genbank);
			ArrayList<Info> blastIntervals = GetFlankingGenes.blastToIntervals(in.blastoutREP);
			ArrayList<Info> geneOverlaps = GetFlankingGenes.getOverlappingGenes(geneTree,
					blastIntervals,false);
			for(int j=0;j<geneOverlaps.size();j++){
				ArrayList<Info> localGeneInfo=new ArrayList<Info>();
				String[] split=geneOverlaps.get(j).info.split("\\|");
				String info=names.get(i)+"_"+split[0];
				localGeneInfo.add(new Info(geneOverlaps.get(j).getStart(),geneOverlaps.get(j).getEnd(),info));
				genes.addAll(Fasta.makeFasta(localGeneInfo, genomeLocal, false));

			}
		}
		return genes;
	}

	
	
	
	
	public static void writeFastaGenes(File outFolder,File inputFile,File flankingGenes,String extension,String refGenome){
		BlastREPs br = new BlastREPs(outFolder, true);
		br.createBlast(inputFile, "tblastn", 1e-30, flankingGenes,true,true);
		ArrayList<String> names = br.getNames();
		HashMap<String, Input> inputInfo = br.inputInfo;
		for (int i = 0; i < names.size(); i++) {
			System.out.println(names.get(i));
			// find all homologous genes, intervals
			Input in = inputInfo.get(names.get(i));
			HashMap<String, StringBuilder> genomeFastaLocal = ReadFasta
					.readFasta(in.genome);
			String genomeLocal = genomeFastaLocal.values().toArray(
					new StringBuilder[0])[0].toString();
			InfoTree geneTree = GetFlankingGenes.makeGeneTree(in.genbank);
			ArrayList<Info> blastIntervals = GetFlankingGenes.blastToIntervals(in.blastoutREP);
			ArrayList<Info> geneOverlaps = GetFlankingGenes.getOverlappingGenes(geneTree,
					blastIntervals,false);
			for(int j=0;j<geneOverlaps.size();j++){
				File folder=new File(outFolder+"/"+names.get(i)+"_"+extension+"/");
				folder.mkdir();
				String[] split1=geneOverlaps.get(j).info.split("\\|");
				String[] split2=split1[1].split("_");
				ArrayList<Fasta> genes=new ArrayList<Fasta>();
				ArrayList<Info> localGeneInfo=new ArrayList<Info>();
				localGeneInfo.add(new Info(geneOverlaps.get(j).getStart(),geneOverlaps.get(j).getEnd(),split1[0]));
				genes.addAll(Fasta.makeFasta(localGeneInfo, genomeLocal, false));
				ArrayList<Info> refGeneInfo=new ArrayList<Info>();
				int  startRef=Integer.parseInt(split2[split2.length-2]);
				int  endRef=Integer.parseInt(split2[split2.length-1]);
				String inf=split1[1].contains("complement")?"Reference_complement":"Reference";
				refGeneInfo.add(new Info(startRef,endRef,inf));
				genes.addAll(Fasta.makeFasta(refGeneInfo, refGenome, false));	
				Fasta.write(genes, new File(folder+"/"+split1[0]+".fas"));
			}
			System.out.println(geneOverlaps.size());
		}
	}

}
