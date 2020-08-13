package util;

import java.io.*;
import java.util.*;


import util.*;
import util.phylogenetics.RunTreePrograms;

//1. randomly select positions of the alignment (90%, 80%, 70%, 60%, 50%, 40%, 30%, 20%, 10%) x 100 for each, 900 files, build tree for each + consensus
//2. depending on the above outcome select the threshold for which the results are still very robust for this threshold analyze gapped positions 
//i.e. randomly select all positions with less than one gap, two gaps, ... per column
public class RefragmentAlignment {
	public static void main(String args[]){
		File alignment=new File(args[0]);
		int repeats=Integer.parseInt(args[1]);
		File masterFolder=new File(args[2]);
		File pathRAxML=new File(args[3]);
		File pathConsense=new File(args[4]);
		boolean codons=Boolean.parseBoolean(args[5]);
		int percent=Integer.parseInt(args[6]);
		File finalResults=new File(masterFolder+"/bootstrapResults.txt");
		printLine(finalResults);
		ArrayList<Fasta> fas=Fasta.readFasta(alignment);
		String name=alignment.getName().split("\\.")[0];
		System.out.println("Start...");
		File outFolder=new File(masterFolder+"/percent"+percent+"/");
		if(!outFolder.exists()){
			outFolder.mkdir();
		}
		double proportion=(percent*1.0)/100.0;
		ArrayList<File> trees=new ArrayList<File>();
		System.out.println("Choosing "+percent+"% of all nucleotide sites...");
		int seqLength=0;
		for(int j=0;j<repeats;j++){
			System.out.println("Repeat "+(j+1));
			File out=new File(outFolder+"/"+name+"_"+percent+"_"+j+".phy");
			ArrayList<Fasta> subAlignment=selectColumns(fas,proportion,codons);
			seqLength=Fasta.writePhylip(subAlignment, out,100);
			String suffix="raxml"+j;
			File RAxMLTree=new File(outFolder+"/RAxML_bestTree."+suffix);
			System.out.println("Running raxml...");
			int rand=(int)(Math.random()*1000000);
			RunTreePrograms.runRAxML(out,pathRAxML,seqLength,suffix,"",!codons,rand);
			trees.add(RAxMLTree);
		}
		File intree=new File(outFolder+"/intree");
		FileHandler.merge(trees,intree);
		RunTreePrograms.runConsense(intree, pathConsense);
		File outfile=new File(outFolder+"/outfile");
		ConsenseOutfile co=new ConsenseOutfile(outfile);
		//File outtree=new File(outFolder+"/outtree");
		//outtree.renameTo(new File(outFolder+"/outtree_"+percent));
		printBootstrapResults(seqLength,percent,0,co.getBranchPointsExcluded(),co.getBranchPointsIncluded(),finalResults);
	}
	

	public static void printLine(File out){
		try{
			if(!out.exists()){
				BufferedWriter bw=new BufferedWriter(new FileWriter(out,true));
				bw.write("The following is the start of a new Analysis:\nNucleotides\tPercent data used\tGaps\tExcluded branchpoints\tIncluded branchpoints\tProportion Included\n");
				bw.close();
			}
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static void printBootstrapResults(int nucs,int percent,int gaps,int ex,int in,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out,true));
			double proportion=(in*1.0)/(ex+in+0.0);
			bw.write(nucs+"\t"+percent+"\t"+gaps+"\t"+ex+"\t"+in+"\t"+proportion+"\n");
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static ArrayList<Integer> initPositions(int size){
		ArrayList<Integer> pos=new ArrayList<Integer>();
		for(int i=0;i<size;i++){
			pos.add(i);
		}
		return pos;
	}
	
	public static ArrayList<Fasta> selectColumns(ArrayList<Fasta> fas,double proportion,boolean codons){
		int size=fas.get(0).getSequence().length();
		int selectCols=(int)(size*proportion);
		ArrayList<Integer> positions=initPositions(size);
		ArrayList<StringBuffer> sequences=initSequences(fas.size());
		ArrayList<String> initSeqs=toStringArray(fas);
		ArrayList<Integer> randomPositions=new ArrayList<Integer>();
		int adjust=codons?3:1;
		for(int i=0;i<selectCols;i+=adjust){
			
			int rand=(int)(Math.random()*(size-(i)));
			int randomPosition=positions.get(rand);
			int frame=codons?randomPosition%3:0;
			randomPositions.add(randomPosition-frame+0);
			positions.remove(rand-frame);
			if(codons){
				randomPositions.add(randomPosition-frame+1);
				randomPositions.add(randomPosition-frame+2);
				positions.remove(rand-frame);
				positions.remove(rand-frame);
			}

			//if(i%1000==0)System.out.println(i);
		}
		sequences=addPositions(fas.size(),randomPositions,initSeqs);
		return makeFasta(sequences,fas);
	}
	
	public static ArrayList<String> toStringArray(ArrayList<Fasta> fas){
		ArrayList<String> seqs=new ArrayList<String>();
		for(int i=0;i<fas.size();i++){
			seqs.add(fas.get(i).getSequence());
			
		}
		return seqs;
	}
	
	public static ArrayList<Fasta> makeFasta(ArrayList<StringBuffer> seqs,ArrayList<Fasta> fas){
		ArrayList<Fasta> newFas=new ArrayList<Fasta>();
		for(int i=0;i<fas.size();i++){
			newFas.add(new Fasta(fas.get(i).getIdent(),seqs.get(i).toString()));
		}
		return newFas;
	}
	
	public static ArrayList<StringBuffer> initSequences(int length){
		ArrayList<StringBuffer> seqs=new ArrayList<StringBuffer>();
		for(int i=0;i<length;i++){
			seqs.add(new StringBuffer());
		}
		return seqs;
	}
	
	public static ArrayList<StringBuffer> addPositions(int size,ArrayList<Integer> pos,ArrayList<String> fas){
		ArrayList<StringBuffer> seqs=new ArrayList<StringBuffer>();
		for(int i=0;i<size;i++){
			int posSize=pos.size();
			StringBuffer temp=new StringBuffer();
			String sequence=fas.get(i);
			for(int j=0;j<posSize;j++){
				temp.append(sequence.charAt(pos.get(j)));
			}
			seqs.add(temp);
		}
		return seqs;
	}
}
