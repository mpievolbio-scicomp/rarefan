package util.phylogenetics;

import java.io.*;
import java.util.*;

import util.*;
import util.phylogenetics.RunTreePrograms;

public class DetermineMasterTrees {
	public static void main(String args[]){
		File phylipPath=new File(args[0]);
		File masterFolder=new File(args[1]);
		String treePuzzlePath=args[2];
//		File neighborPath=new File(phylipPath+"/neighbor");
//		File dnaDistPath=new File(phylipPath+"/dnadist");
		File treedistPath=new File(phylipPath+"/treedist");
		File treedraw=new File(phylipPath+"/drawtree");
		File fontfile=new File(phylipPath+"/font1");

		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(masterFolder+"/treedistances.txt"));
			bw.write("TreePosition\trange\ttreeDif\n");
			for(int i=3;i<args.length;i++){
				ArrayList<Fasta> allFasta=new ArrayList<Fasta>();
				File outFolder=new File(masterFolder+"/"+args[i]);
				if(outFolder.isDirectory()){
					File list2[]=outFolder.listFiles();
					for(int j=0;j<list2.length;j++){
						if(list2[j].getName().endsWith(".fas")){
							allFasta.addAll(Fasta.readFasta(list2[j]));
						}
					}
				}
				File out=new File(outFolder+"/alignment.phy");
				File masterTree=new File(out+".tree");
				Fasta.writePhylip(allFasta, out, 10);
				RunTreePrograms.runTreePuzzle(out, treePuzzlePath);
//				RunTreePrograms.runProgram(dnaDistPath.toString(),out.toString()+"\nT\n1.0\nY\n" , outFolder);
//				File outfile=new File(outFolder+"/outfile");
//				File dist=new File(outFolder+"/master.dist");
//				outfile.renameTo(dist);
//				RunTreePrograms.deleteConsenseFiles(outFolder);
//				RunTreePrograms.runProgram(neighborPath.toString(),dist.toString()+"\nY\n" , outFolder);
//				outfile=new File(outFolder+"/outtree");
//				File masterTree=new File(outFolder+"/master.tree");
//				outfile.renameTo(masterTree);
				String suffix=getSuffix(outFolder);
				File calcTree=new File(masterFolder+"/S1111"+suffix+"/PolySeqOut_NoGenes/polymorphisms_move.phy.tree");
				File merge=new File(outFolder+"/merge.tree");
				FileHandler.merge(calcTree,masterTree, merge);
				RunTreePrograms.deleteConsenseFiles(outFolder);
				RunTreePrograms.runProgram(treedistPath.toString(), merge.toString()+"\nY\n", outFolder);
				File outfile=new File(outFolder+"/outfile");
				parseTreeDistance(bw,outfile,"");
				RunTreePrograms.deleteConsenseFiles(outFolder);
				RunTreePrograms.runProgram(treedraw.toString(), masterTree.toString()+"\n"+fontfile+"\nV\nN\nY\n", outFolder);
				File plotfile=new File(outFolder+"/plotfile");
				plotfile.renameTo(new File(outFolder+"/mastertree.ps"));
				RunTreePrograms.deleteConsenseFiles(calcTree.getParentFile());
				RunTreePrograms.runProgram(treedraw.toString(), calcTree.toString()+"\n"+fontfile+"\nV\nN\nY\n", calcTree.getParentFile());
				plotfile=new File(calcTree.getParent()+"/plotfile");
				plotfile.renameTo(new File(calcTree.getParent()+"/mastertree.ps"));
				
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}

	}
	
	public static void parseTreeDistance(BufferedWriter bw,File in,String dataset)throws IOException{
		BufferedReader br=new BufferedReader(new FileReader(in));
		String line="";
		while((line =br.readLine())!=null){
			if(line.startsWith("Trees 1 and 2:")){
				String split[]=line.split("\\s+");
				String folder[]=in.getParent().split("_");
				bw.write(dataset+"\t"+folder[1]+"\t"+folder[2]+"\t"+split[4]+"\n");
			}
			
		}
		br.close();
	}
	
	
	public static int parseTreeDistance(File in)throws IOException{
		BufferedReader br=new BufferedReader(new FileReader(in));
		String line="";
		while((line =br.readLine())!=null){
			if(line.startsWith("Trees 1 and 2:")){
				br.close();
				String[] split=line.split("\\s+");
				return Integer.parseInt(split[4]);
			}
			
		}
		br.close();
		return -1;
		
	}
	
	public static String getSuffix(File in){
		String name=in.getName();
		String split[]=name.split("_");
		StringBuffer suff=new StringBuffer();
		for(int i=1;i<split.length;i++){
			suff.append("_"+split[i]);
		}
		return suff.toString();
	}
}
