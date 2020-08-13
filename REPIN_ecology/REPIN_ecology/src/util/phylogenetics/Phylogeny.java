package util.phylogenetics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import util.FileHandler;

public class Phylogeny {

	public static File makeIncorrectTree2(File path){
		File out=new File(path+"/inc2.tree");
		if(out.exists()){
			return out;
		}else{
			try{
				BufferedWriter bw=new BufferedWriter(new FileWriter(out));
				bw.write("((S21,S11),S12,S22);");
				bw.close();
			}catch(IOException e){
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		return out; 
	}
	

//	
//	public static HashMap<String,HashMap<String,Double>> getDistanceMatrix(File tree,File Rscript){
//		File distMatrix=new File(tree.getParentFile()+"/distMat.txt");
//		R_functions.makeDistanceMatrix(tree,distMatrix,Rscript);
//		HashMap<String, HashMap<String,Double>> matrix=convertToHash(distMatrix);
//		
//		return matrix;
//	}
//	


	public static HashMap<String,HashMap<String,Double>> convertToHash(File matrix){
		HashMap<String,HashMap<String,Double>> hm=new HashMap<String, HashMap<String,Double>>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(matrix));
			String line;
			int i=0;
			ArrayList<String> leaves=new ArrayList<String>();
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				if(i==0){
					for(int j=0;j<split.length;j++){
						String temp=split[j].replace("\"", "");
						leaves.add(temp);
					}
				}else{
					HashMap<String,Double> temp=new HashMap<String, Double>();
					String leave=split[0].replace("\"", "");
					for(int j=1;j<split.length;j++){
						temp.put(leaves.get(j-1), Double.parseDouble(split[j]));
					}
					hm.put(leave,temp);
				}
				i++;
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return hm;
	}
	
	public static Tree readTree(File newick){
		try{
			BufferedReader br=new BufferedReader(new FileReader(newick));
			TreeParser tp=new TreeParser(br);
			Tree tree=tp.tokenize(1, "whatever", null);
			br.close();
			return tree;
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	public static File makeIncorrectTree1(File path){
		File out=new File(path+"/inc1.tree");
		if(out.exists()){
			return out;
		}else{
			try{
				BufferedWriter bw=new BufferedWriter(new FileWriter(out));
				bw.write("((S22,S11),S12,S21);");
				bw.close();
			}catch(IOException e){
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		return out;
	}
	public static File makeRealTree(File path){
		File out=new File(path+"/realTree.tree");
		if(out.exists()){
			return out;
		}else{
			try{
				BufferedWriter bw=new BufferedWriter(new FileWriter(out));
				bw.write("((S22,S21),S12,S11);");
				bw.close();
			}catch(IOException e){
				e.printStackTrace();
				System.exit(-1);
			}
		}
		return out;
	}
	public static int compare(File tree1,File tree2,File treedistPath)throws IOException{
		
		File outFolder=tree1.getParentFile();
		File merge=new File(outFolder+"/merge.tree");
		FileHandler.merge(tree1,tree2, merge);
		RunTreePrograms.deleteConsenseFiles(outFolder);
		RunTreePrograms.runProgram(treedistPath.toString(), merge.toString()+"\nD\nY\n", outFolder);
		File outfile=new File(outFolder+"/outfile");
		int diff=DetermineMasterTrees.parseTreeDistance(outfile);
		RunTreePrograms.deleteConsenseFiles(outFolder);
		return diff;
	}
}
