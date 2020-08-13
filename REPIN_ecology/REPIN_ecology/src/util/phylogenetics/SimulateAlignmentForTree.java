package util.phylogenetics;

import java.io.*;
import java.util.*;

import util.DNAmanipulations;
import util.Fasta;


public class SimulateAlignmentForTree {
	
	public static void main(String args[]){
		File tree=new File(args[0]);
		int length=Integer.parseInt(args[1]);
		double GC=Double.parseDouble(args[2]);
		File out=new File(tree.getParentFile()+"/simulatedAlignment.fas");
		makeTree(tree,out,length,GC);
	}
	
	public static void makeTree(File treeFile,File out,int length,double GC){
			ArrayList<Fasta> alignment=getAlignment(treeFile, length, GC);
			getAlignment(treeFile,length,GC);
			Fasta.write(alignment, out);

	}

	public static ArrayList<Fasta> getAlignment(File treeFile,int length,double GC){
		ArrayList<Fasta> alignment=new ArrayList<Fasta>();
		try{
			String root=DNAmanipulations.generateRandomSequence(length, GC);
			BufferedReader br=new BufferedReader(new FileReader(treeFile));
			TreeParser tp=new TreeParser(br);
			Tree tree=tp.tokenize(1, "whatever", null);
			recursive_simulateAlignment(tree, 0, 0,alignment,root,GC);
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return alignment;
	}
	
	static void recursive_simulateAlignment (Tree tree,int currkey, int currdepth,ArrayList<Fasta> alignment,String sequence,double GC) {
        TreeNode currNode = tree.getNodeByKey(currkey);
        
        int numChildren = currNode.numberChildren();
        for (int i = 0; i < numChildren; i++) {
            int childkey = currNode.getChild(i).key;
            TreeNode childnode = tree.getNodeByKey(childkey);
            String name=childnode.getName();
            //System.out.println(name);
            Double weight=(double)childnode.getWeight();
            String newSeq=DNAmanipulations.randomizeSequence(sequence, weight, GC);
            if(!name.equals("")){
                alignment.add(new Fasta(name,newSeq));
            }
            recursive_simulateAlignment(tree,childkey, currdepth+1,alignment,newSeq,GC);
        }
    }


}
