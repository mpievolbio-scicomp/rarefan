package util.phylogenetics;

import java.io.*;


public class TreeDeleteEdgeLength {
	
	public static void main(String args[]){
		File tree=new File(args[0]);
		File out=new File(args[1]);
		TreeDeleteEdgeLength tdel=new TreeDeleteEdgeLength(tree);
		tdel.writeTree(out);
	}
	
	Tree tree;
	
	public TreeDeleteEdgeLength(File treeFile){
		annotateTree(treeFile);
	}
	
	public void annotateTree(File treeFile){
		try{
			BufferedReader br=new BufferedReader(new FileReader(treeFile));
			TreeParser tp=new TreeParser(br);
			tree=tp.tokenize(1, "whater", null);
			recursive_ChangeEdges(tree, 0, 0,0.05);
			//System.out.println(info.size()+" of "+total.size());
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	

	

	
	public void writeTree(File out){
		TreeWriter.writeTree(tree, out.toString(),false);
	}
	 void recursive_ChangeEdges (Tree tree,int currkey, int currdepth,double branchlength) {
        TreeNode currNode = tree.getNodeByKey(currkey);
        
        int numChildren = currNode.numberChildren();
        if(numChildren==0){
        	currNode.setWeight(branchlength);
        }else{
        	branchlength=branchlength/2;
        	currNode.setWeight(branchlength);
        }
        for (int i = 0; i < numChildren; i++) {
            int childkey = currNode.getChild(i).key;
            recursive_ChangeEdges(tree,childkey, currdepth+1,branchlength);
        }
    }

}
