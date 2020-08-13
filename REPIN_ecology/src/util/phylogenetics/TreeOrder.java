package util.phylogenetics;

import java.io.*;
import java.util.ArrayList;


public class TreeOrder {
	
	public static void main(String args[]){
		File tree=new File(args[0]);
		File out=new File(args[1]);
		TreeOrder tdel=new TreeOrder(tree);
		tdel.writeTree(out);
	}
	
	Tree tree;
	ArrayList<String> order;
	public TreeOrder(File treeFile){
		calcTreeOrder(treeFile);
	}
	
	public void calcTreeOrder(File treeFile){
		try{
			BufferedReader br=new BufferedReader(new FileReader(treeFile));
			TreeParser tp=new TreeParser(br);
			tree=tp.tokenize(1, "whater", null);
			order=new ArrayList<String>();
			recursive(tree, 0, 0);
			//System.out.println(info.size()+" of "+total.size());
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	

	
	public ArrayList<String> getOrder(){
		return order;
	}
	
	public void writeTree(File out){
		TreeWriter.writeTree(tree, out.toString(),false);
	}
	 void recursive (Tree tree,int currkey, int currdepth) {
        TreeNode currNode = tree.getNodeByKey(currkey);
        int numChildren = currNode.numberChildren();
        if(numChildren==0){
        	order.add(currNode.getName());
        }
        for (int i = 0; i < numChildren; i++) {
            int childkey = currNode.getChild(i).key;
            recursive(tree,childkey, currdepth+1);
        }
    }

}
