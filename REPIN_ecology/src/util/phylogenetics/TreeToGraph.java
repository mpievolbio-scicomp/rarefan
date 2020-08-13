package util.phylogenetics;

import java.io.*;




public class TreeToGraph {
	
	
	public static void main(String[] args){
		File tree=new File(args[0]);
		File out=new File(args[1]);
		TreeToGraph ttg=new TreeToGraph(tree);
		ttg.writeGraph(out, false);
	}
	int currNodeIndex;
	Tree tree;
	StringBuffer graph=new StringBuffer();

	
	
	
	
	public TreeToGraph(File treeFile){
		currNodeIndex=0;
		makeGraph(treeFile);
		
	}
	private void makeGraph(File treeFile){
		try{
			BufferedReader br=new BufferedReader(new FileReader(treeFile));
			TreeParser tp=new TreeParser(br);
			tree=tp.tokenize(1, "whater", null);
			recursiveTranslation(tree, 0, 0,null);
			//System.out.println(info.size()+" of "+total.size());
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public void writeGraph(File graphOut,boolean append){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(graphOut,append));
			bw.write(graph.toString());
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	 void recursiveTranslation (Tree tree,int currkey, int currdepth,TreeNode parent) {
        TreeNode currNode = tree.getNodeByKey(currkey);
        int numChildren = currNode.numberChildren();
    	String name="internal"+currNodeIndex;
        if(currNode.getName().equals("")){
        	currNode.setName(name);
        }
        if(parent!=null){
        	graph.append(currNode.getName()+"\t"+parent.getName()+"\t10000\ttree\n");
        }
        
        for (int i = 0; i < numChildren; i++) {
            int childkey = currNode.getChild(i).key;
            //TreeNode childnode = tree.getNodeByKey(childkey);
            //System.out.println(name);
            currNodeIndex++;
            recursiveTranslation(tree,childkey, currdepth+1,currNode);
        }
    }


}
