package util.phylogenetics;

import java.io.*;
import java.util.*;
import java.util.AbstractMap.SimpleEntry;

import util.Fasta;
/**
 * Does not currently work...since it does not seem to work for our problem
 * @author bertels
 *
 */
public class FitnessProxyScore {
	public class Compound{
		double weight;
		double time;
		public Compound(double weight,double time){
			this.weight=weight;
			this.time=time;
		}
	}
	
	public static void main(String args[]){
		File alignment=new File(args[0]);
		FitnessProxyScore fps=new FitnessProxyScore(alignment);
		fps.writeFPS(new File(alignment.getParent()+"/"+alignment.getName()+".fps"));
	}
	
	Tree tree;
	ArrayList<Fasta>  fasta;
	HashMap<Integer,Compound> treedata=new HashMap<Integer, Compound>(); 
	Alignment alignment;
	String PhymlPath="/Users/bertels/Programs/PhyML-3.1/PhyML-3.1_macOS-MountainLion";
	//defined as distance to leave that is most divergent
	HashMap<Integer,Double> distanceToLeaves=new HashMap<Integer, Double>();
	HashMap<Integer,Double> distanceToRoot=new HashMap<Integer, Double>();
	HashMap<Integer,HashSet<Integer>> ancestorMap=new HashMap<Integer, HashSet<Integer>>();
	HashMap<Integer,HashMap<Integer,Integer>> mrca=new HashMap<Integer, HashMap<Integer,Integer>>();
	ArrayList<SimpleEntry<String,Double>> fpsScores;
	String root;
	double T2=0;
	double tstar;
	public FitnessProxyScore(File alignment){
		determineTopology(alignment);
		
		//UPGMA distMat=new UPGMA(this.alignment);
		//System.out.println("final tree: "+setUPGMADistances(distMat));
		setMaxLeafDistance(tree.getRoot().getKey());
		setMaxRootDistance(tree.getRoot().getKey(), 0.0);
		setAncestorList(tree.getRoot().getKey(),new HashSet<Integer>());
		System.out.println("root: "+tree.getRoot().name);
		initData(tree.getRoot().getKey()); 
		T2=calculateT2root();
		tstar=T2*0.5;
		calculateAllFPS();
		
	}
	
	public void writeFPS(File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<fpsScores.size();i++){
				bw.write(fpsScores.get(i).getKey()+"\t"+fpsScores.get(i).getValue()+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private void setAncestorList(int currkey,HashSet<Integer> ancestorList){
		TreeNode currNode = tree.getNodeByKey(currkey);
		HashSet<Integer> temp=new HashSet<Integer>();
		temp.addAll(ancestorList);
		this.ancestorMap.put(currkey, temp);
		ancestorList.add(currkey);
		int numChildren = currNode.numberChildren();
		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;
			setAncestorList(childkey,ancestorList);


		}
		ancestorList.remove(currkey);
	}
	
	private TreeNode getMRCA(int k1,int k2){
		int rootkey=this.tree.getRoot().getKey();
		TreeNode tn=tree.getNodeByKey(k1);
		int currkey=k1;
		HashSet<Integer> k2set=ancestorMap.get(k2);
		while(currkey!=rootkey){
			if(k2set.contains(currkey)){
				return tn;
			}
			tn=tn.parent;
			currkey=tn.getKey();
		}
		System.err.println("Could not find common ancestor! "+tree.getNodeByKey(k1).getName()+" "+tree.getNodeByKey(k2).getName());
		return null;
	}
	
	public double calculateT2root(){
		LinkedList<TreeNode> leaves=tree.getLeaves(tree.getRoot());
		double sum=0;
		int count=0;
		for(int i=0;i<leaves.size();i++){
			int keyi=leaves.get(i).getKey();
			mrca.put(keyi, new HashMap<Integer,Integer>());
			for(int j=i+1;j<leaves.size();j++){
				int keyj=leaves.get(j).getKey();
				int currMRCA=getMRCA(keyi,keyj).getKey();
				mrca.get(keyi).put(keyj, currMRCA);
				count++;
				sum+=distanceToRoot.get(currMRCA);
				
			}
		}
		return sum/count;
	}
	public double calculateT2leaves(){//TODO!
		return 0.0;
	}
	/**Step 1
	 * Determine the topology of the tree using PHYML.
	 */
	public void determineTopology(File alignment){
		fasta=Fasta.readFasta(alignment);
		this.alignment=new Alignment(fasta);
		this.alignment.deleteDuplicates();
		 root="*"+this.alignment.idents.get(0);
		this.alignment.idents.set(0,		root );
		//this.alignment.renameSequences(1);
		File folder=alignment.getParentFile();
		String name=alignment.getName().split("\\.")[0];
		File phy=new File(folder+"/"+name+".phy");
		fasta=this.alignment.toFasta();
		Fasta.writePhylip(fasta, phy, 100);
		File tree=RunTreePrograms.runPhyml(phy, PhymlPath," -o n");
		this.tree=Phylogeny.readTree(tree);
		setRoot();
		//parseTree(this.tree, this.tree.getRoot().getKey());
	}
	
	private void calculateAllFPS(){
		LinkedList<TreeNode> leaves=tree.getLeaves(tree.getRoot());
		fpsScores=new ArrayList<AbstractMap.SimpleEntry<String,Double>>();
		for(int i=0;i<leaves.size();i++){
			fpsScores.add(new SimpleEntry<String, Double>(leaves.get(i).getName(),calculateFPSLeaves(leaves.get(i))));
		}
	}
	
	public double calculateFPSLeaves(TreeNode tn){
		int currkey=tn.parent.key;
		int rootkey=tree.getRoot().getKey();

		double sum=0;
		int prevKey=tn.getKey();

		while(currkey!=rootkey){
			tn=tn.parent;

			currkey=tn.getKey();
			sum+=phi(distanceToLeaves.get(currkey)/T2)*(treedata.get(currkey).weight-treedata.get(prevKey).weight);
			prevKey=currkey;
		}
		return sum;
		
	}
	
	public double calculateFPSRoot(TreeNode tn){
		int currkey=tn.parent.key;
		int rootkey=tree.getRoot().getKey();

		double sum=0;
		int prevKey=tn.getKey();

		while(currkey!=rootkey){
			tn=tn.parent;

			currkey=tn.getKey();
			sum+=phi(distanceToRoot.get(currkey)/T2)*(treedata.get(currkey).weight-treedata.get(prevKey).weight);
			prevKey=currkey;
		}
		return 1/sum;
		
	}
	
	private double phi(double t){
		System.out.println(tstar+" "+t+" "+T2);
		 return 1/(1+Math.exp(5*(t/tstar-1)));
	}
	
	/**Step 3
	 * Re-weigh the tree based on the given alignment so that the distance from the root to the leaves is equal.
	 * Start with the leaves and then iteratively climb up the tree and assign lengths.
	 */
	public void calculateUPGMATree(){
		
	}
	
	
	
	/**
	 * Step 2: Root is determined in tree building step
	 * Divide strains into two groups via the parsimony method by minimizing the number of polymorphisms that exist in both groups.
	 * */
	
	
	private void printIdents(Alignment alg){
		for(int i=0;i<alg.idents.size();i++){
			System.out.println(alg.idents.get(i));
		}
	}
	
	private void printLeaves(LinkedList<TreeNode> leaves){
		for(int i=0;i<leaves.size();i++){
			System.out.println(leaves.get(i).getName());
		}
	}
	
	private void parseTree(Tree t,int currkey){
		TreeNode currNode = tree.getNodeByKey(currkey);
		int numChildren = currNode.numberChildren();
		System.out.println(currNode.getName());
		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;
			parseTree(t,childkey);
		}
	}
	

	private void setRoot(){
		LinkedList<TreeNode> leaves=tree.getLeaves(tree.getRoot());
		for(int i=0;i<leaves.size();i++){
			if(leaves.get(i).getName().equals(root)){
				TreeNode root=tree.getRoot();
				root.parent=leaves.get(i);
				TreeNode lparent=leaves.get(i).parent;
				leaves.get(i).addChild(root);
				lparent.deleteChild(leaves.get(i));
				leaves.get(i).parent=null;
				tree.setRootNode(leaves.get(i));
				tree.updateTree();
				break;
			}
		}
		setLeafNumbers(tree.getRoot().getKey());
		setNames(tree.getRoot().getKey());
	}

/*	Does not seem to work for our example, I will adjust the method to our particular problem...
 * private String setUPGMADistances(UPGMA distMat){
		LinkedList<TreeNode> leaves=tree.getLeaves(tree.getRoot());
		int count=0;
		for(int i=0;i<leaves.size();i++){
			if(leaves.get(i).getName().equals(root)){
				TreeNode newNode=new TreeNode( );
				TreeNode oldParent=leaves.get(i).parent;

				oldParent.deleteChild(leaves.get(i));
				newNode.addChild(oldParent);
				oldParent.parent=newNode;
				newNode.addChild(leaves.get(i));
				leaves.get(i).parent=newNode;
				tree.setRootNode(newNode);
				tree.updateTree();
				break;
			}
		}
		 leaves=tree.getLeaves(tree.getRoot());
		 System.out.println("rightmost: "+leaves.get(leaves.size()-1).getName());
		while(leaves.size()>1){
			LinkedList<TreeNode> newLeaves=new LinkedList<TreeNode>();
			HashMap<TreeNode,TreeNode> parents=new HashMap<TreeNode,TreeNode>();
			System.out.println(count++);
			for(int i=0;i<leaves.size();i++){
				System.out.println("name: "+leaves.get(i).getName());
				TreeNode parent=leaves.get(i).parent();
				if(parents.containsKey(parent)){
					System.out.println("#children: "+parent.numberChildren());
					TreeNode c1=parent.getChild(0);
					TreeNode c2=parent.getChild(1);
					String key1=c1.getName();
					String key2=c2.getName();
					System.out.println(key1+" "+key2);
					UPGMA.MergeResult mr=distMat.mergeEntries(key1, key2);
					parent.setName(mr.k1k2);
					c1.setWeight(mr.distK1);
					c2.setWeight(mr.distK2);
					newLeaves.add(parent);
					parents.remove(parent);
				}else{
					parents.put(parent, leaves.get(i));
				}
			}
			TreeNode[] remainingLeaves=parents.keySet().toArray(new TreeNode[0]);
			for(int i=0;i<remainingLeaves.length;i++){
				newLeaves.add(parents.get(remainingLeaves[i]));
			}
			leaves=newLeaves;
		}
		return leaves.get(0).getName();
		
	}*/

	
	private void setMaxRootDistance(int currkey,double distance){
		TreeNode currNode = tree.getNodeByKey(currkey);
		distance+=currNode.weight;
		distanceToRoot.put(currNode.getKey(), distance);
		int numChildren = currNode.numberChildren();
		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;
			setMaxRootDistance(childkey,distance);
			
			
		}
	}
	

	
	private void setMaxLeafDistance(int currkey){
		TreeNode currNode = tree.getNodeByKey(currkey);

		int numChildren = currNode.numberChildren();
		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;
			setMaxLeafDistance(childkey);
			
			
		}
		setLeafDistance(currNode);
	}
	
	private void setLeafDistance(TreeNode tn){
		if(tn.isLeaf()){
			distanceToLeaves.put(tn.getKey(), 0.0);
		}else{
			int numChildren=tn.numberChildren();
			double max=Double.NEGATIVE_INFINITY;
			for(int i=0;i<numChildren;i++){
				int key=tn.getChild(i).getKey();
				double currDist=distanceToLeaves.get(key)+tn.getChild(i).getWeight();
				if(max<currDist){
					max=currDist;
				}
				
			}
			distanceToLeaves.put(tn.getKey(),max);
		}
	}
	
	private void setLeafNumbers(int currkey){
		TreeNode currNode = tree.getNodeByKey(currkey);

		int numChildren = currNode.numberChildren();
		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;
			setLeafNumbers(childkey);
			
			
		}
		currNode.setNumberLeaves();
	}
	
	private void setNames(int currkey){
		TreeNode currNode = tree.getNodeByKey(currkey);

		int numChildren = currNode.numberChildren();
		
		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;
			setNames(childkey);
			
			
		}
		if(!currNode.isLeaf()&&!currNode.equals(tree.getRoot())){
			StringBuffer name=new StringBuffer();
			for(int i=0;i<numChildren;i++){
				name.append(currNode.getChild(i).getName());
			}
			currNode.setName(name.toString());
		}
	}
	
	private void initData(int currkey){
		TreeNode currNode = tree.getNodeByKey(currkey);
		double weight=currNode.numberLeaves;
		double time=distanceToLeaves.get(currkey);
		Compound c=new Compound(weight,time);
		System.out.println(currNode.getName()+" "+weight+" "+time);
		treedata.put(currkey, c);
		currNode.getWeight();
		int numChildren = currNode.numberChildren();
		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;
			initData(childkey);
		}
	}
}
