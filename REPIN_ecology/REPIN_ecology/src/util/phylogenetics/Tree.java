package util.phylogenetics;

import java.text.*;
import java.util.*;


public class Tree {

	/** The list of nodes of the tree indexed by their keys, indexed by key */ 
	public ArrayList<TreeNode> nodes; 


	/** 
	 * Most internal nodes don't have names. Do we assign a unique
	 * name to each of them? No! each node has a key and the key is unique
	 * for nodes. 
	 */
	private HashMap<String,TreeNode> nodesByName; 

	/** key should be unique for each tree, set by object that creates trees  */
	private int key;

	/** Leaf counter, for determining grid size, making arrays for tree comparisons */
	private int numLeaves = 0;

	// reference for array of leaves in SC.cullingObject
	/** Split axis reference for leaf recovery (leaves are attached to split line culling objects) */
	//private StaticSplitAxis leafSplitAxis;

	/**
	 * Default tree constructor.  Nodes are created by parser and added in later.
	 *
	 */
	public Tree() {
		root = new TreeNode();
		nodes = new ArrayList<TreeNode>();
		nodesByName = new HashMap<String,TreeNode>();
	}

	/**
	 * Copy constructor used to create versions of trees that are identical to the supplied input tree.
	 * The copying should also make copies of the nodes in the given tree with the node copy constructor.
	 * (This function isn't working properly)
	 * @param treeToCopy Tree used to make a copy.
	 */
	public Tree(Tree treeToCopy)
	{
		// TODO: make this work with copy constructors (this constructor is only used in matrix mode)
		fileName = treeToCopy.fileName;
		height = treeToCopy.height;
		key = treeToCopy.key;
//		leafSplitAxis = new SplitAxis(treeToCopy.leafSplitAxis); // not implemented
		nexusIndex = treeToCopy.nexusIndex;
		numLeaves = treeToCopy.numLeaves;
		
		copyTree(treeToCopy);
		//nodes = copyNodes(treeToCopy.nodes);
		//nodesByName = copyNodesByName(treeToCopy.nodesByName);
		//root = TreetreeToCopy.root;
	}
	public Tree rerootTree(int key){
		Tree newTree=new Tree(this);
		
		TreeNode tn=newTree.getNodeByKey(key);
		if(!tn.isRoot()){
			if(tn.parent.parent!=null){
				tn.parent.parent.deleteChild(tn.parent);
				tn.parent.parent.setWeight(tn.parent.parent.weight+tn.parent.weight);
				tn.parent.addChild(tn.parent.parent);
				tn.parent.weight=0;
			}
			newTree.setRootNode(tn.parent);

			tn.parent.deleteChild(tn);
			updateTree();
		}
		

		return newTree;
	}
	
	public Tree getSubTree(int key){
		Tree newTree=new Tree(this);
		TreeNode tn=newTree.getNodeByKey(key);
		tn.parent=null;
		newTree.setRootNode(tn);
		
		


		

		return newTree;
	}
	
	private void copyTree(Tree t){
		root=new TreeNode(t.getRoot());
		nodes=new ArrayList<TreeNode>();
		nodesByName=new HashMap<String, TreeNode>();
		postProcess();
		setUpNameLists();
		
		
	}

	/**
	 * Clean up method, called when the tree is deleted.
	 * @see TreeNode#close()
	 *
	 */   
	public void close(){
		TreeNode pren = root.leftmostLeaf;		
		for(TreeNode n = pren.preorderNext; n!=null; n=n.preorderNext) {
			n.close();				 
		}
	}

	/**
	 * Calls to #close() when tree is deleted. 
	 */
	protected void finalize() throws Throwable {

		try {
			close();
		}
		finally {
			super.finalize();     
		}
	}

	/**
	 * Returns the number of interior nodes in this tree.  For debugging.
	 * @return Total number of nodes minus the number of leaves.
	 */
	//private int getInteriorCount() { return nodes.size() - numLeaves;}
	/**
	 * Returns the node count, for internal and leaf nodes.
	 * @return Size of the {@link #nodes} array, which contains all nodes.
	 */
	protected int getTotalNodeCount() { return nodes.size();}

	/**
	 * Returns the node indexed by the given key.
	 * @param key Key of the node to retrieve.
	 * @return Treenode referenced by the given key.
	 */
	public TreeNode getNodeByKey(int key){ if (key >= nodes.size()) return null; return (TreeNode) nodes.get(key);}
	/**
	 * Returns the node given by the string.
	 * @param s Name/label of node to retrieve.
	 * @return Treenode referenced by the given name.
	 */
	public TreeNode getNodeByName(String s){ 
		return (TreeNode) nodesByName.get(s);
	}

	/**
	 * Height of tree, which is also the longest path from the root to some leaf node.
	 */
	private int height = 0;
	/**
	 * Accessor for height of tree.  This is also the longest path from the root to some leaf node.
	 * @return value of {@link #height}.
	 */
	public int getHeight() { return height; }

	/** Mutator for key
	 * @param i New value for {@link #key}.
	 */
	public void setKey(int i) {key = i;}
	/** Accessor for key.
	 * @return Value of {@link #key}.
	 */
	public int getKey() {return key;}
	/**
	 * File name accessor.
	 * @return value of {@link #fileName}.
	 */
	public String getName() {return fileName;}
	/** Left most leaf accessor.  This is the "min leaf"
	 * @return root's left most leaf, which is the smallest indexed leaf node in the tree.
	 */
	public TreeNode getLeftmostLeaf() { return root.leftmostLeaf; }
	/** Root accessor.
	 * @return Value of {@link #root}*/
	public TreeNode getRoot() { return root;}
	public void setRootNode(TreeNode newRoot) { this.root = newRoot; }

	/**
	 * File name for this tree.
	 */
	private String fileName = null; // the filename
	/**
	 * Index of tree in nexus file, if this tree is from a nexus file.
	 */
	private int nexusIndex = 0; // the number (>0 for nexus, appended to nexus filename)
	/**
	 * Root node of this tree
	 */
	protected TreeNode root=null;

	/**
	 * Sets the file name.  Copies the value for some reason.
	 * @param tn New value for file name.
	 */
	public void setFileName(String tn) {
		fileName = new String(tn);
	}

	/**
	 * Returns the number of leaves in this tree.
	 * @return value of {@link #numLeaves}.
	 */
	public int getLeafCount() {
		return numLeaves;
	}

	/**
	 * Post processing includes computing size of each node, 
	 * linking nodes in different order, etc.
	 * Sets left and right-most leaves of the tree.
	 * Computes and stores pre- and post-orders for leaves.
	 * Can't do minmax until after linkNodesInPreorder is called 
	 * to set index values!
	 *
	 * @see     TreeNode
	 */
	public void postProcess() {
		preorderPostProcess();
		linkLeaves();
//		System.out.println("progress bar updated: min:" + jpb.getMinimum() + " max:" + jpb.getMaximum() + " value:" + jpb.getValue());
	}

	/**
	 * 
	 * Traverses the tree in pre-order, stores the ordering in the preorderNext field of TreeNodes
	 * Sets node count for the tree.
	 *
	 * @see     TreeNode
	 */
	private void preorderPostProcess()
	{
		// munge names here, names become fully qualified, labels are what names were
		//final char separator = '/'; // separator between name fields
		// arbitrary seen by users in search, no parsing on this is required later
		int index = 0;
		height = 1;
		for(TreeNode n = root; n != null; n = n.preorderNext)
		{
			n.label = n.name;
			n.key = index++;
			nodes.add(n);
			if(n.name != null && n.name.length() > 0) {
				// don't put an empty string in the
				// hash table
				nodesByName.put(n.name, n);
			}
			n.height = (null != n.parent) ? n.parent.height+1 : 1;
			height = (n.height > height) ? n.height : height;
		}

	}

	/**
	 * Traverse the tree and initialize the {@link #nodesByName} and {@link #nodes} data structures.
	 * Used when modifying the names of nodes as well as initialization.
	 *
	 */
	public void setUpNameLists()
	{
		nodes = new ArrayList<TreeNode>();
		nodesByName = new HashMap<String,TreeNode>();
		//final char separator = '/'; // separator between name fields
		for(TreeNode n = root; n != null; n = n.preorderNext)
		{
			System.out.println(n.name);
			n.label = n.name;
			nodes.add(n);
			if(n.name != null && n.name.length() > 0) {
				// don't put an empty string in the
				// hash table
				nodesByName.put(n.name, n);
			}
			n.height = (null != n.parent) ? n.parent.height+1 : 1;
			height = (n.height > height) ? n.height : height;
		}
	}
	
	/**
	 * Wrapper for initiating {@link #linkSubtreeNodesInPreorder(TreeNode)} with the root node.
	 */
//	private void linkNodesInPreorder() {
//
//		linkSubtreeNodesInPreorder(root);
//
//	}

	/**
	 * Traverses the subtree rooted at TreeNode n in pre-order, stores the
	 * ordering in the preorderNext field of TreeNodes. 
	 * @param   n the root of the subtree
	 *
	 * @see     TreeNode
	 */
//	private void linkSubtreeNodesInPreorder(TreeNode n) {
//
//		if(n.isLeaf()) return;
//		for(int i=0; i<n.numberChildren(); i++) {
//			linkSubtreeNodesInPreorder(n.getChild(i));
//		}
//
//		n.preorderNext = n.firstChild();
//		for(int i=0; i<n.numberChildren()-1; i++) {
//			n.getChild(i).rightmostLeaf.preorderNext = n.getChild(i+1);
//		}
//		n.rightmostLeaf.preorderNext = null;
//	}

	public void updateTree(){
		linkNodes(root,new ArrayList<Integer>());
		root.rightmostLeaf.preorderNext=null;
		setUpNameLists();

		preorderPostProcess();
	}
	
	private void linkNodes(TreeNode currNode,ArrayList<Integer> test){

		int numChildren = currNode.numberChildren();
		currNode.key=test.size();
		test.add(0);
		for (int i = 0; i < numChildren; i++) {
			linkNodes(currNode.getChild(i),test);
			
			
		}
		currNode.setExtremeLeaves();
		currNode.linkNodesInPostorder();
		currNode.linkNodesInPreorder();
	}
	
	/**
	 * 
	 * Links leaves of the tree in pre-order,
	 * check to see whether leaves have distinct names.
	 * If leaves have the same name, add a suffix index separated by " "
	 *
	 * @see     #linkNodesInPreorder()
	 * @see     TreeNode
	 * @see     NameComparator
	 * @param jpb Progress bar.
	 */
	private void linkLeaves() {
		//int counter = 0;
		//int percentage = 0;
		TreeNode pren = root.leftmostLeaf;
		Vector<TreeNode> leaves = new Vector<TreeNode>();
		leaves.add(pren);
//		pren.lindex = 0;
		for(TreeNode n = pren.preorderNext; n!=null; n=n.preorderNext)
		{
			//counter++;
			if(n.isLeaf())
			{
				leaves.add(n);
			}
		}
		numLeaves = leaves.size();

		NameComparator myNameComparator = new NameComparator();
		TreeNode[] sortedLeafArray = (TreeNode[])leaves.toArray(new TreeNode[leaves.size()]);
		Arrays.sort(sortedLeafArray, myNameComparator);
		int index = 0;
		TreeNode curr = sortedLeafArray[0];
		TreeNode next;
		for(int i=0; i<leaves.size()-1; i++){
			next = sortedLeafArray[i+1]; // only 1 index lookup per iteration
			boolean compare = myNameComparator.compare(curr, next) == 0; 
			if (compare || index > 0)
			{
				String name = curr.getName();
				nodesByName.remove(curr); // before all nodes with
				// same name were being ignored in search and comparing two identically named
				// leaves was broken, much fewer differences in trees with many leaves that
				// have the same name (imagine: all index.html occurences being marked as
				// different since numbering convention doesn't string match the original node name)
				curr.setName(name+ " " + index); //sb.toString());
				nodesByName.put(name+ " " + index, curr);//sortedLeafArray[i].getName(), sortedLeafArray[i]); // add the node back with number convention
				if (!compare)
					index = 0;
				else
					index++;
			}
			curr = next;
		}
	}

	/** Get the leaf associated with the given leaf index.
	 * @param index A leaf index of interest.
	 * @return The leaf node at the index, or null on error.
	 * */
	public TreeNode getLeaf(int index)
	{
//		System.out.println("getting leaf: " + index );
		return null;
	}

	/** Stub function */
	public float getMinObjectValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	/** Stub function */
	public float getMaxObjectValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	/**
	 * @return Returns the index number of this tree in the nexus file it was found in.
	 * Indexing starts at 0, and 0 for non-nexus.  
	 */
	public int getNexusIndex()
	{
		return nexusIndex;
	}
	/**
	 * @param leafSplitAxis The leafSplitLine to set.
	 */
	//public void setLeafSplitAxis(StaticSplitAxis leafSplitAxis) {
	//	this.leafSplitAxis = leafSplitAxis;
	//}
	
	/**
	 * Get the leaves under this node.  Used for tree to tree comparison, removing leaf nodes from difference calculations when they only appear in one side of the tree.
	 * This operation is constant time per leaf, since it relies on pre-ordered node links and pointers to extreme leaves.
	 * Time complexity of this function is linear in the number of leaves in the subtree under the node.
	 * @param node Node to get leaves under.  The root node will return all leaves in the tree, leaves return a list of just themselves.
	 * @return List of leaves under this node.
	 */
	public LinkedList<TreeNode> getLeaves(TreeNode node)
	{
		LinkedList<TreeNode> leaves = new LinkedList<TreeNode>();
		TreeNode currNode = node.leftmostLeaf;
		while (currNode != node.rightmostLeaf)
		{
			if (!currNode.isLeaf()) // internal node?
				currNode = currNode.leftmostLeaf; // descend to minimal leaf
			leaves.add(currNode);
			currNode = currNode.preorderNext;
		}
		leaves.add(node.rightmostLeaf);
		return leaves;
	}
}

/** Comparator class for Strings */
class NameComparator implements Comparator<TreeNode>{
	/** collator object used for string comparison. */
	Collator myCollator = Collator.getInstance(Locale.US);

	/** String comparator, uses {@link Collator} comparator. */
	public int compare(TreeNode o1, TreeNode o2){
		String s1 = (o1).getName();
		String s2 = (o2).getName();
		return myCollator.compare(s1, s2);
	}

}

