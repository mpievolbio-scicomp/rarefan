package util.phylogenetics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;

import util.Fasta;

public class RejectMethods {
	Tree tree;
	ArrayList<Fasta> fasta;
	private void divideGroups(int currkey,HashMap<Integer,Double> loglikes){
		TreeNode currNode = tree.getNodeByKey(currkey);
		int numChildren = currNode.numberChildren();
		
		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;
			loglikes.put(childkey, getLogLike(childkey));
			divideGroups(childkey,loglikes);
		}
	}
	private double getLogLike(int node){
		//System.out.println(node);
		TreeNode tn=tree.getNodeByKey(node);
		Tree t=new Tree(tree);

		LinkedList<TreeNode> leavesBelow=t.getLeaves(tn);
		Alignment algBelow=getAlignment(leavesBelow);
		Tree tBelow=t.getSubTree(node);
		AlignmentLogLikelihoods allhBelow=new AlignmentLogLikelihoods(tBelow,algBelow);

		t=t.rerootTree(node);

		Tree tAbove=t.getSubTree(t.getRoot().getKey());
		LinkedList<TreeNode> leavesAbove=tAbove.getLeaves(tAbove.getRoot());
		Alignment algAbove=getAlignment(leavesAbove);

		AlignmentLogLikelihoods allhAbove=new AlignmentLogLikelihoods(tAbove,algAbove);	
//		System.out.println("leaves "+tn.getName());
//		printLeaves(tBelow.getLeaves(tBelow.getRoot()));
//
//		System.out.println("alg");
//		printIdents(algBelow);
//		System.out.println("tree");
//		parseTree(tBelow,tBelow.getRoot().getKey());
//		System.out.println("end");
		
		return allhAbove.getTreeLogLike()
				+allhBelow.getTreeLogLike();
	}
	
	private Alignment getAlignment(LinkedList<TreeNode> leaves){
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
		for(int i=0;i<leaves.size();i++){
			int name=Integer.parseInt(leaves.get(i).getName().substring(1));
			fas.add(this.fasta.get(name));
		}
		return new Alignment(fas);
	}
	public void determineRoot(){
		HashMap<Integer,Double> loglikes=new HashMap<Integer, Double>();
		divideGroups(tree.getRoot().getKey(),loglikes);
		int key=getMax(loglikes);
		System.out.println(tree.getNodeByKey(key).getName());

		tree=tree.rerootTree(key);
	}
	
	private int getMax(HashMap<Integer,Double> ll){
		Iterator<Entry<Integer,Double>> it=ll.entrySet().iterator();
		int maxIndex=0;
		double maxLL=Double.NEGATIVE_INFINITY;
		while(it.hasNext()){
			Entry<Integer,Double> e=it.next();
			System.out.println(e.getValue()+" "+e.getKey()+" "+tree.getNodeByKey(e.getKey()).getName());
			if(e.getValue()>maxLL){
				maxLL=e.getValue();
				maxIndex=e.getKey();
			}
		}
		return maxIndex;
	}
}
