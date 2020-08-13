package util.phylogenetics;

import java.io.*;
import java.util.*;
import java.util.Map.*;



import util.Fasta;

public class AnnotateTree {
	
	int strains;
	HashMap<String,Score> Annotation;
	HashMap<String,Score> AnnotationSim;

	ArrayList<Fasta> alignment;
	Tree tree;
	
	public AnnotateTree(File treeFile,HashMap<String,Score> annotationReal,HashMap<String,Score> annotationSim,ArrayList<Fasta> Alignment){
		Annotation=annotationReal;
		AnnotationSim=annotationSim;
		alignment=Alignment;
		strains=alignment.size();
		annotateTree(treeFile);
		System.out.println(AnnotationSim.size());
	}
	
	public void annotateTree(File treeFile){
		try{
			BufferedReader br=new BufferedReader(new FileReader(treeFile));
			TreeParser tp=new TreeParser(br);
			tree=tp.tokenize(1, "whater", null);
			ArrayList<Integer> info=new ArrayList<Integer>();
			ArrayList<Integer> total=new ArrayList<Integer>();
			recursive_simulateAlignment(tree, 0, 0,info,total);
			//System.out.println(info.size()+" of "+total.size());
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private HashMap<String,Score> setDefaultAnnotationSim(){
		HashMap<String,Score> annotationSim=new HashMap<String, Score>();
		Iterator<Entry<String, Score>> it=Annotation.entrySet().iterator();
		while(it.hasNext()){
			Entry<String,Score> e=it.next();
			annotationSim.put(e.getKey(),new Score(1.0,1,1,e.getKey()));
			
		}
		return annotationSim;
	}
	
	public AnnotateTree(File treeFile,HashMap<String,Score> annotation,ArrayList<Fasta> Alignment){
		Annotation=annotation;
		AnnotationSim=setDefaultAnnotationSim();
		alignment=Alignment;
		strains=alignment.size();
		annotateTree(treeFile);
	}
	
	public void writeTree(File out){
		TreeWriter.writeTree(tree, out.toString(),false);
	}
	 void recursive_simulateAlignment (Tree tree,int currkey, int currdepth,ArrayList<Integer> info,ArrayList<Integer> total) {
        TreeNode currNode = tree.getNodeByKey(currkey);
        int numChildren = currNode.numberChildren();
        if(currNode.getName().equals("")){
        	LinkedList<TreeNode> ll=tree.getLeaves(currNode);
        	
        	String pattern;
        	if(ll.size()<=strains/2){
        		pattern=getPattern(ll,'1','0');
        	}else{
            	pattern=getPattern(ll,'0','1');
        	}
        	Score score=Annotation.get(pattern);
        	Score sim=AnnotationSim.get(pattern);
        	
        	if(score==null){
        		score=Annotation.get(inverse(pattern));
        	}
        	if(sim==null){
        		sim=AnnotationSim.get(inverse(pattern));
        	}
        	if(score!=null&&sim!=null){
        		String add=((int)(score.numberSites*1000.0/sim.numberSites)/1000.0)<0.01?((int)(score.numberSites*1000.0/sim.numberSites)/1000.0)+"":"";
        		currNode.setName(((int)((score.score/sim.score)*100)/100.0)+" "+add);
        		info.add(0);
        	}else if(score!=null){
        		currNode.setName(score.score+" noSim "+score.numberSites);
        		info.add(0);
        	}if(score==null&&sim!=null){
        		String add=((int)(1*1000.0/sim.numberSites)/1000.0)<0.01?((int)(1*1000.0/sim.numberSites)/1000.0)+"":"";
        		currNode.setName("noScore "+add);

        	}
        	total.add(0);
        }

        for (int i = 0; i < numChildren; i++) {
            int childkey = currNode.getChild(i).key;
           // TreeNode childnode = tree.getNodeByKey(childkey);
            //String name=childnode.getName();
            //System.out.println(name);
            recursive_simulateAlignment(tree,childkey, currdepth+1,info,total);
        }
    }
	 private String inverse(String pattern){
		 StringBuffer sb=new StringBuffer();
		 for(int i=0;i<pattern.length();i++){
			char a=pattern.charAt(i)=='1'?'0':'1';
			sb.append(a);
		 }
		 return sb.toString();
	 }
	 private String getPattern(LinkedList<TreeNode> idents,char one,char two){
		 StringBuffer sb=new StringBuffer();
		 HashMap<String,Boolean> hm=new HashMap<String, Boolean>();
		 for(int i=0;i<idents.size();i++){
			 hm.put(((TreeNode)idents.get(i)).getName(),true);
		 }

		 for(int i=0;i<alignment.size();i++){
			 if(hm.containsKey(alignment.get(i).getIdent())){
				 sb.append(one);
			 }else{
				 sb.append(two);
			 }
		 }
		
		 return sb.toString();
	 }

}
