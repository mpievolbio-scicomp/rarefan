package util.phylogenetics;

import java.io.File;
import java.util.*;

import util.Fasta;


public class OptimizeLogLike {

	public static void main(String args[]){
		File fas=new File(args[0]);
		File tree=new File(args[1]);
		double precision=Double.parseDouble(args[2]);
		String fasName=fas.getName().split("\\.")[0];
		File out=new File(tree.getParent()+"/"+fasName+"optimized.tree");
		Alignment alg=new Alignment(Fasta.readPhylip(fas));
		Tree topo=Phylogeny.readTree(tree);
		OptimizeLogLike oll=new OptimizeLogLike(topo, alg, precision, 0.5 );
		oll.writeOptimizedTree(out);
	}
	
	Tree tree;
	double[] pi=new double[4];
	int A=0;
	int T=1;
	int C=2;
	int G=3;
	Alignment alignment;

	
	HashMap<String,Integer> sitePatterns=new HashMap<String,Integer>();
	ArrayList<String> sitePatternList=new ArrayList<String>();
	ArrayList<Double> siteClassLogLikes=new ArrayList<Double>();
	int numSitePatterns;
	int K;
	double precision;
	AlignmentLogLikelihoods alllh;
	double scalingFactor=1;
	double likelihood=0;
	public double getLogLikelihood(){
		return likelihood;
	}
	private void checkOrder(){
		String[] order=new String[]{"S11","S12","S21","S22"};
		for(int i=0;i<alignment.numSeqs;i++){
			if(!order[i].equals(alignment.idents.get(i))){
				System.err.println("PROBLEM");
				System.exit(-1);
			}
		}
	}
	public OptimizeLogLike(Tree topology,Alignment alg,double Precision,double initRate){
		//System.err.println("Only tested for four leaf trees!");
		precision=Precision;
		alignment=alg;
		checkOrder();
		K=alignment.length;
		initPi();
		initSitePatterns();
		scalingFactor=calcScalingFactor();
		numSitePatterns=sitePatternList.size();
		tree=new Tree(topology);
		initTree(initRate,0);
		alllh=new AlignmentLogLikelihoods(tree, alg);
		double llold=alllh.getTreeLogLike();
		
		double llnew=optimizeLogLike(0);
	//	System.out.println(llnew+" "+llold);
		while(Math.abs(llold-llnew)>precision){
			llold=llnew;
			llnew=optimizeLogLike(0);
			
		//	System.out.println(llnew);
		}
		likelihood=llnew;
	}
	
	public Tree getTree(){
		return tree;
	}
	private double calcScalingFactor(){
		double sum=1;
		for(int i=0;i<pi.length;i++){
			sum-=Math.pow(pi[i],2);
		}
		return 1/sum;
	}
	public void writeOptimizedTree(File out){
		TreeWriter.writeTree(tree, out.toString(),false);
	}
	
	public void initTree(double initRate,int currkey){
		TreeNode currNode = tree.getNodeByKey(currkey);
		currNode.setWeight(initRate);
		int numChildren = currNode.numberChildren();

		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;
			initTree(initRate,childkey);
		}
	}
	
	private void initSitePatterns(){
		
		for(int i=0;i<alignment.length;i++){
			String col=alignment.getColumn(i);
			if(!sitePatterns.containsKey(col)){
				sitePatterns.put(col,1);
				sitePatternList.add(col);
			}else{
				sitePatterns.put(col, sitePatterns.get(col)+1);
			}
		}
	}
	
	private void initPi(){
		for(int i=0;i<4;i++){
			pi[i]=0.0;
		}
		int count=0;
		for(int i=0;i<alignment.numSeqs;i++){
			for(int j=0;j<alignment.length;j++){
				count++;
				if(alignment.getBase(i, j)=='A'){
					pi[A]+=1;
				}else if(alignment.getBase(i, j)=='T'){
					pi[T]+=1;
				}else if(alignment.getBase(i, j)=='C'){
					pi[C]+=1;
				}else if(alignment.getBase(i, j)=='G'){
					pi[G]+=1;
				}	
			}
		}
		for(int i=0;i<4;i++){
			pi[i]= pi[i]/count;
		}
	}

	private double[] getAlist(int node){
		double[] A=new double[numSitePatterns];
		for(int i=0;i<numSitePatterns;i++){
			double loglike=0;
			alllh.setTree(tree);
			double[] ll1=alllh.getLogLikeSites(sitePatternList.get(i), node);
			Tree reroot=tree.rerootTree(node);
			alllh.setTree(reroot);
			double[] ll2=alllh.getLogLikeSites(sitePatternList.get(i), reroot.getRoot().getKey());
			for(int j=0;j<4;j++){
				
				double mult1=pi[j]*ll1[j];
				
				
				double mult2=ll2[j];
				loglike+=mult1*mult2;
			}
			A[i]=loglike;
		}
		return A;
	}
	
	private double calcLogLikeB(double[] loglikes){
		double loglike=0;
		for(int j=0;j<4;j++){
			loglike+=loglikes[j]*pi[j];
		}
		return loglike;
	}
	
	private double[] getBlist(int node){
		double[] B=new double[numSitePatterns];
		for(int i=0;i<numSitePatterns;i++){
			alllh.setTree(tree);
			double loglike1=calcLogLikeB(alllh.getLogLikeSites(sitePatternList.get(i), node));
			Tree reroot=tree.rerootTree(node);
			alllh.setTree(reroot);
			double loglike2=calcLogLikeB(alllh.getLogLikeSites(sitePatternList.get(i), reroot.getRoot().getKey()));
			
			B[i]=loglike1*loglike2;
		}
		return B;
	}
	private double getLogLike(double[] A,double[] B,double rate){
		double sum=0;
		double q=Math.pow(Math.E,-rate*scalingFactor);
		double p=1-q;
		for(int i=0;i<numSitePatterns;i++){
			sum+=sitePatterns.get(sitePatternList.get(i))*Math.log(A[i]*q+B[i]*p);
		}
		return sum;
	}
	private double optimizeLogLike(int currkey) {
		TreeNode currNode = tree.getNodeByKey(currkey);
		int numChildren = currNode.numberChildren();
		double loglike=0;
		if(currNode.getKey()!=0){
			double[] Alist=getAlist(currNode.getKey());
			double[] Blist=getBlist(currNode.getKey());
			double rate=currNode.getWeight();
			double ratenew=iterate(rate,Alist,Blist);
			currNode.setWeight(ratenew);
			loglike=getLogLike(Alist,Blist,ratenew);
		}else{
			currNode.setWeight(0);
		}
		for (int i = 0; i < numChildren; i++) {
			int childkey = currNode.getChild(i).key;        
			loglike=optimizeLogLike(childkey);
		}
		return loglike;
	}
	
	private double calcP(double p,double q,double[] A,double[] B){
		double pnew=0;
		
		for(int i=0;i<numSitePatterns;i++){
			String sp=sitePatternList.get(i);
			pnew+=sitePatterns.get(sp)*(B[i]*p/(A[i]*q+B[i]*p));
		}
		pnew=pnew/K;
		return pnew;
	}
	 
	public double iterate(double rate,double[] A,double[] B){
		double q=Math.pow(Math.E,-rate*scalingFactor);
		double p=1-q;
		double pnew=calcP(p,q, A, B);

		while(Math.abs(pnew-p)>precision){
			p=calcP(pnew,1-pnew, A, B);
			double s=p;
			pnew=p;
			p=s;
		}
		
		return -Math.log(1-pnew)/scalingFactor;
	}
}
		
		


		

		
