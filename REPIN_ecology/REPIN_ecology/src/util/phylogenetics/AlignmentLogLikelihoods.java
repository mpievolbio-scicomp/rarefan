package util.phylogenetics;

import java.io.File;
import java.util.*;

import util.Fasta;





public class AlignmentLogLikelihoods {
	Tree tree;
	double[] pi=new double[4];
	int A=0;
	int T=1;
	int C=2;
	int G=3;
	Alignment alignment;
	HashMap<String,Integer> sitePatterns=new HashMap<String,Integer>();
	ArrayList<String> sitePatternList=new ArrayList<String>();
	ArrayList<Double> sitePatternLikelihoods=new ArrayList<Double>();
	double overallLogLike=0;
	double scalingFactor=1;

	public static void main(String args[]){
		Tree t=Phylogeny.readTree(new File("/home/frederic/Basel/randomAlignments/problematicTS8_-1_3/trial_71_-1_3_0/sequences_7/optimizedREALFull.tree"));
		Alignment alg=new Alignment(Fasta.readPhylip(new File("/home/frederic/Basel/randomAlignments/problematicTS8_-1_3/trial_71_-1_3_0/sequences_7/alignment.phy")));
		AlignmentLogLikelihoods allh=new AlignmentLogLikelihoods(t, alg);
		System.out.println(allh.getTreeLogLike());
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
	
	
	public AlignmentLogLikelihoods(Tree phylogeny,Alignment alg){
		tree=phylogeny;
		alignment=alg;
		//checkOrder();
		initPi();
		initSitePatterns();
		
		scalingFactor=calcScalingFactor();

	}
	
	public void setScalingFactor(double s){
		scalingFactor=s;
	}
	private double calcScalingFactor(){
		double sum=1;
		for(int i=0;i<pi.length;i++){
			sum-=Math.pow(pi[i],2);
		}
		return 1/sum;
	}
	
	public void setTree(Tree t){
		overallLogLike=0;
		tree=t;
		
	}
	
	public double changeTree(Tree t){
		tree = t;
		sitePatternLikelihoods=new ArrayList<Double>();
		calculatesitePatternLikelihoods();
		calculateTreeLogLike();
		return overallLogLike;
	}
	
//	public AlignmentLogLikelihoods(Tree phylogeny,double[] pi){
//		tree=phylogeny;
//		initPi(pi);
//	}
	
	public double getTreeLogLike(){
		if(overallLogLike==0){
			calculatesitePatternLikelihoods();
			calculateTreeLogLike();
		}
		return overallLogLike;
	}
	
	public ArrayList<String> getSitePatterns(){
		return sitePatternList;
	}
	
	public ArrayList<Double> getSitePatternLikelihood(){
		if(sitePatternLikelihoods.size()==0){
			calculatesitePatternLikelihoods();
		}
		return sitePatternLikelihoods;
	}
	
	private void calculateTreeLogLike(){
		double sum=0.0;
		for(int i=0;i<sitePatternLikelihoods.size();i++){
			sum+=sitePatterns.get(sitePatternList.get(i))*Math.log(sitePatternLikelihoods.get(i));
		}
		overallLogLike=sum;
		
		
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
	
//	private void initPi(double[] Pi){
//		pi=List_Array.toDouble(Pi);
//	}
	
	private void initPi(){
		for(int i=0;i<4;i++){
			pi[i]=0.0;
		}
		int count=0;
		for(int i=0;i<alignment.numSeqs;i++){
			for(int j=0;j<alignment.length;j++){
				count++;
				if(alignment.getBase(i, j)=='A'){
					pi[A]=pi[A]+1;
				}else if(alignment.getBase(i, j)=='T'){
					pi[T]=pi[T]+1;
				}else if(alignment.getBase(i, j)=='C'){
					pi[C]=pi[C]+1;
				}else if(alignment.getBase(i, j)=='G'){
					pi[G]=pi[G]+1;
				}	
			}
		}
		for(int i=0;i<4;i++){
			pi[i]=pi[i]/count;
		}
	}
	

	
	private void calculatesitePatternLikelihoods(){
		for(int i=0;i<sitePatternList.size();i++){
			sitePatternLikelihoods.add(getLogLike(sitePatternList.get(i)));
		}
	}
	

	private double getLogLike(String sitePattern){
		HashMap<String,Character> initStates=initStates(sitePattern);
		double[] loglike=LogLike(tree.getRoot().getKey(),initStates);
		double sum=0;
		for(int i=0;i<4;i++){
			sum+=loglike[i]*pi[i];
		}
		return sum;
	}
	

	
	public double[] getLogLikeSites(String sitePattern,int nodekey){
		HashMap<String,Character> initStates=initStates(sitePattern);
		double[] loglike=LogLike(nodekey,initStates);
		return loglike;
	}
	

	
	public double[] getPi(){
		return pi;
	}
	
	public double getLogLike(String sitePattern,int nodekey){
		HashMap<String,Character> initStates=initStates(sitePattern);
		double[] loglike=LogLike(nodekey,initStates);
		double sum=0;
		for(int i=0;i<4;i++){
			sum+=loglike[i]*pi[i];
		}
		return sum;
	}
	
	private HashMap<String,Character> initStates(String siteClass){
		HashMap<String,Character> states=new HashMap<String, Character>();
		for(int i=0;i<alignment.numSeqs;i++){
			states.put(alignment.getIdent(i),siteClass.charAt(i));
		}
		return states;
	}
	
	private int baseToInt(char c){
		if(c=='A'){
			return A;
		}
		if(c=='T'){
			return T;
		}
		if(c=='C'){
			return C;
		}
		if(c=='G'){
			return G;
		}
		return -1;
	}
	
	

	
	private double[] LogLike(int currkey,HashMap<String,Character> initStates) {
		TreeNode currNode = tree.getNodeByKey(currkey);
       //System.out.println(currNode.key+" "+currNode.weight+" "+currNode.name);

		int numChildren = currNode.numberChildren();

		double[] loglike=initLogLike(0.0);
		if(numChildren==0){
			loglike[baseToInt(initStates.get(currNode.getName()))]=1.0;
		}else{
			loglike=initLogLike(1.0);
			for (int i = 0; i < numChildren; i++) {
				int childkey = currNode.getChild(i).key;
				double[] llchild=LogLike(childkey,initStates);
				//System.out.println(llchild);
				double rate=tree.getNodeByKey(childkey).getWeight();
				double eut=Math.pow(Math.E,-rate*scalingFactor);
				multiply(loglike,llchild,eut);
			}
			//System.out.println(loglike);
		}
		return loglike;
	}
	 

	 
	 private void multiply(double[] loglike,double[] child,double eut){
		
		 for(int i=0;i<4;i++){
			 double sum=0;
			 for(int j=0;j<4;j++){
				 sum+=child[j]*P(i,j,eut);
			 }
			 loglike[i]=loglike[i]*sum;
		 }
	 }
	 
	 private double P(int i,int j,double eut){
		 double delta=i!=j?0:1;
		 
		 double p=eut*delta+(1-eut)*pi[j];
		 return p;
	 }
	 
	 private double[] initLogLike(double init){
		 double[] loglike=new double[4];
		 for(int i=0;i<4;i++){
			 loglike[i]=init;
		 }
		 return loglike;
	 }

}
