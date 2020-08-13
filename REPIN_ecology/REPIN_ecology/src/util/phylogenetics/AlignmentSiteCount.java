package util.phylogenetics;

import java.io.*;
import java.util.*;

import statistics.Stats;
import util.*;
import util.phylogenetics.*;

public class AlignmentSiteCount {
	
	public static void main(String args[])throws IOException{
		File alignmentFolder=new File(args[0]);
		File treedistpath=new File(args[1]);
		String treeshape=args[2];
		File Rscript=new File(args[3]);
		File RAxMLPath=new File(args[4]);

		System.setOut(new PrintStream(new OutputStream(){
			public void write(int b) {
			}
		}));
		String[] order=new String[]{"S11","S12","S21","S22"};
		File[] files=alignmentFolder.listFiles();
		ArrayList<Tree> treeCorrectRPML=new ArrayList<Tree>();
		ArrayList<Tree> treeIncorrect1RPML=new ArrayList<Tree>();
		ArrayList<Tree> treeIncorrect2RPML=new ArrayList<Tree>();


		
		ArrayList<File> treeCorrectRPRAxML=new ArrayList<File>();
		ArrayList<File> treeIncorrect1RPRAxML=new ArrayList<File>();
		ArrayList<File> treeIncorrect2RPRAxML=new ArrayList<File>();

		
		ArrayList<File> refAlignments=new ArrayList<File>();
		ArrayList<File> fullAlignments=new ArrayList<File>();
		
		ArrayList<Tree> treeCorrectFullML=new ArrayList<Tree>();
		ArrayList<Tree> treeIncorrect1FullML=new ArrayList<Tree>();
		ArrayList<Tree> treeIncorrect2FullML=new ArrayList<Tree>();
		
		ArrayList<File> treeCorrectFullRAxML=new ArrayList<File>();
		ArrayList<File> treeIncorrect1FullRAxML=new ArrayList<File>();
		ArrayList<File> treeIncorrect2FullRAxML=new ArrayList<File>();

		
		ArrayList<AlignmentSiteCount> exact=new ArrayList<AlignmentSiteCount>();
		ArrayList<AlignmentSiteCount> realphy=new ArrayList<AlignmentSiteCount>();
		File realTree=Phylogeny.makeRealTree(alignmentFolder);
		File inc1=Phylogeny.makeIncorrectTree1(alignmentFolder);
		File inc2=Phylogeny.makeIncorrectTree2(alignmentFolder);
		Tree realT=Phylogeny.readTree(realTree);
		Tree inc1T=Phylogeny.readTree(inc1);
		Tree inc2T=Phylogeny.readTree(inc2);
		ArrayList<Integer> bestScoreRAxML=new ArrayList<Integer>();
		ArrayList<Integer> bestScoreML=new ArrayList<Integer>();

		for(int i=0;i<files.length;i++){
			if(files[i].getName().startsWith("trial")){
				File alignmentRef=new File(files[i]+"/S11_"+treeshape+"/PolySeqOut_NoGenes/polymorphisms_move.phy");
				File alignmentFull=new File(files[i]+"/sequences_"+treeshape+"/alignment.phy");
				//File tree=new File(files[i]+"/S11_"+treeshape+"/PolySeqOut_NoGenes/RAxML_bestTree.raxml");
				ArrayList<Fasta> fas=Fasta.readPhylip(alignmentRef,order);
				ArrayList<Fasta> fasFull=Fasta.readPhylip(alignmentFull, order);
				AlignmentSiteCount alscRef=new AlignmentSiteCount(fas);
				AlignmentSiteCount alscFull=new AlignmentSiteCount(fasFull);
				realphy.add(alscRef);
				exact.add(alscFull);
				File alignmentRefOrder=new File(files[i]+"/S11_"+treeshape+"/PolySeqOut_NoGenes/polymorphisms_move_order.phy");
				File alignmentFullOrder=new File(files[i]+"/sequences_"+treeshape+"/alignment_order.phy");

				bestScoreRAxML.add(getBestTree(realTree,inc1,inc2,new File(files[i]+"/S11_"+treeshape+"/PolySeqOut_NoGenes/RAxML_bestTree.raxml"),treedistpath));
				
				Fasta.writePhylip(fas, alignmentRefOrder, 100);
				Fasta.writePhylip(fasFull, alignmentFullOrder, 100);
				refAlignments.add(alignmentRefOrder);
				fullAlignments.add(alignmentFullOrder);
				File rpFolder=new File(files[i]+"/S11_"+treeshape+"/PolySeqOut_NoGenes");
				File fullFolder=new File(files[i]+"/sequences_"+treeshape+"/");
				Tree correctRPML=getOptimizedTree(realT,fas,new File(rpFolder+"/optimizedREALRP.tree"),false);
				treeCorrectRPML.add(correctRPML);
				Tree inc1RPML=getOptimizedTree(inc1T,fas,new File(rpFolder+"/optimizedInc1RP.tree"),false);
				treeIncorrect1RPML.add(inc1RPML);
				Tree inc2RPML=getOptimizedTree(inc2T,fas,new File(rpFolder+"/optimizedInc2RP.tree"),false);
				treeIncorrect2RPML.add(inc2RPML);
				bestScoreML.add(getBestTree(correctRPML,inc1RPML,inc2RPML,new Alignment(fas)));				
				
				treeCorrectFullML.add(getOptimizedTree(realT,fasFull,new File(fullFolder+"/optimizedREALFull.tree"),false));				
				treeIncorrect1FullML.add(getOptimizedTree(inc1T,fasFull,new File(fullFolder+"/optimizedInc1Full.tree"),false));
				treeIncorrect2FullML.add(getOptimizedTree(inc2T,fasFull,new File(fullFolder+"/optimizedInc2Full.tree"),false));		
				
				treeCorrectRPRAxML.add(RunTreePrograms.runRAxMLGivenTopo(realTree,alignmentRefOrder, RAxMLPath,"realRP", "S11", 1234,false,-1));
				treeIncorrect1RPRAxML.add(RunTreePrograms.runRAxMLGivenTopo(inc1,alignmentRefOrder, RAxMLPath, "inc1RP", "S11", 1234,false,-1));
				treeIncorrect2RPRAxML.add(RunTreePrograms.runRAxMLGivenTopo(inc2,alignmentRefOrder, RAxMLPath, "inc2RP", "S11", 1234,false,-1));
				
				treeCorrectFullRAxML.add(RunTreePrograms.runRAxMLGivenTopo(realTree,alignmentFullOrder, RAxMLPath,"realFull", "S11", 1234,false,-1));
				treeIncorrect1FullRAxML.add(RunTreePrograms.runRAxMLGivenTopo(inc1,alignmentFullOrder, RAxMLPath, "inc1Full", "S11", 1234,false,-1));
				treeIncorrect2FullRAxML.add(RunTreePrograms.runRAxMLGivenTopo(inc2,alignmentFullOrder, RAxMLPath, "inc2Full", "S11", 1234,false,-1));
				

			}
		}
		//first level 15 site classes
		//second level 100 trials
		ArrayList<ArrayList<LogLikeDiff>> loglikeRAxMLInc1=initloglike();
		ArrayList<ArrayList<LogLikeDiff>> loglikeRAxMLInc2=initloglike();
		
		ArrayList<ArrayList<LogLikeDiff>> loglikeMLInc1=initloglike();
		ArrayList<ArrayList<LogLikeDiff>> loglikeMLInc2=initloglike();
		
		writeSiteClasses(alignmentFolder,realphy);

		ArrayList<Integer> correctRAxML=new ArrayList<Integer>();
		ArrayList<Integer> inc1RAxML=new ArrayList<Integer>();
		ArrayList<Integer> inc2RAxML=new ArrayList<Integer>();
		
		ArrayList<Integer> correctML=new ArrayList<Integer>();
		ArrayList<Integer> inc1ML=new ArrayList<Integer>();
		ArrayList<Integer> inc2ML=new ArrayList<Integer>();

		//for(int i=0;i<treeCorrectRPRAxML.size();i++){
		for(int i=0;i<treeCorrectFullML.size();i++){
			System.out.println("Correct tree "+i+" of "+treeCorrectRPRAxML.size());
			
			HashMap<String,Double> freqTreeCorrectRAxMLRP=SiteClasses.getLogLikelihoodsRAxML(treeCorrectRPRAxML.get(i),refAlignments.get(i), order,RAxMLPath,"correct",false);
			HashMap<String,Double> freqTreeIncorrect1RAxMLRP=SiteClasses.getLogLikelihoodsRAxML(treeIncorrect1RPRAxML.get(i),refAlignments.get(i), order,RAxMLPath,"inc1",false);
			HashMap<String,Double> freqTreeIncorrect2RAxMLRP=SiteClasses.getLogLikelihoodsRAxML(treeIncorrect2RPRAxML.get(i),refAlignments.get(i), order,RAxMLPath,"inc2",false);
			
			fillIndexList(bestScoreRAxML,correctRAxML,inc1RAxML,inc2RAxML,i);
			fillIndexList(bestScoreML,correctML,inc1ML,inc2ML,i);

			HashMap<String,Double> freqTreeCorrectRAxMLFull=SiteClasses.getLogLikelihoodsRAxML(treeCorrectFullRAxML.get(i),fullAlignments.get(i), order,RAxMLPath,"correct",false);
			HashMap<String,Double> freqTreeIncorrect1RAxMLFull=SiteClasses.getLogLikelihoodsRAxML(treeIncorrect1FullRAxML.get(i),fullAlignments.get(i), order,RAxMLPath,"inc1",false);
			HashMap<String,Double> freqTreeIncorrect2RAxMLFull=SiteClasses.getLogLikelihoodsRAxML(treeIncorrect2FullRAxML.get(i),fullAlignments.get(i), order,RAxMLPath,"inc2",false);
			
			
			HashMap<String,Double> freqTreeCorrectMLRP=SiteClasses.getLogLikelihoodsML(treeCorrectRPML.get(i),refAlignments.get(i), order);			
			HashMap<String,Double> freqTreeIncorrect1MLRP=SiteClasses.getLogLikelihoodsML(treeIncorrect1RPML.get(i),refAlignments.get(i), order);
			HashMap<String,Double> freqTreeIncorrect2MLRP=SiteClasses.getLogLikelihoodsML(treeIncorrect2RPML.get(i),refAlignments.get(i), order);
			
			HashMap<String,Double> freqTreeCorrectMLFull=SiteClasses.getLogLikelihoodsML(treeCorrectFullML.get(i),fullAlignments.get(i), order);			
			HashMap<String,Double> freqTreeIncorrect1MLFull=SiteClasses.getLogLikelihoodsML(treeIncorrect1FullML.get(i),fullAlignments.get(i), order);
			HashMap<String,Double> freqTreeIncorrect2MLFull=SiteClasses.getLogLikelihoodsML(treeIncorrect2FullML.get(i),fullAlignments.get(i), order);
						
			add(loglikeMLInc1,getLogLikelihoods(exact.get(i).getSiteClassFreqs(),realphy.get(i).getSiteClassFreqs(),freqTreeCorrectMLRP,freqTreeIncorrect1MLRP,freqTreeCorrectMLFull,freqTreeIncorrect1MLFull));
			add(loglikeMLInc2,getLogLikelihoods(exact.get(i).getSiteClassFreqs(),realphy.get(i).getSiteClassFreqs(),freqTreeCorrectMLRP,freqTreeIncorrect2MLRP,freqTreeCorrectMLFull,freqTreeIncorrect2MLFull));

			add(loglikeRAxMLInc1,getLogLikelihoods(exact.get(i).getSiteClassFreqs(),realphy.get(i).getSiteClassFreqs(),freqTreeCorrectRAxMLRP,freqTreeIncorrect1RAxMLRP,freqTreeCorrectRAxMLFull,freqTreeIncorrect1RAxMLFull));
			add(loglikeRAxMLInc2,getLogLikelihoods(exact.get(i).getSiteClassFreqs(),realphy.get(i).getSiteClassFreqs(),freqTreeCorrectRAxMLRP,freqTreeIncorrect2RAxMLRP,freqTreeCorrectRAxMLFull,freqTreeIncorrect2RAxMLFull));
//			if(loglikeMLInc1.get(0).get(i).FullCorr>-0.18){
//				//getLogLikelihoods(exact.get(i).getSiteClassFreqs(),realphy.get(i).getSiteClassFreqs(),freqTreeCorrectRAxMLRP,freqTreeIncorrect1RAxMLRP,freqTreeCorrectRAxMLFull,freqTreeIncorrect1RAxMLFull, true);
//				System.err.println(treeCorrectFullRAxML.get(i));
//			}
		}
		
//		printLoglikeStats(loglikeRAxMLInc1,new File(alignmentFolder+"/loglikeStatsRAxML_inc1.txt"),Rscript);
//		printLoglikeStats(loglikeRAxMLInc2,new File(alignmentFolder+"/loglikeStatsRAxML_inc2.txt"),Rscript);
//
//		
//		printLoglikeStats(loglikeMLInc1,new File(alignmentFolder+"/loglikeStatsML_inc1.txt"),Rscript);
//		printLoglikeStats(loglikeMLInc2,new File(alignmentFolder+"/loglikeStatsML_inc2.txt"),Rscript);


	
	}
	
	private static void fillIndexList(ArrayList<Integer> bestTrees,ArrayList<Integer> correct,ArrayList<Integer> inc1,ArrayList<Integer> inc2,int i){
		if(bestTrees.get(i)==0){
			correct.add(i);
		}else if(bestTrees.get(i)==1){
			inc1.add(i);
		}else if(bestTrees.get(i)==2){
			inc2.add(i);
		}
	}
	
	public static Tree getOptimizedTree(Tree tree,ArrayList<Fasta> fas,File out,boolean forceOptimize){
		if(out.exists()&&!forceOptimize){
			Tree t1= Phylogeny.readTree(out);
//			OptimizeLogLike oll=new OptimizeLogLike(tree,new Alignment(fas),0.1,0.5);
//			Tree direct=oll.getTree();
//			compareTrees(direct,t1,out,new Alignment(fas));
//			System.err.println("______");
			return t1;
		}else{
			OptimizeLogLike oll=new OptimizeLogLike(tree,new Alignment(fas),0.1,0.5);
			oll.writeOptimizedTree(out);
//			Tree direct=oll.getTree();
//			Tree indirect=Phylogeny.readTree(out);
//			compareTrees(direct,indirect,out);
			return oll.getTree();
		}
	}
	static void recursive(Tree t1,int key){
		TreeNode tn1=t1.getNodeByKey(key);
		int numchildren=tn1.numberChildren();
		System.err.println(tn1.getKey()+" "+tn1.getName()+" "+tn1.getWeight());
		for(int i=0;i<numchildren;i++){
			TreeNode child=tn1.getChild(i);
			recursive(t1,child.getKey());
		}
	}
	
	static void recursive(Tree t1,Tree t2,int key,File file){
		TreeNode tn1=t1.getNodeByKey(key);
		TreeNode tn2=t2.getNodeByKey(key);
		int numchildren=tn1.numberChildren();
		if(tn1.getWeight()!=tn2.getWeight()){
			System.err.println("Problem!!! "+tn1.getWeight()+" "+tn2.getWeight()+" "+file);
		}
		if(!tn1.getName().equals(tn2.getName())){
			System.err.println("Problem!!! "+tn1.getName()+" "+tn2.getName()+" "+file);
		}
			//tn.setWeight(key/10.0);
		
		for(int i=0;i<numchildren;i++){
			TreeNode child=tn1.getChild(i);
			recursive(t1,t2,child.getKey(),file);
		}
	}
	private static void compareTrees(Tree t1,Tree t2){
		recursive(t1,t1.getRoot().key);
		recursive(t2,t2.getRoot().key);
//		AlignmentLogLikelihoods allh=new AlignmentLogLikelihoods(t1,alg);
//		System.err.println(allh.getTreeLogLike());
//		allh=new AlignmentLogLikelihoods(t2,alg);
//		System.err.println(allh.getTreeLogLike());
		//print(SiteClasses.getLogLikelihoodsML(t1, alg),SiteClasses.getLogLikelihoodsML(t2, alg));

	}
	private static void print(HashMap<String,Double> sc1,HashMap<String,Double> sc2){
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			if(sc1.get(SiteClasses.siteClasses.get(i)).equals(sc2.get(SiteClasses.siteClasses.get(i)))){
				System.err.println(SiteClasses.siteClasses.get(i)+" "+sc1.get(SiteClasses.siteClasses.get(i))+" "+sc2.get(SiteClasses.siteClasses.get(i)));
			}
		}
	}
	
	
	public static int getBestTree(Tree t1,Tree t2,Tree t3,Alignment alg){
		AlignmentLogLikelihoods all=new AlignmentLogLikelihoods(t1, alg);
		ArrayList<Double> ll=new ArrayList<Double>();
		ll.add(all.getTreeLogLike());
		 all=new AlignmentLogLikelihoods(t2, alg);
		ll.add(all.getTreeLogLike());
		 all=new AlignmentLogLikelihoods(t3, alg);
		ll.add(all.getTreeLogLike());
		int maxIndex=getMaxIndex(ll);
		return maxIndex;
	}
	
	public static int getBestTree(File t1,File t2,File t3,File calc,File treedistPath){
		try{
		if(Phylogeny.compare(t1, calc, treedistPath)==0)return 0;
		if(Phylogeny.compare(t2, calc, treedistPath)==0)return 1;
		if(Phylogeny.compare(t3, calc, treedistPath)==0)return 2;
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return -1;
	}
	
	private static int getMaxIndex(ArrayList<Double> list){
		int maxIndex=-1;
		double max=-Double.MAX_VALUE;
		for(int i=0;i<list.size();i++){
			if(max<list.get(i)){
				max=list.get(i);
				maxIndex=i;
			}
		}
		return maxIndex;
	}
	public static void writeSiteClasses(File alignmentFolder,ArrayList<AlignmentSiteCount> realphy) throws IOException{
		File siteclasses=new File(alignmentFolder+"/siteclassesRealPhy.txt");
		BufferedWriter bw=new BufferedWriter(new FileWriter(siteclasses));
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			String sc=SiteClasses.siteClasses.get(i);
			bw.write(sc);
			for(int j=0;j< realphy.size();j++){
				bw.write("\t"+realphy.get(j).getSiteClassFreqs().get(sc)+"__"+realphy.get(j).numSites);
			}
			bw.write("\n");
		}
		bw.close();
	}
	
	private static HashMap<Integer,Boolean> getHM(ArrayList<Integer> list){
		HashMap<Integer,Boolean> indexHM=new HashMap<Integer, Boolean>();
		for(int i=0;i<list.size();i++){
			indexHM.put(list.get(i), true);
		}
		return indexHM;
	}
	
	public static ArrayList<ArrayList<LogLikeDiff>> getItems(ArrayList<ArrayList<LogLikeDiff>> list,ArrayList<Integer> indeces){
		ArrayList<ArrayList<LogLikeDiff>> newList=new ArrayList<ArrayList<LogLikeDiff>>();
		HashMap<Integer,Boolean> indexHM=getHM(indeces);
		for(int i=0;i<list.size();i++){
			
			newList.add(new ArrayList<LogLikeDiff>());
			for(int j=0;j<list.get(i).size();j++){
				if(indexHM.containsKey(j)){
					newList.get(i).add(list.get(i).get(j));
				}
			}
		}
		return newList;
	}
	
//	public static void printLoglikeStats(ArrayList<ArrayList<LogLikeDiff>> loglike,File out,File scriptPath){
//		try{
//			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
//			File Rout=new File(out.getParent()+"/"+out.getName().split("\\.")[0]+".jpg");
//			File Rcodeout=new File(out.getParent()+"/"+out.getName().split("\\.")[0]+".R");
//			int size=SiteClasses.siteClasses.size();
//			RCode rc=R_functions.plot_InitJpg(Rout,15.0, 7.5);
//			rc.addRCode("par(cex=0.5,mar=c(4,4,1,1)+0.1)");
//			//ArrayList<ArrayList<LogLikeDiff>> loglikeRed=getItems(loglike,red);
//			//ArrayList<ArrayList<LogLikeDiff>> loglikeGreen=getItems(loglike,green);
//			R_functions.produceRcode_BoxplotCategoriesLogLike(rc, loglike,SiteClasses.siteClasses, "","green","red",-0.004,0.005,-0.2,0,-0.2,0);
//			//R_functions.produceRcode_BoxplotCategories(rc, loglike,SiteClasses.siteClasses, "");
//			
//			bw.write("SiteClass\tlikelihoodFullCorr\tstderr\tlikelihoodFullInc\tstderr\tlikelihoodRPCorr\tstderr\tlikelihoodRPInc\tstderr\tTotal\tstderr\n");
//			for(int i=0;i<size;i++){
//				ArrayList<Stats> s=getStats(loglike.get(i));
//				bw.write(SiteClasses.siteClasses.get(i));
//				for(int j=0;j<s.size();j++){
//					bw.write("\t"+s.get(j).getAverage()+"\t"+s.get(j).getStandardError()+"\n");
//				}
//				bw.write("\n");
//			}
//			R_functions.runRCode(rc, scriptPath);
//
//			bw.close();
//			BufferedWriter rcode=new BufferedWriter(new FileWriter(Rcodeout));
//			rcode.write(rc.getCode().toString());
//			rcode.close();
//		}catch(IOException e){
//			e.printStackTrace();
//		}
//	}
	
//	public static void printLoglikeStats(ArrayList<ArrayList<LogLikeDiff>> loglike,ArrayList<Integer> correct,ArrayList<Integer> incorrect,File out,File scriptPath){
//		try{
//			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
//			File Rout=new File(out.getParent()+"/"+out.getName().split("\\.")[0]+".jpg");
//			File Rcodeout=new File(out.getParent()+"/"+out.getName().split("\\.")[0]+".R");
//			int size=SiteClasses.siteClasses.size();
//			RCode rc=R_functions.plot_InitJpg(Rout,15.0, 15.0);
//			rc.addRCode("par(cex=0.5)");
//			//ArrayList<ArrayList<LogLikeDiff>> loglikeRed=getItems(loglike,red);
//			//ArrayList<ArrayList<LogLikeDiff>> loglikeGreen=getItems(loglike,green);
//			R_functions.produceRcode_BoxplotCategoriesLogLike(rc, loglike,correct,incorrect,SiteClasses.siteClasses, "","green","red",-0.004,0.005,-0.15,0,-0.2,0);
//			//R_functions.produceRcode_BoxplotCategories(rc, loglike,SiteClasses.siteClasses, "");
//			
//			bw.write("SiteClass\tlikelihoodFullCorr\tstderr\tlikelihoodFullInc\tstderr\tlikelihoodRPCorr\tstderr\tlikelihoodRPInc\tstderr\tTotal\tstderr\n");
//			for(int i=0;i<size;i++){
//				ArrayList<Stats> s=getStats(loglike.get(i));
//				bw.write(SiteClasses.siteClasses.get(i));
//				for(int j=0;j<s.size();j++){
//					bw.write("\t"+s.get(j).getAverage()+"\t"+s.get(j).getStandardError()+"\n");
//				}
//				bw.write("\n");
//			}
//			R_functions.runRCode(rc, scriptPath);
//
//			bw.close();
//			BufferedWriter rcode=new BufferedWriter(new FileWriter(Rcodeout));
//			rcode.write(rc.getCode().toString());
//			rcode.close();
//		}catch(IOException e){
//			e.printStackTrace();
//		}
//	}
	
	public static ArrayList<Stats> getStats(ArrayList<LogLikeDiff> list){
		ArrayList<Stats> s=new ArrayList<Stats>();
		double[] fullCorr=new double[list.size()];
		double[] fullInc=new double[list.size()];
		double[] rpCorr=new double[list.size()];
		double[] rpInc=new double[list.size()];
		double[] total=new double[list.size()];
		for(int i=0;i<list.size();i++){
			fullCorr[i]=list.get(i).FullCorr;
			fullInc[i]=list.get(i).FullInc;
			rpCorr[i]=list.get(i).RPCorr;
			rpInc[i]=list.get(i).RPInc;
			total[i]=list.get(i).Total;
		}
		s.add(new Stats(fullCorr));
		s.add(new Stats(fullInc));
		s.add(new Stats(rpCorr));
		s.add(new Stats(rpInc));
		s.add(new Stats(total));
		return s;
	}
	
	public static ArrayList<ArrayList<LogLikeDiff>> initloglike(){
		ArrayList<ArrayList<LogLikeDiff>> list=new ArrayList<ArrayList<LogLikeDiff>>();
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			list.add(new ArrayList<LogLikeDiff>());
		}
		return list;
	}
	public static void add(ArrayList<ArrayList<LogLikeDiff>> loglikes,ArrayList<LogLikeDiff> loglikeSiteClasses){
		for(int i=0;i<loglikeSiteClasses.size();i++){
			loglikes.get(i).add(loglikeSiteClasses.get(i));
		}
		
	}
	
	public static HashMap<String,Double> getSiteClassFreqs(File tree, String[] order){
		ArrayList<Fasta> treeFasIncorrect=Fasta.sort(SimulateAlignmentForTree.getAlignment(tree, 1000000, 0.5),order);


		AlignmentSiteCount alscTreeIncorrect=new AlignmentSiteCount(treeFasIncorrect);

		HashMap<String,Double> freqTreeIncorrect=alscTreeIncorrect.siteClassFreqs;

		return freqTreeIncorrect;

	}
	
	public static double getLike(File in){
		double like=0;
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			
			while((line=br.readLine())!=null){
				if(line.startsWith("Final GAMMA-based Score of best tree")){
					like=Double.parseDouble(line.split("\\s+")[6]);
				}
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return like;
	}
	
	public static void printLogLikelihoods(HashMap<String,Double> freqFull,HashMap<String,Double> freqRef,HashMap<String,Double> freqTreeCorrect,HashMap<String,Double> freqTreeIncorrect,File out,int count,boolean append){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out,append));
			bw.write("siteClasses\tfreqFull\tfreqRef\tfreqTreeCorrect\tfreqTreeIncorrect\tlogLike\n");
			for(int i=0;i<SiteClasses.siteClasses.size();i++){
				String sc=SiteClasses.siteClasses.get(i);
				double fFull=freqFull.get(sc)/count;
				double fRef=freqRef.get(sc)/count;
				double ftreeCorrect=freqTreeCorrect.get(sc);
				double ftreeIncorrect=freqTreeIncorrect.get(sc);
				double logLike=Math.log10(ftreeCorrect/ftreeIncorrect)*(fFull-fRef);
				//double correctLogLike=fRef*(Math.log10(ftreeCorrect));
				//double incorrectLogLike=fRef*(Math.log10(ftreeIncorrect));
				bw.write(sc+"\t"+fFull+"\t"+fRef+"\t"+ftreeCorrect+"\t"+ftreeIncorrect+"\t"+logLike+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}

	public static ArrayList<Double> getLogLikelihoods(HashMap<String,Double> freqFull,HashMap<String,Double> freqRef,HashMap<String,Double> freqTreeCorrectRP,HashMap<String,Double> freqTreeIncorrectRP,HashMap<String,Double> freqTreeCorrectFull,HashMap<String,Double> freqTreeIncorrectFull,boolean print){
		ArrayList<Double> loglikes=new ArrayList<Double>();
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			String sc=SiteClasses.siteClasses.get(i);
			double fFull=freqFull.get(sc);
			double fRef=freqRef.get(sc);
			double ftreeCorrectRP=freqTreeCorrectRP.get(sc);
			double ftreeIncorrectRP=freqTreeIncorrectRP.get(sc);
			double ftreeCorrectFull=freqTreeCorrectFull.get(sc);
			double ftreeIncorrectFull=freqTreeIncorrectFull.get(sc);

			double logLike=fFull*Math.log(ftreeCorrectFull/ftreeIncorrectFull)-Math.log(ftreeCorrectRP/ftreeIncorrectRP)*fRef;
			//double logLike=(ftreeCorrect-ftreeIncorrect)*(fFull-fRef);
			System.err.println(sc+" log10("+ftreeCorrectFull+"/"+ftreeIncorrectFull+")*"+fFull+"-"+fRef+"*log10("+ftreeCorrectRP+"/"+ftreeIncorrectRP+")");
			loglikes.add(logLike);
		}
		return loglikes;
	}

	
	public static ArrayList<LogLikeDiff> getLogLikelihoods(HashMap<String,Double> freqFull,HashMap<String,Double> freqRef,HashMap<String,Double> freqTreeCorrectRP,HashMap<String,Double> freqTreeIncorrectRP,HashMap<String,Double> freqTreeCorrectFull,HashMap<String,Double> freqTreeIncorrectFull){
		ArrayList<LogLikeDiff> loglikes=new ArrayList<LogLikeDiff>();
		//double sum=0;
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			String sc=SiteClasses.siteClasses.get(i);
			double fFull=freqFull.get(sc);
			double fRef=freqRef.get(sc);
			double ftreeCorrectRP=freqTreeCorrectRP.get(sc);
			double ftreeIncorrectRP=freqTreeIncorrectRP.get(sc);
			double ftreeCorrectFull=freqTreeCorrectFull.get(sc);
			double ftreeIncorrectFull=freqTreeIncorrectFull.get(sc);
			//sum+=ftreeIncorrectRP;
			//if(i==0)System.err.println(sc+" log("+ftreeCorrectFull+"/"+ftreeIncorrectFull+")*"+fFull+"-"+fRef+"*log("+ftreeCorrectRP+"/"+ftreeIncorrectRP+")");
			double fullCorr=fFull*Math.log(ftreeCorrectFull);
			double fullInc=fFull*Math.log(ftreeIncorrectFull);
			double rpCorr=fRef*Math.log(ftreeCorrectRP);
			double rpInc=fRef*Math.log(ftreeIncorrectRP);
			double logLike=fullCorr-fullInc-(rpCorr-rpInc);

			loglikes.add(new LogLikeDiff(logLike, fullCorr, fullInc, rpCorr, rpInc));
		}//System.out.println(sum);
		
		return loglikes;
	}
	
	public static ArrayList<Double> getLogLikelihoods(HashMap<String,Double> freqFull,HashMap<String,Double> freqRef,HashMap<String,Double> freqTreeCorrect,HashMap<String,Double> freqTreeIncorrect,int count){
		ArrayList<Double> loglikes=new ArrayList<Double>();
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			String sc=SiteClasses.siteClasses.get(i);
			double fFull=freqFull.get(sc)/count;
			double fRef=freqRef.get(sc)/count;
			double ftreeCorrect=freqTreeCorrect.get(sc);
			double ftreeIncorrect=freqTreeIncorrect.get(sc);
			double logLike=Math.log10(ftreeCorrect/ftreeIncorrect)*(fFull-fRef);
			loglikes.add(logLike);
		}
		return loglikes;
	}
	
	
	public static void printSiteClassFreqs(HashMap<String,Double> s1,HashMap<String,Double> s2,File out,int count){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<SiteClasses.siteClasses.size();i++){
				String sc=SiteClasses.siteClasses.get(i);
				double f1=s1.get(sc)/count;
				double f2=s2.get(sc)/count;
				double diff=f1-f2;
				double logdiff=Math.log10(f1)-Math.log10(f2);
				bw.write(sc+"\t"+f1+"\t"+f2+"\t"+diff+"\t"+logdiff+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static HashMap<String,Double> addSiteClassFreqs(HashMap<String,Double> sc1,HashMap<String,Double> sc2){
		HashMap<String,Double> hm=new HashMap<String, Double>();
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			String sc=SiteClasses.siteClasses.get(i);
			hm.put(sc, sc1.get(sc)+sc2.get(sc));
		}
		return hm;
	}
	
	public static void printListComparison(ArrayList<Double> l1,ArrayList<Double> l2,HashMap<String,Double> freq1,HashMap<String,Double> freq2,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<l1.size();i++){
				bw.write(SiteClasses.siteClasses.get(i)+"\t"+l1.get(i)+"\t"+l2.get(i)+"\t"+(Math.log10(l1.get(i))-Math.log10(l2.get(i)))+"\t"+freq1.get(SiteClasses.siteClasses.get(i))+"\t"+freq2.get(SiteClasses.siteClasses.get(i))+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}

	
	HashMap<String,Double> siteClassFreqs=new HashMap<String, Double>();
	int numSites;
	public AlignmentSiteCount divide(AlignmentSiteCount div){
		AlignmentSiteCount alsc=new AlignmentSiteCount();
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			String sc=SiteClasses.siteClasses.get(i);
			alsc.siteClassFreqs.put(sc,this.siteClassFreqs.get(sc)/div.siteClassFreqs.get(sc));
		}
		return alsc;
		
	}
	

	
	public AlignmentSiteCount(ArrayList<Fasta> alignment){
		initSiteClassFreqs();
		calculateSiteClassFrequencies(alignment);
	}
	
	public AlignmentSiteCount(){
		initSiteClassFreqs();
	}
	public HashMap<String,Double> getSiteClassFreqs(){
		return siteClassFreqs;
	}
	
	public ArrayList<Double> getSiteClassFreqsList(){
		ArrayList<Double> list=new ArrayList<Double>();
		for(int i=0;i<SiteClasses.size;i++){
			list.add(siteClassFreqs.get(SiteClasses.siteClasses.get(i)));
		}
		return list;
	}
	public void printSiteClassFreqs(File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<SiteClasses.siteClasses.size();i++){
				bw.write(SiteClasses.siteClasses.get(i)+"\t"+siteClassFreqs.get(SiteClasses.siteClasses.get(i))+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	

	
	
	
	
	private void calculateSiteClassFrequencies(ArrayList<Fasta> fas){
		ArrayList<String> seqs=Fasta.getSequences(fas);
		int length=seqs.get(0).length();
		numSites=length;
		int seqNum=fas.size();
		for(int i=0;i<length;i++){
			StringBuffer col=new StringBuffer();
			for(int j=0;j<seqNum;j++){
				col.append(seqs.get(j).charAt(i));
			}
			String trans=SiteClasses.translate(col.toString());
			siteClassFreqs.put(trans, (1.0/numSites)+siteClassFreqs.get(trans));
		
		}

	}
	
	private void initSiteClassFreqs(){
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			siteClassFreqs.put(SiteClasses.siteClasses.get(i), 0.0);
		}
	}
	static HashMap<String,Double> init(){
		HashMap<String,Double> hm=new HashMap<String, Double>();
		for(int i=0;i<SiteClasses.siteClasses.size();i++){
			hm.put(SiteClasses.siteClasses.get(i), 0.0);
		}
		return hm;
	}
	

	
}
