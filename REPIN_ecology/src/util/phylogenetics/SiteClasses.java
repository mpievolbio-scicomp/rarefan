package util.phylogenetics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import util.Fasta;
import util.List_Array;





public class SiteClasses {
	static int A=0;
	static int T=1;
	static int C=2;
	static int G=3;
	//15 site classes
	public static ArrayList<String> siteClasses=new ArrayList<String>();
	static {
		siteClasses.add("AAAA");
		siteClasses.add("ATTT");
		siteClasses.add("ATAA");
		siteClasses.add("AATA");
		siteClasses.add("AAAT");
		siteClasses.add("AATT");
		siteClasses.add("ATTA");
		siteClasses.add("ATAT");
		siteClasses.add("AATC");
		siteClasses.add("ATTC");
		siteClasses.add("ATAC");
		siteClasses.add("ATCC");
		siteClasses.add("ATCA");
		siteClasses.add("ATCT");
		siteClasses.add("ATCG");
	}
	
	public static HashMap<String,Integer> siteClassNumbers=new HashMap<String, Integer>();
	static {
		for(int i=0;i<siteClasses.size();i++){
			siteClassNumbers.put(siteClasses.get(i), i);
		}
	}
	
	public static int size=15;
	public static void main(String args[]){
		File tree=new File(args[0]);
		Tree phylo=Phylogeny.readTree(tree);
		String[] order=new String[]{"S11","S12","S21","S22"};
		HashMap<String,Double> calculatedSiteClassFreqs=getLogLikelihoods(phylo, order);
		HashMap<String,Double> simulated=AlignmentSiteCount.getSiteClassFreqs(tree, order);
		double sCalc=0;
		double sSim=0;
		for(int i=0;i<siteClasses.size();i++){
			double calc=calculatedSiteClassFreqs.get(siteClasses.get(i));
			double sim=simulated.get(siteClasses.get(i));
			sCalc+=calc;
			sSim+=sim;
			System.out.println(calc+" "+sim);
		}
		System.out.println(sCalc+" "+sSim);
	}
	
	public static ArrayList<Integer> getSiteClassNumbers(String[] sites){
		ArrayList<Integer> list=new ArrayList<Integer>();
		for(int i=0;i<sites.length;i++){
			list.add(siteClassNumbers.get(sites[i]));
		}
		return list;
	}
	
	public static HashMap<String,Double> getLogLikelihoods(Tree tree,String[] idents){
		HashMap<String,Alignment> sitePatterns=translateFromSiteClasses(idents);
		HashMap<String,Double> loglike=new HashMap<String,Double>();
		
		for(int i=0;i<siteClasses.size();i++){
			AlignmentLogLikelihoods allh=new AlignmentLogLikelihoods(tree, sitePatterns.get(siteClasses.get(i)));
			loglike.put(siteClasses.get(i),allh.getTreeLogLike());
		}
		return loglike;
	}

	public static ArrayList<Double> readLikelihoods(File in, int lineNumber){
		ArrayList<Double> loglikes=new ArrayList<Double>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line;
			int i=-1;
			while((line=br.readLine())!=null){
				i++;
				if(i!=lineNumber){
					continue;
				}
				String split[]=line.split("\\s+");
				for(int j=1;j<split.length;j++){
					loglikes.add(Math.exp(Double.parseDouble(split[j])));
				}
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return loglikes;
	}
	
	public static HashMap<String,Double> getLogLikelihoodsRAxML(File tree,File alignment,String[] idents,File RAxMLPath,String suffix,boolean force){
		HashMap<String,Double> loglikePattern=new HashMap<String,Double>();
		//File sitePatternPhyl=new File(tree.getParent()+"/sitePatternRAxML.txt");
		//Alignment alg=makeSitePatternPhy(sitePatternPhyl,idents);
		Alignment alg=new Alignment(Fasta.readPhylip(alignment,idents));
//		System.err.println(alignment);
//		for(int i=0;i<alg.numSeqs;i++){
//			System.err.println(alg.getIdent(i));
//		}
		File likelihoods=RunTreePrograms.runRAxMLSiteLikelihood(tree, alignment, RAxMLPath, "sitell"+suffix,"", 1234,force);
		ArrayList<Double> raxmlLoglike=readLikelihoods(likelihoods,1);
		
		for(int i=0;i<alg.getLength();i++){
			String sc=alg.getColumn(i);
			double siteLoglike=raxmlLoglike.get(i);
			loglikePattern.put(sc,siteLoglike);

		}
		HashMap<String,Double> loglike=new HashMap<String,Double>();
		initLoglike(loglike,alg.length);
		Iterator<Entry<String,Double>> it=loglikePattern.entrySet().iterator();
		while(it.hasNext()){
			Entry<String,Double> e=it.next();
			String sc=translate(e.getKey());
			loglike.put(sc,loglike.get(sc) +e.getValue());
		}
		return loglike;
	}
	public static HashMap<String,Double> getLogLikelihoodsML(Tree tree,Alignment alg){
		AlignmentLogLikelihoods alll=new AlignmentLogLikelihoods(tree,alg);
		ArrayList<Double> mlLoglike=alll.getSitePatternLikelihood();
		ArrayList<String> sitePatterns=alll.getSitePatterns();
		HashMap<String,Double> loglike=new HashMap<String,Double>();
		initLoglike(loglike,alg.length);
		for(int i=0;i<mlLoglike.size();i++){
			String sc=translate(sitePatterns.get(i));
			loglike.put(sc,loglike.get(sc) +mlLoglike.get(i));
		}
		return loglike;
	}
	public static HashMap<String,Double> getLogLikelihoodsML(Tree tree,File alignment,String[] idents){
		Alignment alg=new Alignment(Fasta.readPhylip(alignment,idents));
		AlignmentLogLikelihoods alll=new AlignmentLogLikelihoods(tree,alg);
		ArrayList<Double> mlLoglike=alll.getSitePatternLikelihood();
		ArrayList<String> sitePatterns=alll.getSitePatterns();
		HashMap<String,Double> loglike=new HashMap<String,Double>();
		initLoglike(loglike,alg.length);
		for(int i=0;i<mlLoglike.size();i++){
			String sc=translate(sitePatterns.get(i));
			loglike.put(sc,loglike.get(sc) +mlLoglike.get(i));
		}
		return loglike;
	}

	
	public static void initLoglike(HashMap<String,Double> ll,int algLength){
		for(int i=0;i<siteClasses.size();i++){
			ll.put(siteClasses.get(i), 1.0/algLength);
		}
	}
	
	private static char intToChar(int i){
		if(i==A){
			return 'A';
		}
		if(i==T){
			return 'T';
		}
		if(i==C){
			return 'C';
		}
		if(i==G){
			return 'G';
		}
		return '-';
	}
	
	static HashMap<String,Alignment> initTrans(String[] order){
		HashMap<String,Alignment> hm=new HashMap<String, Alignment>();
		for(int i=0;i<siteClasses.size();i++){
			Alignment temp=new Alignment();
			ArrayList<String> idents=new ArrayList<String>();
			for(int j=0;j<order.length;j++){
				idents.add(order[j]);
			}
			
			temp.changeIdents(idents);
			hm.put(siteClasses.get(i),temp);
		}
		return hm;
	}
	
	public static Alignment makeSitePatternPhy(File out,String idents[]){
		Alignment alg=makeSitePatternAlignment(idents);
		Fasta.writePhylip(alg.toFasta(), out,100);
		return alg;
	}
	
	public static Alignment makeSitePatternAlignment(String idents[]){
		Alignment alg=new Alignment();
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					for(int l=0;l<4;l++){
						alg.addColumn(intToChar(i)+""+intToChar(j)+""+intToChar(k)+""+intToChar(l));
					}	
				}	
			}	
		}
		alg.changeIdents(List_Array.toString(idents));
		return alg;
	}
	
	public static HashMap<String,Alignment> translateFromSiteClasses(String[] idents){
		HashMap<String,Alignment> translation=initTrans(idents);
		Alignment alg=makeSitePatternAlignment(idents);
		for(int i=0;i<alg.length;i++){
			translation.get(translate(alg.getColumn(i))).addColumn(alg.getColumn(i));
		}
		return translation;
		
	}
	
	public static String translate(String col){
		//translate first letter to A next letter that is different to T next to C and last to G
		HashMap<Character,Character> chars=new HashMap<Character, Character>();
		String ATGC="ATCG";
		StringBuffer trans=new StringBuffer();
		for(int i=0;i<col.length();i++){
			char c=col.charAt(i);
			if(!chars.containsKey(c)){
				int size=chars.size();
				chars.put(c, ATGC.charAt(size));
			}
			trans.append(chars.get(c));
			
		}
		return trans.toString();
	}
}
