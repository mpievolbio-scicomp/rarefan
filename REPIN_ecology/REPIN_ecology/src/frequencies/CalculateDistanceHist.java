package frequencies;

import java.io.*;
import java.util.*;

//import rcaller.*;
import util.*;

public class CalculateDistanceHist {
	static File scriptPath=new File("/usr/local/bin/Rscript");
	public static void main(String args[]){
		File inFolder=new File(args[0]);
		CalculateDistanceHist cdh=new CalculateDistanceHist();
		cdh.writeAll(inFolder,6);
		//cdh.plot();
	}
	HashMap<String,HashMap<Integer,Integer>> hists=new HashMap<String,HashMap<Integer,Integer>>();
	ArrayList<File> outFiles=new ArrayList<File>();
	File inFolder;
	public void writeAll(File inFolder,int maxMuts){
		this.inFolder=inFolder;
		ArrayList<File> files=getFiles(inFolder);
		for(int i=0;i<files.size();i++){
			calculateHist(files.get(i),maxMuts);
		}
		File out=new File(inFolder+"/distHists.txt");
		writeAllHistograms(out);
	}
	
	public void write(File in,File out){
		hists=new HashMap<String,HashMap<Integer,Integer>>();
		calculateHist(in,-1);
		writeAllHistograms(out);
	}
	public void writeMCH(File in,File out){
		hists=new HashMap<String,HashMap<Integer,Integer>>();
		calculateHistMCH(in,-1);
		writeAllHistograms(out);
	}
	private void writeAllHistograms(File out){
		try{
			outFiles.add(out);

			BufferedWriter bw=new BufferedWriter(new FileWriter(out));

			String[] names=hists.keySet().toArray(new String[0]);
			for(int i=0;i<names.length;i++){
				writeHistogram(hists.get(names[i]),bw,names[i]);
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
//	public void plot(){
//		File out=new File(inFolder+"/allHists.pdf");
//		RCode rc=R_functions.plot_InitPdf(out, 50, 50);
//		int num=outFiles.size();
//		rc.addRCode("par(mfrow=c(3,"+(num/3+1)+"))");
//		for(int i=0;i<num;i++){
//			rc.addRCode("t<-read.table(\""+outFiles.get(i)+"\",header=FALSE)");
//			rc.addRCode("barplot(t[,2],names.arg=t[,1],main=\""+outFiles.get(i).getName()+"\")");
//		}
//		R_functions.writeRCode(rc, new File(out+".R"));
//		R_functions.runRCode(rc, scriptPath);
//	}
	
	private void writeHistogram(HashMap<Integer,Integer> hist,BufferedWriter bw,String name) throws IOException{
			Integer[] keys=hist.keySet().toArray(new Integer[0]);
			for(int i=0;i<keys.length;i++){
				String names[]=name.split("\\.");
				String seql=names[1];
				String n=names[0];
				bw.write(n+"\t"+keys[i]+"\t"+hist.get(keys[i])+"\t"+seql+"\n");
			}

	}

	
	private void calculateHist(File in,int maxMuts){
		String split[]=in.getName().split("\\.");
		String name=split[split.length-2];
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		HashMap<Integer,Integer> histogram=getHistogram(fas,maxMuts);
		int seql=calculateSequenceLength(fas);
		hists.put(name+"."+seql, histogram);
	}
	
	private int calculateSequenceLength(ArrayList<Fasta> fas){
		String template=fas.get(0).getSequence();
		int l=template.length();
		String last=template.substring(l/2);
		for(int i=0;i<fas.size();i++){
			if(!fas.get(i).getSequence().substring(l/2).equals(last)){
				return l;
			}
		}
		return l/2;
	}
	
	private void calculateHistMCH(File in,int maxMuts){
		ArrayList<ArrayList<Fasta>> mutationClasses=getMutationClasses(in,maxMuts);
		for(int i=0;i<mutationClasses.size();i++){
			HashMap<Integer,Integer> histogram=getHistogramMCH(mutationClasses.get(i));
			hists.put(i+"", histogram);

		}
	}
	private ArrayList<ArrayList<Fasta>> getMutationClasses(File in,int maxMuts){
		ArrayList<ArrayList<Fasta>> mutclasses=new ArrayList<ArrayList<Fasta>>();
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		int maxIndex=getMaxIndex(fas);
		String masterSeq=fas.get(maxIndex).getSequence();

		for(int i=0;i<fas.size();i++){
			int numDiff=getNumDiff(fas.get(i).getSequence(),masterSeq);
			if(numDiff>maxMuts&&maxMuts!=-1){
				numDiff=maxMuts;
			}
			while(numDiff>=mutclasses.size()){
				mutclasses.add(new ArrayList<Fasta>());
			}
			mutclasses.get(numDiff).add(fas.get(i));
		}

		return mutclasses;
	}
	
	private HashMap<Integer,Integer> getHistogramMCH(ArrayList<Fasta> fas){
		HashMap<Integer,Integer> hist=new HashMap<Integer,Integer>();
		for(int i=0;i<fas.size();i++){
			int freq=RAYTfrequencyGraph.getFreq(fas.get(i).getIdent());
			if(!hist.containsKey(freq)){
				hist.put(freq, 0);
			}
			hist.put(freq,hist.get(freq)+1);
		}
		return hist;
	}
	
	private HashMap<Integer,Integer> getHistogram(ArrayList<Fasta> fas,int maxMuts){
		HashMap<Integer,Integer> hist=new HashMap<Integer,Integer>();
		int maxIndex=getMaxIndex(fas);
		String masterSeq=fas.get(maxIndex).getSequence();
		for(int i=0;i<fas.size();i++){
			int numDiff=getNumDiff(fas.get(i).getSequence(),masterSeq);
			if(numDiff>maxMuts&&maxMuts!=-1){
				numDiff=maxMuts;
			}
			int freq=RAYTfrequencyGraph.getFreq(fas.get(i).getIdent());
			if(!hist.containsKey(numDiff)){
				hist.put(numDiff, 0);
			}
			hist.put(numDiff,hist.get(numDiff)+freq);
		}
		return hist;
	}
	
	public static int getNumDiff(String a,String b){
		int diffs=0;
		for(int i=0;i<a.length();i++){
			if(a.charAt(i)!=b.charAt(i)){
				diffs++;
			}
		}
		return diffs;
	}
	
	
	private int getMaxIndex(ArrayList<Fasta> fas){
		int max=0;
		int maxIndex=-1;
		for(int i=0;i<fas.size();i++){
			int occ=RAYTfrequencyGraph.getFreq(fas.get(i).getIdent());
			if(max<occ){
				max=occ;
				maxIndex=i;
			}
		}
		return maxIndex;
	}
	
	private ArrayList<File> getFiles(File in){
		File[] list=in.listFiles();
		ArrayList<File> selection=new ArrayList<File>();
		for(int i=0;i<list.length;i++){
			String name=list[i].getName();
			if(name.endsWith("largestCluster.ss")){
				selection.add(list[i]);
			}
		}
		return selection;
	}
	
	
}
