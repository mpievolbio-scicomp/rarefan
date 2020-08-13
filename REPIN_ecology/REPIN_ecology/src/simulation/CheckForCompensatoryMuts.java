package simulation;

import java.io.*;
import java.util.*;

import frequencies.RAYTfrequencyGraph;
import frequencies.REPINProperties;
import frequencies.util.SimilarityNetwork;
import util.*;

public class CheckForCompensatoryMuts {
	public static void main(String args[]){
		File inFolder=new File(args[0]);
		File out=new File(inFolder+"/compensatory.out");
		CheckForCompensatoryMuts cfcm=new CheckForCompensatoryMuts(inFolder);
		cfcm.writeResults(out);

	}
	int numMuts=2;
	ArrayList<Fasta> seqs=new ArrayList<Fasta>();
	HashMap<String,Integer> doubles=new HashMap<String,Integer>();
	HashMap<String,Integer> compensatory=new HashMap<String, Integer>();
	HashMap<String,Integer> noncompensatory=new HashMap<String, Integer>();
	String master;
	
	private void write(BufferedWriter bw,HashMap<String,Integer> comp,String name)throws IOException {
		String[] keys=comp.keySet().toArray(new String[0]);
		int number=0;
		for(int i=0;i<keys.length;i++){

			bw.write(name+"\t"+keys[i]+"\t"+comp.get(keys[i])+"\n");
		}
	}
	
	public void writeResults(File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write("name\tstrain\tseqs\n");
			write(bw,compensatory,"compensatory");
			write(bw,noncompensatory,"noncompensatory");
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public CheckForCompensatoryMuts(File inFolder){
		File[] files=inFolder.listFiles();
		for(int i=0;i<files.length;i++){
			if(files[i].getName().endsWith("largestCluster.ss")){
				File sequences=files[i];
				ArrayList<Fasta> fas=Fasta.readFasta(sequences);
				String master=REPINProperties.getMaxSequence(sequences);
				String ident=files[i].getName().split("\\.")[0];
				runIndividual(fas, master,ident);
			}
		}
	}
	
	public void runIndividual(ArrayList<Fasta> fas,String master,String ident){
		seqs=fas;
		this.master=master;
		
		doubles=getMutations();
		countCompensatory(ident);
	}
	
	private void countCompensatory(String ident){
			ArrayList<String> comps=new ArrayList<String>();
			String[] d=doubles.keySet().toArray(new String[0]);
			int comp=0;
			int count=0;
			for(int j=0;j<d.length;j++){


				ArrayList<Character> muts=getMutations(master,d[j]);
				int number=doubles.get(d[j]);

				if(isComplement(muts.get(0),muts.get(1))){
					comp+=number;
					comps.add(d[j]+""+number);

				}else{

				}
				count+=number;

			}
			compensatory.put(ident, comp);
			noncompensatory.put(ident, count-comp);


	}
	
	private char getComplement(char c){
		String correct="ATCG";
		String complement="TAGC";
		int pos=correct.indexOf(c);
		return complement.charAt(pos);
	}
	
	private boolean isComplement(char c1,char c2){
		return (c1==getComplement(c2));
	}
	
	private ArrayList<Character> getMutations(String s1,String s2){
		ArrayList<Character> list=new ArrayList<Character>();
		for(int i=0;i<s1.length();i++){
			char c1=s1.charAt(i);
			char c2=s2.charAt(i);
			if(c1!=c2){
				list.add(c2);
			}
		}
		return  list;
	}
	
	private HashMap<String,Integer> getMutations(){
		HashMap<String,Integer>  muts=new HashMap<String,Integer> ();
		for(int i=0;i<seqs.size();i++){
			String orig=seqs.get(i).getSequence();
			if(orig.equals(master))continue;
			int diff=SimilarityNetwork.getNumDifferences(orig,master);
			if(diff==2){
				muts.put(orig,RAYTfrequencyGraph.getFreq(seqs.get(i).getIdent()));
			}
		}
		return muts;
	}
	

	
}
