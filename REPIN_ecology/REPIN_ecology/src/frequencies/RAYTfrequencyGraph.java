package frequencies;

import java.io.*;
import java.util.*;

import util.*;

public class RAYTfrequencyGraph {
	public static void main(String args[]){
		File alignment=new File(args[0]);
		String group[]=alignment.getName().split("\\.");
		File out=new File(alignment.getParentFile()+"/mutfreq_"+group[0]+".txt");
		RAYTfrequencyGraph rfg=new RAYTfrequencyGraph(alignment);
		rfg.writeMutationFrequencies(out);
	}
	
	public static HashMap<String,Double> read(File in){
		HashMap<String,Double> mutFreq=new HashMap<String,Double>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String split[]=line.split("\\s+");
				mutFreq.put(split[0], Double.parseDouble(split[1]));
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return mutFreq;
	}
	
	HashMap<String,Double> mutFreq=new HashMap<String,Double>();
	ArrayList<Fasta> alignment;
	
	public RAYTfrequencyGraph(File alg){
		alignment=Fasta.readFasta(alg);
		calculateMutationFrequencies();
	}
	
	public HashMap<String,Double> getMutationFrequencies(){
		return mutFreq;
	}
	
	private void calculateMutationFrequencies(){
		String anc=getAnc(alignment);
		mutFreq=determineMutationFreqs(alignment,anc);

	}
	
	public void writeMutationFrequencies(File out){
		printMutFreq(mutFreq,out);
	}
	
	public static void printMutFreq(HashMap<String,Double> mf,File out){
		try{
			if(mf.size()>0){
				BufferedWriter bw=new BufferedWriter(new FileWriter(out));
				String[] keys=mf.keySet().toArray(new String[0]);
				for(int i=0;i<keys.length;i++){
					bw.write(keys[i]+"\t"+mf.get(keys[i])+"\n");
				}
				bw.close();
			}
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private static HashMap<String,Double> determineMutationFreqs(ArrayList<Fasta> fas,String anc){
		HashMap<String,Double> mutfreq=new HashMap<String,Double>();
		int popsize=getPopSize(fas);
		//System.out.println(popsize);
		for(int j=0;j<anc.length();j++){
			char a=anc.charAt(j);
			for(int i=0;i<fas.size();i++){
				double freq=getFreq(fas.get(i).getIdent());
				String seq=fas.get(i).getSequence();

				char s=seq.charAt(j);
				if(a!=s){
					String mut=a+""+j+""+s;
					if(!mutfreq.containsKey(mut)){
						mutfreq.put(mut, 0.0);
					}
					mutfreq.put(mut,mutfreq.get(mut)+(freq/popsize));
				}
			}
			
		}
		return mutfreq;
	}
	
	private static int getPopSize(ArrayList<Fasta> fas){
		int popsize=0;
		for(int i=0;i<fas.size();i++){
			int freq=getFreq(fas.get(i).getIdent());
			popsize+=freq;
		}
		return popsize;
	}

	
	private static String getAnc(ArrayList<Fasta> fas){
		String anc=fas.get(0).getSequence();
		int max=getFreq(fas.get(0).getIdent());
		for(int i=0;i<fas.size();i++){
			int freq=getFreq(fas.get(i).getIdent());
			if(max<freq){
				anc=fas.get(i).getSequence();
				max=freq;
			}
		}
		return anc;
	}
	
	public static int getFreq(String ident){
		String ns= ident.replaceAll("[A-Za-z]", "");
		return Integer.parseInt(ns.split("\\s+")[0]);
	}
	
}
