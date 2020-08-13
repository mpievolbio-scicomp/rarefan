package blastTools;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

import pairwiseAlignment.NeedlemanWunsch;

import util.Fasta;
import util.ReadGenbank;
import util.phylogenetics.Phylogeny;

public class Reciprocal_BestHits {
	static String blastall="/usr/local/bin/blastall";
	static String formatdb="/usr/local/bin/formatdb";
	public static void main(String args[]){
		File genome1=new File(args[0]);
		File genome2=new File(args[1]);
		ReadGenbank rgb=new ReadGenbank(genome1, "CDS");
		File CDS1=new File(genome1+"_CDS.fas");
		Fasta.write(rgb.getFeaturesSequences(),CDS1);
		rgb=new ReadGenbank(genome2, "CDS");
		File CDS2=new File(genome2+"_CDS.fas");
		Fasta.write(rgb.getFeaturesSequences(),CDS2);
		File out=new File(args[2]);
		File blast=new File(out.getParent()+"/blast1.out");
		PerformBlast.blast("blastall","formatdb", "blastn", 1e-5, blast, CDS1, CDS2, true,true, true, true);
		ReadBlast rb1=new ReadBlast(blast);
		HashMap<String,String> bestHits1=getBestHits(rb1);
		blast=new File(out.getParent()+"/blast2.out");
		PerformBlast.blast("blastall","formatdb", "blastn", 1e-5, blast, CDS2, CDS1, true,true, true, true);
		ReadBlast rb2=new ReadBlast(blast);
		HashMap<String,String> bestHits2=getBestHits(rb2);
		printReport(bestHits1,bestHits2,rb1,out,CDS1,CDS2);
	}
	
	public static HashMap<String,String> getBestHits(ReadBlast rb){
		HashMap<String,String> bestHits=new HashMap<String, String>();
		String prev="";
		for(int i=0;i<rb.getDatabase().size();i++){
			String query=rb.getQuery().get(i);
			if(prev.equals(query))continue;

			bestHits.put( query,rb.getDatabase().get(i));
			prev=query;
		}
		return bestHits;
	}
	
	public static HashMap<String,ArrayList<String>> reciprocalBestHits(File[] fas,File outFolder,boolean DNA,double covT){
		HashMap<String,ArrayList<String>> recBestHits=new HashMap<String, ArrayList<String>>();
		for(int i=0;i<fas.length;i++){
			for(int j=i+1;j<fas.length;j++){	
				putAll(recBestHits,reciprocalBestHits(fas[i], fas[j], outFolder, DNA,covT));
			}
		}
		return recBestHits;
	}
	
	public static HashMap<String,ArrayList<String>> reciprocalBestHitsReference(File[] fas,File outFolder,boolean DNA,double covT){
		HashMap<String,ArrayList<String>> recBestHits=new HashMap<String, ArrayList<String>>();
		for(int i=1;i<fas.length;i++){
				putAll(recBestHits,reciprocalBestHits(fas[0], fas[i], outFolder, DNA,covT));
		}
		return recBestHits;
	}
	
	public static HashMap<String,ArrayList<String>> reciprocalBestHitsNW(File[] fas,File outFolder,double covT,File matrix){
		HashMap<String,ArrayList<String>> recBestHits=new HashMap<String, ArrayList<String>>();
		for(int i=0;i<fas.length;i++){
			for(int j=i+1;j<fas.length;j++){	
				putAll(recBestHits,reciprocalBestHitsNW(fas[i], fas[j], outFolder, covT,matrix));
			}
		}
		//print(recBestHits);
		return recBestHits;
	}
	
	private static void putAll(HashMap<String,ArrayList<String>> bestHits,HashMap<String,ArrayList<String>> add){
		Iterator<Entry<String,ArrayList<String>>> it=add.entrySet().iterator();
		while(it.hasNext()){
			Entry<String, ArrayList<String>> e=it.next();
			if(!bestHits.containsKey(e.getKey())){
				bestHits.put(e.getKey(),e.getValue());
			}else{
				//not sure if this is ok
				bestHits.get(e.getKey()).addAll(e.getValue());
			}
		}	
	}
	
	private static void print(HashMap<String,ArrayList<String>> hm){
		Iterator<Entry<String,ArrayList<String>>> it=hm.entrySet().iterator();
		while(it.hasNext()){
			Entry<String, ArrayList<String>> e=it.next();
			System.out.println(e.getKey()+" "+e.getValue());
		}
	}
	
	private static void init(int[] index){
		for(int i=0;i<index.length;i++){
			index[i]=-1;
		}
	}
	
	public static HashMap<String,ArrayList<String>> reciprocalBestHitsNW(File fas1,File fas2,File outFolder,double covT,File matrix){
		HashMap<String,ArrayList<String>> recBestHits=new HashMap<String, ArrayList<String>>();
		ArrayList<Fasta> fasta1=Fasta.readFasta(fas1);
		ArrayList<Fasta> fasta2=Fasta.readFasta(fas2);
		double[] max1=new double[fasta1.size()];
		double[] max2=new double[fasta2.size()];
		int[] index1=new int[fasta1.size()];
		int[] index2=new int[fasta2.size()];
		init(index1);
		init(index2);
		for(int i=0;i<fasta1.size();i++){
			max1[i]=Double.NEGATIVE_INFINITY;
			for(int j=0;j<fasta2.size();j++){
				if(i==0){
					max2[j]=Double.NEGATIVE_INFINITY;
				}
				String seq1=fasta1.get(i).getSequence();
				String seq2=fasta2.get(j).getSequence();
				if(seq1.length()/seq2.length()>covT&&seq2.length()/seq1.length()>covT){
					NeedlemanWunsch nw=new NeedlemanWunsch(seq1, seq2,NeedlemanWunsch.readSimilarityMatrix(matrix),3.0,1.0);
					double score=nw.getScore();
					if(max1[i]<score){
						max1[i]=score;
						index1[i]=j;
					}
					if(max2[j]<score){
						max2[j]=score;
						index2[j]=i;
					}
				}
			}
		}
		for(int i=0;i<index1.length;i++){
			if(index1[i]!=-1&& i==index2[index1[i]]){
				String[] idents1=fasta1.get(i).getIdent().split("\\s+");
				String[] idents2=fasta2.get(index1[i]).getIdent().split("\\s+");
//				String ident1=idents1[1]+"_"+idents1[0];
//				String ident2=idents2[1]+"_"+idents2[0];
				String ident1=idents1[3];
				String ident2=idents2[3];
				if(!recBestHits.containsKey(ident1)){
					recBestHits.put(ident1, new ArrayList<String>());
				}
				if(!recBestHits.containsKey(ident2)){
					recBestHits.put(ident2, new ArrayList<String>());
				}
				recBestHits.get(ident1).add(ident2);
				recBestHits.get(ident2).add(ident1);
			}
		}
		
			return recBestHits;
	}
	
	private static String getOrg(File fas) {
		String split[]=fas.getName().split("_");
		StringBuffer org=new StringBuffer(split[0]);
		for(int i=1;i<split.length-1;i++) {
			org.append(split[i]);
		}
		return org.toString();
	}
	
	public static HashMap<String,ArrayList<String>> reciprocalBestHits(File fas1,File fas2,File outFolder,boolean DNA,double covT){
		HashMap<String,ArrayList<String>> recBestHits=new HashMap<String, ArrayList<String>>();
		File blast=new File(outFolder+"/blast.out");
		String org1=getOrg(fas1);
		String org2=getOrg(fas2);
		String program=DNA?"blastn":"blastp";
		PerformBlast.blast(blastall,formatdb, program, 1e-5, blast, fas1, fas2, true,true, DNA, true);
		ReadBlast rb1=new ReadBlast(blast);
		HashMap<String,String> bh1=getBestHits(rb1);
		PerformBlast.blast(blastall,formatdb, program, 1e-5, blast, fas2, fas1, true,true, DNA, true);
		ReadBlast rb2=new ReadBlast(blast);
		HashMap<String,String> bh2=getBestHits(rb2);
			ArrayList<Fasta> CDS1fas=Fasta.readFasta(fas1);
			ArrayList<Fasta> CDS2fas=Fasta.readFasta(fas2);
			HashMap<String,String> fasta=Fasta.fasToHash(CDS1fas,false);
			HashMap<String,Fasta> CDS1hm=Fasta.fasToFastaHash(CDS1fas, false);
			HashMap<String,Fasta> CDS2hm=Fasta.fasToFastaHash(CDS2fas, false);
			String prev="";
			for(int i=0;i<rb1.getDatabase().size();i++){
				String query=rb1.getQuery().get(i);
				String db=rb1.getDatabase().get(i);
				if(prev.equals(query))continue;
				if(bh2.containsKey(db)&&bh1.containsKey(query)){
					if(bh2.get(db).equals(query)&&bh1.get(query).equals(db)){
						double coverage=(1.0*rb1.getLength().get(i))/fasta.get(query).length();
						if(covT<coverage){
							String idents1[]=CDS1hm.get(query).getIdent().split("\\s+");
							String idents2[]=CDS2hm.get(db).getIdent().split("\\s+");
							String ident1=org1+"_"+idents1[2];
							String ident2=org2+"_"+idents2[2];
							if(!recBestHits.containsKey(ident1)){
								recBestHits.put(ident1, new ArrayList<String>());
							}
							if(!recBestHits.containsKey(ident2)){
								recBestHits.put(ident2, new ArrayList<String>());
							}
							recBestHits.get(ident1).add(ident2);
							recBestHits.get(ident2).add(ident1);

						}
					}
				}
				prev=query;
			}
			return recBestHits;
	}
	
//	private static void print(HashMap<String,String> test){
//		Iterator<Entry<String,String>> it=test.entrySet().iterator();
//		while(it.hasNext()){
//			Entry<String,String> e=it.next();
//			System.out.println(e.getKey()+" "+e.getValue());
//			
//		}
//	}
	
	public static void printReport(HashMap<String,String> bh1,HashMap<String,String> bh2,ReadBlast rb,File out,File CDS1, File CDS2){
		try{
			ArrayList<Fasta> CDS1fas=Fasta.readFasta(CDS1);
			ArrayList<Fasta> CDS2fas=Fasta.readFasta(CDS2);
			HashMap<String,String> fasta=Fasta.fasToHash(CDS1fas,false);
			HashMap<String,Fasta> CDS1hm=Fasta.fasToFastaHash(CDS1fas, false);
			HashMap<String,Fasta> CDS2hm=Fasta.fasToFastaHash(CDS2fas, false);
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write("Query\tDatabase\tidentities\tcoverage\n");
			String prev="";
			for(int i=0;i<rb.getDatabase().size();i++){
				String query=rb.getQuery().get(i);
				String db=rb.getDatabase().get(i);
				if(prev.equals(query))continue;
				if(bh2.containsKey(query)&&bh1.containsKey(db)){
					if(bh2.get(query).equals(db)&&bh1.get(db).equals(query)){
						double coverage=(1.0*rb.getLength().get(i))/fasta.get(query).length();
						bw.write(getLocusTag(CDS1hm.get(query).getIdent())+"\t"+getLocusTag(CDS2hm.get(db).getIdent())+"\t"+(rb.getScore().get(i)/2)+"\t"+coverage+"\n");
					}
				}
				prev=query;
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static String getLocusTag(String ident){
		String split[]=ident.split("\\s+");
		for(int i=0;i<split.length;i++){
			if(split[i].equals("locus_tag=")){
				return split[i+1];
			}
		}
		return "noLocusTag";
	}
	
	
}
