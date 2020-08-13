package blastTools;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;




import util.*;

public class BuildBlastAlignment {
	ArrayList<Integer> alLength = new ArrayList<Integer>();
	ArrayList<String> alName = new ArrayList<String>();
	HashMap<String,int[]> coverage;
	HashMap<String,int[]> coveragePos;
	HashMap<String,int[]> coverageNeg;

	HashMap<String,String> fasta;
	HashMap<String,Double> avgCov=new HashMap<String, Double>();
	HashMap<String,Double> percentCov=new HashMap<String,Double>();
	HashMap<String,InfoTree> matchedRegions=new HashMap<String, InfoTree>();
	

	public HashMap<String,int[]> getCoverage(){
		return coverage;
	}
	public HashMap<String,int[]> getCoveragePos(){
		return coveragePos;
	}
	public HashMap<String,int[]> getCoverageNeg(){
		return coverageNeg;
	}
	public int[] getCoverage(String id){
		return coverage.get(id);
	}
	public int[] getCoveragePos(String id){
		return coveragePos.get(id);
	}
	public int[] getCoverageNeg(String id){
		return coverageNeg.get(id);
	}
	public String getGene(String id){
		return fasta.get(id);
	}

	
	
	public BuildBlastAlignment(File RefSeq,int flank,File blast,double eValue,double overlap,HashMap<String,Boolean> exclude,double perBaseScore,HashMap<String,Boolean> excludeGenes,HashMap<String,Boolean> showOverlap,HashMap<String,Integer> geneHash,ArrayList<ArrayList<String>> groups){
		readSeq(RefSeq, flank);
		fasta=readFasta(RefSeq);
		coverage = initCoverage(alName,alLength);
		coveragePos = initCoverage(alName,alLength);
		coverageNeg = initCoverage(alName,alLength);
		readBlast(blast, eValue,flank,overlap,exclude,perBaseScore,excludeGenes,showOverlap,geneHash,groups);

		setAvgCoverage();
		setPercentCoverage();
	}
	
	public boolean checkPercentCoverage(String id,double t){
		
		return percentCov.get(id)>t;
	}
	public double getPercentCoverage(String id){
		
		return percentCov.get(id);
	}
	public boolean checkAvgCoverage(String id,double t){
		return avgCov.get(id)>t;
	}
	
	public void setAvgCoverage(){
		for(int i=0;i<alName.size();i++){
			int sum=0;
			for(int j=0;j<alLength.get(i);j++){
				sum+=coverage.get(alName.get(i))[j];
			}
			avgCov.put(alName.get(i), (sum/(alLength.get(i)*1.0)));
		}
		
	}
	
	public void setPercentCoverage(){
		for(int i=0;i<alName.size();i++){
			int sum=0;
			for(int j=1;j<alLength.get(i);j++){
				if(coverage.get(alName.get(i))[j]>0){
					sum++;
				}
			}
			percentCov.put(alName.get(i), (sum/((alLength.get(i)-1.0))));
		}
		
	}
	
	public static HashMap<String,String> readFasta(File in){
		HashMap<String,String> hm=new HashMap<String, String>();
 		ArrayList<Fasta> fas=Fasta.readFasta(in);
		for(int i=0;i<fas.size();i++){
			String id=fas.get(i).getIdent().startsWith(">")?fas.get(i).getIdent().substring(1).split("\\s+")[0]:fas.get(i).getIdent().split("\\s+")[0];
			hm.put(id,fas.get(i).getSequence());
		}
		return hm;
	}
	

	
	public static HashMap<String,int[]> initCoverage(ArrayList<String> alName,ArrayList<Integer> alLength){
		HashMap<String,int[]> coverage=new HashMap<String, int[]>();
		for(int i=0;i<alName.size();i++){
			coverage.put(alName.get(i), new int[alLength.get(i)+1]);
		}
		return coverage;
	}
	

	public static String getStrain(String id){
		StringBuffer strain=new StringBuffer();
		String split[]=id.split("-");
		for(int i=0;i<split.length-2;i++){
			strain.append(split[i]+"-");
		}
		strain.append(split[split.length-2]);
		return strain.toString();
	}
	
	
	
	public void readBlast(File blast, double eValue,int flank,double overlap,HashMap<String,Boolean> exclude,double perBaseScore,HashMap<String,Boolean> excludeGenes,HashMap<String,Boolean> showOverlap,HashMap<String,Integer> geneHash,ArrayList<ArrayList<String>> groups) {
		ReadBlast rb=new ReadBlast(blast);
		
		for(int i=0;i<rb.getDatabase().size();i++){
			String id=rb.getDatabase().get(i);
			String strain=getStrain(id);
			String idQ=rb.getQuery().get(i);
			int startDB=rb.getStartDB().get(i);
			int endDB=rb.getEndDB().get(i);
			int startQ=rb.getStartQuery().get(i);
			int endQ=rb.getEndQuery().get(i);
			int start=startQ>endQ?endQ:startQ;
			int end=endQ>startQ?endQ:startQ;
			double score=rb.getScore().get(i);
			int length = Math.abs(startDB-endDB)+1;

			if(eValue<rb.getEvalue().get(i))continue;
			if(exclude.containsKey(strain))continue;
			if(perBaseScore>(score/length))continue;
			if(excludeGenes.containsKey(id))continue;
			
			

			boolean cont=false;
			int sover=(int)(start+length*overlap);
			int eover=(int)(end-length*overlap);
			if(sover<eover){//CAUTION! OVERLAPS OF UP TO 'overlap' BP TOLERATED
				Info info=new Info(sover,eover,score+"\t"+id);
				if(showOverlap.containsKey(id)){
					int group=geneHash.get(id);
					System.out.println(id+" "+group+" "+groups.get(group).size());
					showOverlap(info,idQ,geneHash,groups);
				}

				if(matchedRegions.containsKey(idQ)){
					InfoTree it=matchedRegions.get(idQ);

					if(it.checkOverlap(info)){
						ArrayList<Info> al=new ArrayList<Info>();
						it.search(info, al);
						for(int j=0;j<al.size();j++){
							if(Double.parseDouble(al.get(j).info.split("\t")[0])>1.1*score){//parametrize 1.1???
								cont=true;
								break;
							}
						}
					}else{
						matchedRegions.get(idQ).insert(info);
					}
				}else{
					InfoTree temp=new InfoTree();
					temp.insert(info);
					matchedRegions.put(idQ, temp);
				}
			}else{
				continue;
			}
			if(cont)continue;
			char orientation=startDB>endDB?'+':'-';
			int pos = startDB>endDB?endDB:startDB;
			if(pos<=flank){
				length=(pos+length)-flank;
				pos=0;
			}else{
				pos=pos-flank;
			}
			if(!coverage.containsKey(id)){
				System.err.println(id+" not found in database.");
				continue;
			}
			for (int j = pos; j < pos + length &&j<coverage.get(id).length; j++) {
				if (coverage.get(id).length > j)
					coverage.get(id)[j]++;
				if(orientation=='+'){
					coveragePos.get(id)[j]++;
				}else{
					coverageNeg.get(id)[j]++;
				}
			}
		}



	}

	public void showOverlap(Info inf,String id,HashMap<String,Integer> geneHash,ArrayList<ArrayList<String>> groups){
		if(matchedRegions.containsKey(id)){
			InfoTree it=matchedRegions.get(id);
			ArrayList<Info> al=new ArrayList<Info>();
			it.search(inf, al);
			for(int i=0;i<al.size();i++){
				Info temp=al.get(i);
				String split[]=temp.info.split("\t");
				String gene=split[1];
				String score=split[0];
				int group=geneHash.get(gene);
				int gs=groups.get(group).size();
				System.out.println(temp.getStart()+" "+temp.getEnd()+" "+score+" "+gene+" "+group+" "+gs);
			}
		}else{
			System.out.println("none");
		}
	}
	

	public void readSeq(File ref, int flank) {
		try {
			BufferedReader is = new BufferedReader(new FileReader(ref));
			String line = "";
			StringBuilder sb = new StringBuilder();
			String name = "";
			while ((line = is.readLine()) != null) {
				if (!line.startsWith(">"))
					sb.append(line.trim());
				else {
					if (name.length() > 0) {
						alName.add(name.split("\\s+")[0]);
						alLength.add(sb.length()-2*flank);
						sb = new StringBuilder();
					}
					name = line.substring(1);

				}
			}
			alName.add(name.split("\\s+")[0]);
			alLength.add(sb.length());
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}


}
