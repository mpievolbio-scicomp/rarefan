package repinEvolution;

import java.io.*;
import java.util.*;

public class AnalyseREPINClusters {
	
	//output: genomeID1(ref) genomeID2 REPINtype1 REPINType2 isREPIN1 isREPIN2 count

	REPINCluster rc;
	HashMap<REPINTransition,Integer> rts;
	String reference;
	HashSet<String> queries;
	public AnalyseREPINClusters(REPINCluster rc,String ref,HashSet<String> queries) {
		this.rc=rc;
		this.reference=ref;
		this.queries=queries;
		rts=new HashMap<REPINTransition,Integer>();
		calculateREPINTransition();
	}
	
	public void writeTransitionStats(File out) {
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			REPINTransition[] rt=rts.keySet().toArray(new REPINTransition[0]);
			for(int i=0;i<rt.length;i++) {
				bw.write(rt[i].toString()+"\t"+rts.get(rt[i])+"\n");
			}
			bw.close();
		}catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	//input is a single cluster
	//calculate transition: genomeID1(ref) genomeID2 REPINtype1 REPINType2 isREPIN1 isREPIN2
	private void calculateREPINTransition() {
		Integer[] pos=rc.getStart();
		for(int i=0;i<pos.length;i++) {
			ArrayList<REPIN> cluster=rc.getREPINCluster(pos[i]);
			REPIN refREPIN=getReference(cluster);
			ArrayList<REPIN> queryREPINs=getQueries(cluster);
			for(int j=0;j<queryREPINs.size();j++) {	
				addREPINTransition(new REPINTransition(new REPINSignature(refREPIN),new REPINSignature(queryREPINs.get(j)),reference));
			}
		}
		
	}
	
	private void addREPINTransition(REPINTransition rt) {
		if(!rts.containsKey(rt)) {
			rts.put(rt, 0);
		}
		rts.put(rt,rts.get(rt)+1);
		REPINTransition rts[]=this.rts.keySet().toArray(new REPINTransition[0]);

	}
	
	private ArrayList<REPIN> getQueries(ArrayList<REPIN> cluster){
		ArrayList<REPIN> queries=new ArrayList<REPIN>();
		HashSet<String> queriesHS=new HashSet<String>(this.queries);
		//System.out.println("clustersize "+cluster.size());
		for(int i=0;i<cluster.size();i++) {
			if(!cluster.get(i).genomeID.equals(reference)) {
				queries.add(cluster.get(i));
				queriesHS.remove(cluster.get(i).genomeID);
			}
		}
		String[] remaining=queriesHS.toArray(new String[0]);
		for(int i=0;i<remaining.length;i++) {
			queries.add(new REPIN(remaining[i]));
		}
		//System.out.println("queries size "+queries.size());
		return queries;
	}
	
	private REPIN getReference(ArrayList<REPIN> cluster) {
		for(int i=0;i<cluster.size();i++) {
			if(cluster.get(i).genomeID.equals(reference)) {
				return cluster.get(i);
			}
		}
		return new REPIN(reference);
	}
	


}
