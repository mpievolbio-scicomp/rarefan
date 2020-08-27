package repinEvolution;

import java.io.*;
import java.util.*;

import REPINpopulations.REPINposition;
import frequencies.*;
import util.*;
//THIS HAS TO BE FIXED, SO IT CAN DEAL WITH MULTIPLE SEQUENCES PER FASTA FILE!
public class ClusterREPINs {
	int minREPINSize=30;
	REPINCluster rc;
	int proximity=130;
	String reference;
	String genome;
	String[] reps;
	public ClusterREPINs(REPINpopulations  refRepinPopulations,String reference,String refGenomeSeq,String[] reps) {
//		String[] ids=atg.getIDs();
		this.reps=reps;
		rc=new REPINCluster(reference, proximity);
		this.reference=reference;
		addReferenceREPINs(refRepinPopulations,refGenomeSeq);
	}
	

	
	private void addReferenceREPINs(REPINpopulations rp,String refGenomeSeq) {
		addREPINsToCluster(true, rp,null, reference,refGenomeSeq);

	}
	
	public void addREPINsToCluster(REPINpopulations rp,AlignTwoGenomes atg,String genomeID) {
		System.out.println(genomeID);
		addREPINsToCluster(false, rp,atg, genomeID,atg.getGenomeSequence(genomeID));
	}

	private void addREPINsToCluster(boolean isReference,REPINpopulations rp,AlignTwoGenomes atg,String genomeID,String genomeSeq) {
		
		
		for(int i=0;i<reps.length;i++) {
			REPINProperties props=rp.getREPINProperties(reps[i]);
			ArrayList<REPIN> repins=convertToRepin(toPosHashMap(props.getREPINPositions()),toPosHashMap(props.getLargestClusterREPINPositions()),i,genomeSeq,genomeID);
			for(int j=0;j<repins.size();j++) {
				int refPos=isReference?repins.get(j).start:atg.getPositionQuery(repins.get(j).start, genomeID);
				rc.addREPIN(reference, refPos, repins.get(j));
			}
		}
	}
	
	private ArrayList<REPINposition> toPosHashMap(HashMap<String,ArrayList<REPINposition>> repinPos){
		ArrayList<REPINposition> pos=new ArrayList<REPINposition>();
		String[] keys=repinPos.keySet().toArray(new String[0]);
		for(int i=0;i<keys.length;i++) {
			pos.addAll(repinPos.get(keys[i]));
		}
		return pos;
	}
	
	private ArrayList<REPIN> convertToRepin(ArrayList<REPINposition> pos,HashMap<Integer,Integer> largestCluster,int repintype,String genome,String genomeID){
		ArrayList<REPIN> repins=new ArrayList<REPIN>();
		
		for(int i=0;i<pos.size();i++) {
			int start=pos.get(i).start;
			int end=pos.get(i).end;
			int isLargestCluster;
			if(largestCluster.containsKey(start)) {
				isLargestCluster=1;
			}else {
				isLargestCluster=0;
			}
			if(start>end) {
				int mem=start;
				start=end;
				end=mem;
			}
			
			int size=end-start;
			int isREPIN;
			if(size<minREPINSize) {
				isREPIN=0;
			}else {
				isREPIN=1;
			}
			REPIN repin=new REPIN(start,end , repintype, new Fasta("repin_"+repintype+"_"+i,genome.subSequence(start, end)+""), genomeID,isREPIN,isLargestCluster);
			repins.add(repin);
		}
		return repins;
	}
	
	public HashMap<String,ArrayList<Fasta>> getREPINAlignments() {
		Integer[] repinPos=rc.getStart();

		HashMap<String,ArrayList<Fasta>> fasMap=new HashMap<String,ArrayList<Fasta>>();
		for(int i=0;i<repinPos.length;i++) {
			ArrayList<REPIN> repins=rc.getREPINCluster(repinPos[i]);
			ArrayList<Fasta> fas=toFasta(repins);
			if(fas.size()>0) {
				String name=fas.get(0).getSequence().length()+"_"+i;
				fasMap.put(name, fas);

			}
		}
		return fasMap;

	}
	
	public void writeREPINClusters(File out) {
		Integer[] repinPos=rc.getStart();
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write("cluster\t"+REPIN.getHeading());
			for(int i=0;i<repinPos.length;i++) {
				ArrayList<REPIN> repins=rc.getREPINCluster(repinPos[i]);
				for(int j=0;j<repins.size();j++) {
					bw.write(i+"\t"+repins.get(j).toString());
				}
			}
			bw.close();
		}catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	ArrayList<Fasta> toFasta(ArrayList<REPIN> repins){
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
		for(int i=0;i<repins.size();i++) {
			if(repins.get(i).isREPIN==1)fas.add(repins.get(i).toFasta());
		}
		return fas;
	}
	
	
}
