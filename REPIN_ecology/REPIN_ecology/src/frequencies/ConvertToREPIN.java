package frequencies;

import java.io.*;
import java.util.*;

import REPINpopulations.REPINposition;
import util.*;



public class ConvertToREPIN {
	public class OrientationSequence{
		boolean orientation;
		String sequence;
		public OrientationSequence(boolean orientation,String seq){
			this.orientation=orientation;
			this.sequence=seq;
		}
	}
	ArrayList<Fasta> fas;
	ArrayList<Fasta> seeds;
	ArrayList<Fasta> REPINSeeds;
	int maxDist;
	String masterSeq;
	ArrayList</*list of positions in fasta seqs, each element in list corresponds to sequence in fasta file*/TreeMap<Integer/*pos*/,OrientationSequence/*orientation and sequence*/>> positions;
	HashMap<String,ArrayList<REPINpopulations.REPINposition>> REPINposition=new  HashMap<String,ArrayList<REPINpopulations.REPINposition>>();
	
	public static void main(String args[]) {
		File folder=new File("/Users/bertels/Documents/Arne/mutationFrequencies/chlororaphis/REPINs21/");
		String id="chl50083";
		for(int i=0;i<3;i++) {
			File seedSequences=new File(folder+"/"+id+"_"+i+"/"+id+"_"+i+".ss.REP");
			File genome=new File(folder+"/"+id+".fas");
			int maxDist=130;
			String masterSeq=getMasterSequence(folder, i, id);
			System.out.println(DNAmanipulations.reverse(masterSeq));
			ConvertToREPIN ctr=new ConvertToREPIN(seedSequences, genome, maxDist, masterSeq);
			ctr.write(new File(folder+"/"+id+"_"+i+".seqpos"));
		}
	}
	
	public static String getMasterSequence(File folder,int i,String id) {
		File in=new File(folder+"/"+id+"_"+i+"/"+id+"_"+i+".mw");
		String masterSeq="";
		try {
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line=br.readLine();
			masterSeq=line.split("\\s+")[0];
			br.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return masterSeq;
	}
	
	public ConvertToREPIN(File seedSequences,File genome,int maxDist,String masterSeq){
		fas=Fasta.readFasta(genome);
		this.masterSeq=masterSeq;
		System.out.println("Read REP seeds...");
		seeds=Fasta.readFasta(seedSequences);
		
		this.maxDist=maxDist;
		System.out.println("Calculate REP Positions...");
		calculatePositionsOrientation();
		System.out.println("Calculate REPINs...");
		REPINSeeds=getREPINs();
	} 
	
	public HashMap<String,ArrayList<REPINpopulations.REPINposition>> getPositions(){
		return REPINposition;
		
	}
	
	private ArrayList<TreeMap<Integer,OrientationSequence>> getClusters(TreeMap<Integer,OrientationSequence> repinPos){
		ArrayList<TreeMap<Integer,OrientationSequence>> clusters=new ArrayList<TreeMap<Integer,OrientationSequence>>();
		Integer[] pos=repinPos.keySet().toArray(new Integer[0]);
		TreeMap<Integer,OrientationSequence> currentCluster=new TreeMap<Integer, ConvertToREPIN.OrientationSequence>();
		for(int i=0;i<pos.length;i++){
			currentCluster.put(pos[i], repinPos.get(pos[i]));
			if(i<pos.length-1 && pos[i+1]-pos[i]<maxDist){
			}else{
				clusters.add(currentCluster);
				currentCluster=new TreeMap<Integer, OrientationSequence>();
			}

		}	
		return clusters;
	}
	
	private HashMap<Integer/*fasta id*/,ArrayList<TreeMap<Integer,OrientationSequence>>> getClusters(){
		HashMap<Integer,ArrayList<TreeMap<Integer,OrientationSequence>>> clusters=new HashMap<Integer,ArrayList<TreeMap<Integer,OrientationSequence>>>();
		for(int i=0;i<positions.size();i++){
			clusters.put(i,getClusters(positions.get(i)));
		}
		return clusters;
	}
	
	private ArrayList<Fasta> getREPINs(){
		HashMap<Integer/*sequence position in fasta file*/,ArrayList<TreeMap<Integer,OrientationSequence>>> clusters=getClusters();
	
		HashMap<String,ArrayList<String>> REPINfreq= makeREPINs(clusters);

		return toFasta(REPINfreq);
	}
	
	private ArrayList<Fasta> toFasta(HashMap<String,ArrayList<String>> repinFreq){
		ArrayList<Fasta> repins=new ArrayList<Fasta>();
		String[] keys=repinFreq.keySet().toArray(new String[0]);
		for(int i=0;i<keys.length;i++){
			int freq=repinFreq.get(keys[i]).size();
			StringBuffer positions=new StringBuffer(freq+"");
			for(int j=0;j<repinFreq.get(keys[i]).size();j++) {
				positions.append("\t"+repinFreq.get(keys[i]).get(j));
			}
			repins.add(new Fasta("wordOccurs"+positions.toString(),keys[i]));
			

		}
		return repins;
	}
	
	private HashMap<String,ArrayList<String>> makeREPINs(HashMap<Integer,ArrayList<TreeMap<Integer,OrientationSequence>>> clusters){
		HashMap<String,ArrayList<String>> repinFreq=new HashMap<String,ArrayList<String>>();
		Integer[] ids=clusters.keySet().toArray(new Integer[0]);
		for(int i=0;i<ids.length;i++){
			ArrayList<TreeMap<Integer,OrientationSequence>> positions=clusters.get(ids[i]);
			for(int j=0;j<positions.size();j++) {
				add(repinFreq,calculateREPINsFromClusters(positions.get(j), ids[i]));
			}
		}
		return repinFreq;
	}
	
	private void add(HashMap<String,ArrayList<String>> repinFreq,HashMap<String,ArrayList<String>> repinHM){
		String[] repins=repinHM.keySet().toArray(new String[0]);
		for(int i=0;i<repins.length;i++){
			if(!repinFreq.containsKey(repins[i])){
				repinFreq.put(repins[i],repinHM.get(repins[i]));
			}else{
				repinFreq.get(repins[i]).addAll( repinHM.get(repins[i]));
			}
		}
	}
	
	
	private HashMap<String,ArrayList<String>> calculateREPINsFromClusters(TreeMap<Integer,OrientationSequence> clusters,int id){
		HashMap<String,ArrayList<String>> REPINs=new HashMap<String,ArrayList<String>>();
		Integer[] pos=clusters.keySet().toArray(new Integer[0]);
		for(int i=0;i<pos.length;i++){
			String repin;
			int pos1;
			int pos2;
			if(i<pos.length-1){
				if(clusters.get(pos[i]).orientation!=clusters.get(pos[i+1]).orientation){
					
					//always construct REPINs so the positive orientation comes first and then the negative
					
					String s1=clusters.get(pos[i]).sequence;
					String s2=clusters.get(pos[i+1]).sequence;
					 repin=clusters.get(pos[i]).orientation==true?s1+s2:s2+s1;
					

					pos1=pos[i];
					pos2=pos[i+1]+clusters.get(pos[i+1]).sequence.length();

					i++;
				}else{
					repin=clusters.get(pos[i]).sequence+repeat("A",clusters.get(pos[i]).sequence.length());
					pos1=pos[i];
					pos2=pos1+clusters.get(pos[i]).sequence.length();

				}
			}else{
				repin=clusters.get(pos[i]).sequence+repeat("A",clusters.get(pos[i]).sequence.length());
				pos1=pos[i];
				pos2=pos1+clusters.get(pos[i]).sequence.length();
			}
			if(!REPINs.containsKey(repin)) {
				REPINs.put(repin, new ArrayList<String>());
			}
			if(!REPINposition.containsKey(repin)) {
				REPINposition.put(repin, new ArrayList<REPINposition>());
			}
			REPINs.get(repin).add(id+"_"+pos1+"_"+pos2);
			REPINposition.get(repin).add(new REPINposition(pos1,pos2, id));
		}
		return REPINs;
	}
	
	private String repeat(String s,int rep){
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<rep;i++){
			sb.append(s);
		}
		return sb.toString();
	}
	
	private  void initPositions(){
		positions=new ArrayList<TreeMap<Integer,OrientationSequence>>();

		for(int i=0;i<fas.size();i++) {
			positions.add(new TreeMap<Integer, ConvertToREPIN.OrientationSequence>());
		}
	}
	
	private void calculatePositionsOrientation(){
		initPositions();
		
		addAll(positions,find(true));
		addAll(positions,find(false));
	}
	
	private void addAll(ArrayList<TreeMap<Integer/*pos*/,OrientationSequence/*orientation and sequence*/>> posOrient,ArrayList<TreeMap<Integer/*pos*/,OrientationSequence/*orientation and sequence*/>> add){
		for(int i=0;i<posOrient.size();i++) {
			posOrient.get(i).putAll(add.get(i));
		}
	}
	
	private ArrayList<TreeMap<Integer,OrientationSequence>> find(boolean orientation){
		ArrayList<TreeMap<Integer,OrientationSequence>> posOrient=new ArrayList< TreeMap<Integer,OrientationSequence>>();

		for(int i=0;i<fas.size();i++){
			TreeMap<Integer,OrientationSequence> tm=new TreeMap<Integer, ConvertToREPIN.OrientationSequence>();
			for(int j=0;j<seeds.size();j++){
				String gseq=fas.get(i).getSequence();
				String sseq=orientation?seeds.get(j).getSequence():DNAmanipulations.reverse(seeds.get(j).getSequence());
				int k=0;
				while(k>-1){
					
					k=gseq.indexOf(sseq, k);

					if(k>-1){
						tm.put(k, new OrientationSequence(orientation,seeds.get(j).getSequence()));
						k++;
					}
					
				}
			}
			posOrient.add(tm);
		}

		return posOrient;
	}
	//turn the arraylist of treelists positions into a linear info file that contains the start,end and name of each REP sequence
	//it also contains the sequence, which we do not need here
	private ArrayList<Info> getInfo(ArrayList<Boolean> orientation){
		ArrayList<Info> REPpositions=new ArrayList<Info>();
		for(int i=0;i<positions.size();i++) {
			TreeMap<Integer,OrientationSequence> tm=positions.get(i);
			ArrayList<Integer> pos=new ArrayList<>(tm.descendingKeySet());
			for(int j=0;j<pos.size();j++) {
				OrientationSequence os=tm.get(pos.get(j));
				String name=os.sequence;
				int  start=pos.get(j);
				int  end=pos.get(j)+name.length();
				orientation.add(os.orientation);
				REPpositions.add(new Info(start,end,name));
			}
		}
		return REPpositions;
	}
	
	public void writeREPAnnotation(File outprefix,String seqid) {
		ArrayList<Boolean> orientation=new ArrayList<Boolean>();
		ArrayList<Info> infoREP=getInfo(orientation);
		WriteGenomeAnnotation.writeTab(infoREP, new File(outprefix+".tab"));
		WriteGenomeAnnotation.writeGff(infoREP,new File(outprefix+".gff3"),seqid,orientation,true,"REP");
	}
	
	public void write(File out){
		Fasta.write(REPINSeeds, out);
		String id=out.getName().split("\\.")[0];
		File outFolder=new File(out.getParent()+"/"+id+"/");
		outFolder.mkdir();
		writeWholeREPINs(outFolder,id);
	}
	
	private void writeWholeREPINs(File outFolder,String name) {
		for(int i=0;i<REPINSeeds.size();i++) {
			String[] pos=REPINSeeds.get(i).getIdent().split("\\s+");
			File out=new File(outFolder+"/"+name+"_"+REPINSeeds.get(i).getSequence()+"_"+(pos.length-1)+".fas");
			ArrayList<Fasta> fas=new ArrayList<Fasta>();
			for(int j=1;j<pos.length;j++) {
				fas.add(new Fasta(pos[j],getSequence(pos[j])));
			}
			Fasta.write(fas, out);
		}
	}
	
	private String getSequence(String pos) {
		int id=Integer.parseInt(pos.split("_")[0]);
		int start=Integer.parseInt(pos.split("_")[1]);
		int end=Integer.parseInt(pos.split("_")[2]);
		return fas.get(id).getSequence().substring(start,end);
	}
	
}
