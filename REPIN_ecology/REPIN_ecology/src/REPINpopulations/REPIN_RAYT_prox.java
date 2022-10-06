package REPINpopulations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import frequencies.REPINProperties;
import identifyRAYTs.BlastRAYTs;
import util.DNAmanipulations;
import util.Fasta;
import util.Info;

public class REPIN_RAYT_prox {
	int vicinityDistance;
	ArrayList<ProxStats> stats=new ArrayList<ProxStats>();
	File outFolder;
	ArrayList<Fasta> allRAYTs=new ArrayList<Fasta>();
	ArrayList<Info> allRAYTsPos=new ArrayList<Info>();
	ArrayList<Vicinity> allVicinity=new ArrayList<Vicinity>();
	int repinGroups;
	public class Vicinity{
		ArrayList<Integer> REPINgroups=new ArrayList<Integer>();
		public void add(int group) {
			REPINgroups.add(group);
		}
		public String toString() {
			StringBuffer sb=new StringBuffer();
			if(REPINgroups.size()>0) {
				sb.append(REPINgroups.get(0));
			}
			for(int i=1;i<REPINgroups.size();i++) {
				sb.append(","+REPINgroups.get(i));
			}
			return sb.toString();
		}
		public ArrayList<Integer> get(){
			return REPINgroups;
		}
	}

	public REPIN_RAYT_prox(File outFolder,int repinGroups,int distanceRAYTGene) {
		this.vicinityDistance=distanceRAYTGene;
		this.outFolder=outFolder;
		this.repinGroups=repinGroups;
	}
	private void addRAYTs(ArrayList<Fasta> rayts,ArrayList<Info> pos,ArrayList<Vicinity> vicinity) {
		allRAYTs.addAll(rayts);
		allRAYTsPos.addAll(pos);
		allVicinity.addAll(vicinity);
	}

	public void addRAYT(ArrayList<Info> raytPos,File genome,ArrayList<REPINGenomePositions> rgp) {
		ArrayList<Fasta> rayts=makeRAYTFasta(raytPos,genome);
		ArrayList<Vicinity> vicinity=getVicinityInformation(raytPos,rgp,RAREFAN_MAIN.getGenomeID(genome));
		addRAYTs(rayts,raytPos,vicinity);
	}
	
	public void write(File out) {
		Fasta.write(allRAYTs, new File(out+".fas"));
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write("Genome\tRAYT\tREPINgroups\n");
			for(int i=0;i<allRAYTs.size();i++) {
				String split[]=allRAYTs.get(i).getIdent().split("_"); 
				String raytNumber=split[split.length-1];
				String genome=getGenome(allRAYTs.get(i).getIdent());
				bw.write(genome+"\t"+raytNumber+"\t"+allVicinity.get(i).toString()+"\n");
			}
			bw.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private String getGenome(String ident) {
		String split[]=ident.split("_");
		String genome=split[0];
		for(int j=1;j<split.length-1;j++) {
		  genome=genome+"_"+split[j];
		}
		return genome;
	}
	
	private HashMap<String,String> addRAYTsToREPINTypeDS(HashMap<String,String> repintypeDS){
		for(int i=0;i<allRAYTs.size();i++) {
			ArrayList<Integer> repintypes=allVicinity.get(i).get();
			String genome=getGenome(allRAYTs.get(i).getIdent());
			for(int j=0;j<repintypes.size();j++) {
				String genomeRtype=genome+"\t"+repintypes.get(j);
				String rayt=allRAYTs.get(i).getIdent();
				if(repintypeDS.containsKey(genomeRtype) && repintypeDS.get(genomeRtype).length()>0) {
					repintypeDS.put(genomeRtype, repintypeDS.get(genomeRtype)+","+rayt);
				}else {
					repintypeDS.put(genomeRtype,rayt);
				}
			}
		}
		return repintypeDS;
	}
	
	private HashMap<String/*repintype\tgenome*/,/*rayt1,rayt2...*/String> initREPINTypes(ArrayList<String> genomeIDs,int maxTypes){
		HashMap<String/*repintype\tgenome*/,/*rayt1,rayt2...*/String> rtypes=new HashMap<String, String>();
		for(int i=0;i<genomeIDs.size();i++) {
			for(int j=0;j<maxTypes;j++) {
				rtypes.put(genomeIDs.get(i)+"\t"+j, "");
			}
		}
		return rtypes;
	}
	
	
	public void writeREPINType(File out,ArrayList<String> genomeIDs,int maxTypes) {
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write("genome\trepintype\trayts\traytIDs\n");
			
			HashMap<String/*genome\trepintype*/,/*rayt1,rayt2...*/String> repintypeDS=initREPINTypes(genomeIDs,maxTypes);
			repintypeDS=addRAYTsToREPINTypeDS(repintypeDS);
			String[] repintypes=repintypeDS.keySet().toArray(new String[0]);
			for(int i=0;i<repintypes.length;i++) {
				String rayts=repintypeDS.get(repintypes[i]);
				int raytNumber=0;
				if(rayts.length()>0) {
					String split[]=rayts.split(",");
					raytNumber=split.length;
				}else {
					rayts="NA";
				}
				bw.write(repintypes[i]+"\t"+raytNumber+"\t"+rayts+"\n");
			}
			bw.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}
	
	
	private ArrayList<Vicinity> getVicinityInformation(ArrayList<Info> raytPos,ArrayList<REPINGenomePositions> rgp,String genomeID){
		ArrayList<Vicinity> repinVicinity=new ArrayList<REPIN_RAYT_prox.Vicinity>();
		for(int i=0;i<raytPos.size();i++) {
			Vicinity vic=new Vicinity();
			for(int j=0;j<repinGroups;j++) {
				String gID=genomeID+"_"+j;
				File localOutFolder=new File(outFolder+"/"+gID);
				Info rayt=raytPos.get(i);
				if(isInVicinityInformation(rayt, gID, localOutFolder, rgp.get(j))) {
					vic.add(j);
				}
			}
			repinVicinity.add(vic);
		}
		return repinVicinity;
	}
	
	private boolean isInVicinityInformation(Info raytPos,String genomeID,File outFolder,REPINGenomePositions rgp){
		File ss=new File(outFolder+"/"+genomeID+".ss");
		if(ss.exists()) {
			ArrayList<String> allREPINs=toArrayList(getAllREPINs(ss));
			HashMap<String,ArrayList<REPINposition>> repinPos=rgp.get();
			return(isNear(allREPINs,raytPos,repinPos));
		}else {
			return false;
		}
	}
	
	private static ArrayList<Fasta> makeRAYTFasta(ArrayList<Info> inf/*blast result*/,File in/*genome Fasta sequence*/) {
		ArrayList<Fasta> seqs=new ArrayList<Fasta>();/*result sequences*/
		String genomeID=RAREFAN_MAIN.getGenomeID(in);
		for(int i=0;i<inf.size();i++) {

			HashMap<String,String> fas=Fasta.fasToHash(Fasta.readFasta(in), false);
			String fasIdent=inf.get(i).info.split("---")[0];

			

			String name=genomeID+"_"+i;
			int start=inf.get(i).getStart();
			int end=inf.get(i).getEnd();
			boolean rev=false;
			if(start>end) {
				int temp=start;
				start=start<end?start:end;
				end=end>start?end:temp;
				rev=true;
			}
			System.out.println(fasIdent+" "+in);
			System.out.println(start+" "+end);
			String seq=getSeq(fas.get(fasIdent),start,end,rev);
		
			seqs.add(new Fasta(name,seq));


		}
		return seqs;
	}
	private static String getSeq(String genomeSeq,int start,int end,boolean rev) {
		String seq=genomeSeq.substring(start,end);

		if(rev) {
			seq=DNAmanipulations.reverse(seq);
			if(seq.startsWith("T")) {
				seq=genomeSeq.substring(start,end+1);
				seq=DNAmanipulations.reverse(seq);
			}
		}else {
			seq=genomeSeq.substring(start-1,end);

		}
		return seq;
	}
	
//We may need to replace this function by a function that actually sorts RAYTs into groups depending on which REPINs they are associated with	
	public void addRAYTREPINProximity(int focalSeedGroup,String genomeID,File outFolder,REPINProperties rp,ArrayList<Info> raytPos) {
		File mcl=new File(outFolder+"/"+genomeID+".mcl");
		File ss=new File(outFolder+"/"+genomeID+".ss");
		File largestClusterF=new File(outFolder+"/"+genomeID+"_largestCluster.ss");
		if(mcl.exists()) {
			HashSet<String> all=getAllREPINs(ss);
			HashSet<String> largestCluster=getAllREPINs(largestClusterF);
			ArrayList<ArrayList<String>> clusters=getClusters(mcl);
			remove(all,clusters);
			HashSet<String> rest=all;
			clusters.add(toArrayList(rest));
			HashMap<String,ArrayList<REPINposition>> repinPos=rp.getREPINPositions();
			try {
				File out=new File(outFolder+"/"+genomeID+"_rayt_repin_prox.txt");
				BufferedWriter bw=new BufferedWriter(new FileWriter(out));
				int firstREPINCluster=-1;
				for(int j=0;j<clusters.size();j++) {
					if(j<clusters.size()-1) {
						String isREPINC=isREPIN(clusters.get(j).get(0));
						bw.write("\t"+j+isREPINC);
						if(isREPINC.equals("repin")&&firstREPINCluster==-1) {
							firstREPINCluster=j;
						}
					}
					else bw.write("\toutsideClusters");
				}
				if(firstREPINCluster==-1) {
					firstREPINCluster=0;
				}
				bw.write("\n");
				HashSet<Integer> clustermap=new HashSet<Integer>();
				int lc=0;  
				for(int i=0;i<raytPos.size();i++) {
					int start=raytPos.get(i).getStart();
					int end=raytPos.get(i).getEnd();
					bw.write(raytPos.get(i).info+"_"+start+"_"+end);
					boolean present=false;
					for(int j=0;j<clusters.size();j++) {
						if(isNear(clusters.get(j),raytPos.get(i),repinPos)) {
							bw.write("\t1");
							if(!present) {
								present=true;
								clustermap.add(j);
								if(largestCluster.contains(clusters.get(j).get(firstREPINCluster))) {
									lc++;
								}
							}else {
								System.err.println("Two clusters flanking RAYT in "+genomeID);
							}
						}else {
							bw.write("\t0");
						}
					}
					bw.write("\n");

				}
				stats.add(new ProxStats(genomeID,focalSeedGroup, raytPos.size(), clustermap.size(),lc));
				bw.close();
			}catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
	}
	
	private static HashSet<String> getAllREPINs(File in){
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		HashSet<String> allREPINs=new HashSet<String>();
		for(int i=0;i<fas.size();i++) {
			allREPINs.add(fas.get(i).getSequence());
		}
		return allREPINs;
	}
	
	private void remove(HashSet<String> all,ArrayList<ArrayList<String>> clusters){
		for(int i=0;i<clusters.size();i++) {
			for(int j=0;j<clusters.get(i).size();j++) {
				all.remove(clusters.get(i).get(j));
			}
		}
	}
	
	private static ArrayList<String> toArrayList(HashSet<String> rest){
		String[] list=rest.toArray(new String[0]);
		ArrayList<String> al=new ArrayList<String>();
		for(int i=0;i<list.length;i++) {
			al.add(list[i]);
		}
		return al;
	}
	
	public void writeStats(File out) {
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write(ProxStats.printHeading());
			for(int i=0;i<stats.size();i++) {
				bw.write(stats.get(i).print());
			}
			bw.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private String isREPIN(String repin) {
		if(repin.endsWith(makeAs(repin.length()/2))) {
			return "rep";
		}else {
			return "repin";
		}
	}
	
	private String makeAs(int length) {
		StringBuffer s=new StringBuffer();
		for(int i=0;i<length;i++) {
			s.append("A");
		}
		return s.toString();
	}
	
	private ArrayList<ArrayList<String>> getClusters(File in){
		ArrayList<ArrayList<String>> clusters=new ArrayList<ArrayList<String>>();
		try {
			
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null) {
				String repins[]=line.split("\\s+");
				ArrayList<String> cluster=new ArrayList<String>();
				for(int i=0;i<repins.length;i++) {
					cluster.add(repins[i].split("_")[0]);
				}
				clusters.add(cluster);
			}
			br.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return clusters;
	}
	private boolean isNear(ArrayList<String> repins,Info raytPos,HashMap<String,ArrayList<REPINposition>> repinPos) {
		ArrayList<Info> repinInfo=getREPINPos(repinPos,repins);
		for(int i=0;i<repinInfo.size();i++) {
			int repinstart=repinInfo.get(i).getStart();
			int repinend=repinInfo.get(i).getEnd();
			if(repinstart>repinend) {
				int swap=repinend;
				repinend=repinstart;
				repinstart=swap;
			}
			Info bigREPIN=new Info(repinstart-vicinityDistance,repinend+vicinityDistance,"repin");
			if(raytPos.overlapsWith(bigREPIN)) {
				return true;
			}
		}
		return false;
	}
	
	private ArrayList<Info> getREPINPos(HashMap<String,ArrayList<REPINposition>> repinPos,ArrayList<String> repins){
		ArrayList<Info> repinInf=new ArrayList<Info>();
		for(int i=0;i<repins.size();i++) {
			repinInf.addAll(getREPINPos(repinPos,repins.get(i)));
		}
		return repinInf;
	}
	
	private ArrayList<Info> getREPINPos(HashMap<String,ArrayList<REPINposition>> repinPos,String repin){
		ArrayList<Info> repinInf=new ArrayList<Info>();
		ArrayList<REPINposition> repList=repinPos.get(repin);
		for(int i=0;i<repList.size();i++) {
			int start=repList.get(i).start;
			int end=repList.get(i).end;
			repinInf.add(new Info(start,end,repin+"_fastaPos_"+repList.get(i).id));
		}
		return repinInf;
	}

}
