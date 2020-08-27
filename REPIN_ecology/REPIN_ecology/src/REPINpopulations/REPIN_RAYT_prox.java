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
import util.Fasta;
import util.Info;

public class REPIN_RAYT_prox {
	int vicinity=200;
	ArrayList<ProxStats> stats=new ArrayList<ProxStats>();

	
	
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
				for(int j=0;j<clusters.size();j++) {
					if(j<clusters.size()-1) {
						bw.write("\t"+j+isREPIN(clusters.get(j).get(0)));
						
					}
					else bw.write("\toutsideClusters");
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
								if(largestCluster.contains(clusters.get(j).get(0))) {
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
	
	private HashSet<String> getAllREPINs(File in){
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
	
	private ArrayList<String> toArrayList(HashSet<String> rest){
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
			Info bigREPIN=new Info(repinstart-vicinity,repinend+vicinity,"repin");
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
