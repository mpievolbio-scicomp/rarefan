package identifyRAYTs;

import java.io.*;
import java.util.*;

import blastTools.*;
import util.*;


public class BlastRAYTs {
	static int add=0;
	public static void main (String args[]) {
		File inFolder=new File(args[0]);
		File query=new File(args[1]);
		File outFolder=new File(args[2]);
		String e="1e-100";
		String program=args[3];
		String nameSeqs=args[4];
		String[] repType=Arrays.copyOfRange(args,5,args.length);

		runProgram(inFolder,query,outFolder,e,program,repType,nameSeqs);
	}
	static int minClusterSize=10;
	public static void runProgram(File inFolder,File query,File outFolder,String e,String program,String[] repType,String nameSeqs) {
		for(int k=0;k<repType.length;k++) {
			ArrayList<Fasta> seqs=new ArrayList<Fasta>();
			File out=new File(outFolder+"/"+nameSeqs);
			File presAbs=new File(outFolder+"/presAbs_"+repType[k]+".txt");
			File maxREPINOut=new File(outFolder+"/maxREPIN_"+repType[k]+".txt");
			HashMap<String,String> presAbsHash=new HashMap<String,String>();
			File[] dbs=inFolder.listFiles();

			for(int i=0;i<dbs.length;i++) {
				if(dbs[i].getName().endsWith("fas")||dbs[i].getName().endsWith("fna")) {
					File db=dbs[i];
					ArrayList<Info> bi=blastQuery(db, query, outFolder, e,program);
					String seqName=dbs[i].getName().split("\\.")[0];
					String maxREPIN=getMaxREPIN(seqName,inFolder,repType[k]);
					System.out.println(seqName+" "+inFolder+" "+repType[k]);
					if(!maxREPIN.equals("-1")){
						System.out.println(maxREPIN);
						int maxREPINNum=Integer.parseInt(maxREPIN.split("_|\\s+")[1]);
						
						int repins=Integer.parseInt(maxREPIN.split("_")[2]);
						int allREPs=getREPNumbers(seqName, inFolder,repType[k]);
						int numREPINClusters=getNumREPINClusters(seqName,inFolder,repType[k]);
						presAbsHash.put(seqName, bi.size()+"\t"+repins+"\t"+maxREPIN.split("_")[0]+"\t"+maxREPINNum+"\t"+allREPs+"\t"+numREPINClusters);
						print(bi,dbs[i],seqs);
					}
				}
			}
			printHM(presAbsHash,presAbs,maxREPINOut);
			Fasta.write(seqs, out);
		}

	}
	
	public static int getNumREPINClusters(String name,File folder,String repType) {
		try {
			File in=new File(folder+"/"+name+"_"+repType+"/"+name+"_"+repType+".mcl");
			//System.out.println(in);
			int numClusters=0;
			if(in.exists()) {
				BufferedReader br=new BufferedReader(new FileReader(in));
				String line="";
				while((line=br.readLine())!=null) {
					int sum=0;
					String[] split=line.split("\\s+");
					for(int i=0;i<split.length;i++) {
						String[] splitREPIN=split[i].split("_");
						String REPIN=splitREPIN[0];
						int num=Integer.parseInt(splitREPIN[1]);
						if(isREPIN(REPIN)) {
							sum+=num;
						}else {
							break;
						}
					}
					if(sum>=minClusterSize) {
						numClusters++;
					}
				}
				br.close();
				return numClusters;
			}else {
				return -1;
			}
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return -1;
	}

	
	public static String getMaxREPIN(String name,File folder,String repType) {
		try {
			File in=new File(folder+"/"+name+"_"+repType+"/"+name+"_"+repType+"_largestCluster.nodes");
			if(!in.exists())System.err.println(in);
			if(in.exists()) {
				BufferedReader br=new BufferedReader(new FileReader(in));
				String line="";
				int max=0;
				//String maxrepin="-1";
				String all="";
				int sum=0;
				while((line=br.readLine())!=null) {

					String[] split=line.split("_|\\s+");
					
					int num= Integer.parseInt(split[1]);
					String curr=split[0];
					if(isREPIN(curr)) {
						sum+=num;
						if(num>max) {
							max=num;
							//maxrepin=curr;
							all=split[0]+"_"+split[1];
						}
					}
				}
				br.close();
				if(all.length()>1 && sum>0) {
					return all+"_"+sum;
				}else return "-1";
			}else {
				return "-1";
			}
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return "-1";
	}
	
	private static boolean isREPIN(String repin) {
		int length=repin.length();
		if(repin.substring(length/2).equals(makeREPIN('A',length/2))) {
			return false;
		}
		return true;
	}
	
	private static String makeREPIN(char nuc,int rep) {
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<rep;i++) {
			sb.append(nuc);
		}
		return sb.toString();
	}
	
	public static int getREPNumbers(String name,File folder,String repType) {
		try {
			File in=new File(folder+"/"+name+"_"+repType+"/"+name+"_"+repType+".nodes");
			//System.out.println(in);
			if(in.exists()) {
				BufferedReader br=new BufferedReader(new FileReader(in));
				String line="";
				int reps=0;
				while((line=br.readLine())!=null) {
				  reps+=Integer.parseInt(line.split("\\s+")[1]);
				}
				br.close();
				return reps;
			}else {
				return -1;
			}
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return -1;
	}
	
	public static void printHM(HashMap<String,String> hm,File out,File maxREPINOut) {
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			BufferedWriter bwMaxREPIN=new BufferedWriter(new FileWriter(maxREPINOut));
			String[] keys=hm.keySet().toArray(new String[0]);
			for(int i=0;i<keys.length;i++) {
				bw.write(keys[i]+"\t"+hm.get(keys[i])+"\n");
				bwMaxREPIN.write(">"+keys[i]+"\n"+hm.get(keys[i]).split("\t")[2]+"\n");
			}
			bw.close();
			bwMaxREPIN.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	public static String getName(String fasIdent) {
		String split[]=fasIdent.split("\\s+");
		String name=split[2].substring(0,3);
	    if(fasIdent.contains("chromosome")){
	        name=name+split[split.length-4];
	      }else{
		        name=name+split[split.length-3];
	      }
	      name=name.replace(".","");
	      name=name.replace("*","");
	      name=name.replace(",","");
	      
	      return name;
	}
	
	public static void print(ArrayList<Info> inf,File in,ArrayList<Fasta> seqs) {
		for(int i=0;i<inf.size();i++) {
			//System.out.println(inf.get(i).toString());
			ArrayList<Fasta> fas=Fasta.readFasta(in);
			String fasIdent=fas.get(0).getIdent();
			//String shortName=getName(fasIdent);
			

			String name=fasIdent.split("\\s+")[0]+"_"+inf.get(i).getStart()+"_"+inf.get(i).getEnd();
			int start=inf.get(i).getStart();
			int end=inf.get(i).getEnd();
			boolean rev=false;
			if(start>end) {
				int temp=start;
				start=start<end?start:end;
				end=end>start?end:temp;
				rev=true;
			}
			String seq=getSeq(fas,start,end,rev);
		
			seqs.add(new Fasta(name+" "+inf.get(i).toString(),seq));
			

		}
	}
	
	private static String getSeq(ArrayList<Fasta> fas,int start,int end,boolean rev) {
		String seq=fas.get(0).getSequence().substring(start-add,end+add);

		if(rev) {
			seq=DNAmanipulations.reverse(seq);
			if(seq.startsWith("T")) {
				seq=fas.get(0).getSequence().substring(start-add,end+add+1);
				seq=DNAmanipulations.reverse(seq);
			}
		}else {
			seq=fas.get(0).getSequence().substring(start-add-1,end+add);

		}
		return seq;
	}
	
	public static ArrayList<Info> blastQuery(File db, File query,File outFolder,String e,String program){
		return blastQuery(db, query, outFolder, e, program, "/usr/local/share/");
	}
	
	public static ArrayList<Info> blastQuery(File db, File query,File outFolder,String e,String program,String legacyBlastPerlLocation){
		File blastFolder=new File(outFolder+"/blastout/");
		blastFolder.mkdir();
		File out=new File(blastFolder+"/"+db.getName()+".blast");
		ArrayList<Info> blastIntervals=new ArrayList<Info>();
		//out.deleteOnExit();
	
		PerformBlast.blast(legacyBlastPerlLocation+"legacy_blast.pl blastall",legacyBlastPerlLocation+"legacy_blast.pl formatdb",program, Double.parseDouble(e), out, query, db, false,false,true,false);
		ReadBlast rb=new ReadBlast(out);
		int querylength=Fasta.readFasta(query).get(0).getSequence().length();
		for(int i=0;i<rb.getDatabase().size();i++){
			int start=rb.getStartDB().get(i);
			int end=rb.getEndDB().get(i);
			//int temp=start;
			//start=start<end?start:end;
			//end=end>start?end:temp;
			
			if(Math.abs(start-end)>240) {
				int multi=program.equals("blastn")?1:3;
				int adjstart=end>start?start-rb.getStartQuery().get(i)*multi:start+rb.getStartQuery().get(i)*multi;
				int adjend=end>start?end+(querylength-rb.getEndQuery().get(i))*multi:end-(querylength-rb.getEndQuery().get(i))*multi;
				System.out.println(adjstart-adjend);
				blastIntervals.add(new Info(adjstart,adjend,rb.getDatabase().get(i)+"---"+rb.getEvalue().get(i)));
			}
		}
		return blastIntervals;
	}

}