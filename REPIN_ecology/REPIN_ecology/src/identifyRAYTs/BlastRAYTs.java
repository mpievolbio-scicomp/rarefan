package identifyRAYTs;

import java.io.*;
import java.util.*;

import REPINpopulations.RAREFAN_MAIN;
import blastTools.*;
import util.*;


public class BlastRAYTs {
	static int add=0;
	static int masterDist=3;
	public static void main (String args[]) {
		File inFolder=new File(args[0]);
		File query=new File(args[1]);
		File outFolder=new File(args[2]);
		String e="1e-100";
		String program=args[3];
		String nameSeqs=args[4];
		String[] repType=Arrays.copyOfRange(args,5,args.length);

		runProgram(inFolder,query,outFolder,e,program,repType,nameSeqs,true);
	}
	static int minClusterSize=10;
	

	
	public static void runProgram(File inFolder,File query,File outFolder,String e,String program,String[] repType,String nameSeqs,boolean analyseREPIN) {
			
		System.out.println("Running program " + program + ".");
		for(int k=0;k<repType.length;k++) {
			ArrayList<Fasta> seqs=new ArrayList<Fasta>();
			File out=new File(outFolder+"/"+nameSeqs);
			File presAbs=new File(outFolder+"/presAbs_"+repType[k]+".txt");
			File maxREPINOut=new File(outFolder+"/maxREPIN_"+repType[k]+".txt");
			HashMap<String,String> presAbsHash=new HashMap<String,String>();
			File[] dbs=inFolder.listFiles();

			for(int i=0;i<dbs.length;i++) {
				if(RAREFAN_MAIN.hasCorrectExtension(dbs[i])) {
					File db=dbs[i];
					ArrayList<Info> bi=blastQuery(db, query, outFolder, e,program);
					String seqName=getName(dbs[i]);
					String maxREPIN=getMaxREPIN(seqName,outFolder,repType[k],analyseREPIN);
					//System.out.println(seqName+" "+outFolder+" "+repType[k]);
					if(!maxREPIN.equals("-1")){
						//System.out.println(maxREPIN);
						int maxREPINNum=Integer.parseInt(maxREPIN.split("_|\\s+")[1]);

						int masterSeqs=Integer.parseInt(maxREPIN.split("_")[2]);
						int allREPs=getOnlyREPINNumbers(seqName, outFolder,repType[k],false);
						int allREPINs=getOnlyREPINNumbers(seqName, outFolder,repType[k],true&&analyseREPIN);
						int numREPINClusters=getNumREPINClusters(seqName,outFolder,repType[k],analyseREPIN);
						int numREPINDist=getNumREPINDist(seqName, outFolder,repType[k],masterDist,maxREPIN.split("_")[0],analyseREPIN);
						presAbsHash.put(seqName, bi.size()+"\t"+masterSeqs+"\t"+maxREPIN.split("_")[0]+"\t"+maxREPINNum+"\t"+allREPs+"\t"+numREPINClusters+"\t"+allREPINs+"\t"+numREPINDist);
						print(bi,dbs[i],seqs);
					}
				}
			}
			printHM(presAbsHash,presAbs,maxREPINOut);
			Fasta.write(seqs, out);
		}

	}
	
	private static String getName(File in) {
		String parts=RAREFAN_MAIN.getGenomeID(in);
		return parts;
	}
	
	public static int getNumREPINClusters(String name,File folder,String repType,boolean analyseREPIN) {
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
						if(isREPIN(REPIN)||!analyseREPIN) {
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
	
	public static String getMaxREPIN(String name,File folder,String repType,boolean analyseREPIN) {
		try {
			File in=new File(folder+"/"+name+"_"+repType+"/"+name+"_"+repType+"_largestCluster.nodes");
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
					if(isREPIN(curr)||!analyseREPIN) {
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
				System.err.println("Cannot find file "+in+". There is probably no REPIN in the corresponding submitted genome.");
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

	
	public static int getNumREPINDist(String name,File folder,String repType,int masterDist,String master,boolean analyseREPINs) {
		try {
			File in=new File(folder+"/"+name+"_"+repType+"/"+name+"_"+repType+".nodes");
			//System.out.println(in);
			if(in.exists()) {
				BufferedReader br=new BufferedReader(new FileReader(in));
				String line="";
				int reps=0;
				while((line=br.readLine())!=null) {
					String split[]=line.split("\\s+");
					String curr=split[0];
					String REPIN=curr.split("_")[0];
					int occ=Integer.parseInt(split[1]);
					if((isREPIN(REPIN)|| !analyseREPINs) && maxDistMaster(REPIN,master,masterDist)) {
						reps+=occ;
					}
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
	
	private static boolean maxDistMaster(String REPIN,String master,int maxDist) {
		String a=REPIN;
		String b=master;
		int differences=0;
		if(a.length()!=b.length()) {
			System.err.println("REPIN ("+REPIN+")is not the same length as master sequence ("+master+")");
			System.exit(-1);
		}
		for(int i=0;i<a.length();i++){
			if(a.charAt(i)!=b.charAt(i)){
				differences++;
			}
		}
		return differences<=maxDist;
	}
	
	public static int getOnlyREPINNumbers(String name,File folder,String repType,boolean analyseOnlyREPINs) {
		try {
			File in=new File(folder+"/"+name+"_"+repType+"/"+name+"_"+repType+".nodes");
			//System.out.println(in);
			if(in.exists()) {
				BufferedReader br=new BufferedReader(new FileReader(in));
				String line="";
				int reps=0;
				while((line=br.readLine())!=null) {
					String split[]=line.split("\\s+");
					String curr=split[0];
					String REPIN=curr.split("_")[0];
					int occ=Integer.parseInt(split[1]);
					if(isREPIN(REPIN)||!analyseOnlyREPINs) {
						reps+=occ;
					}
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
			bw.write("strain\tnumRAYTs\tNumREPINsLargestCluster\tmasterSequence\tmastersequenceFreq\tallREP\tnumberOfRepinClusters\tallREPINs\tREPINNumMasterDist_"+masterDist+"\n");
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
	
	public static void print(ArrayList<Info> inf/*blast result*/,File in/*genome Fasta sequence*/,ArrayList<Fasta> seqs/*result sequences*/) {
		for(int i=0;i<inf.size();i++) {

			HashMap<String,String> fas=Fasta.fasToHash(Fasta.readFasta(in), false);
			String fasIdent=inf.get(i).info.split("---")[0];

			

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
			//System.out.println(fasIdent+" "+in);
			//System.out.println(start+" "+end);
		
			String seq=getSeq(fas.get(fasIdent),start,end,rev);
		
			seqs.add(new Fasta(name+" "+inf.get(i).toString(),seq));

		}
	}
	
	private static String getSeq(String genomeSeq,int start,int end,boolean rev) {
		String seq=genomeSeq.substring(start-add,end+add);

		if(rev) {
			seq=DNAmanipulations.reverse(seq);
			if(seq.startsWith("T")) {
				seq=genomeSeq.substring(start-add,end+add+1);
				seq=DNAmanipulations.reverse(seq);
			}
		}else {
			seq=genomeSeq.substring(start-add-1,end+add);

		}
		return seq;
	}

	public static ArrayList<Info> blastQuery(File db, File query,File outFolder,String e,String program){
		return blastQuery(db, query, outFolder, e, program, "");
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
			HashMap<String,Fasta> fas=Fasta.fasToFastaHash(Fasta.readFasta(db), false);
			int seqlength=fas.get(rb.getDatabase().get(i)).getSequence().length();
			int start=rb.getStartDB().get(i);
			int end=rb.getEndDB().get(i);
			//int temp=start;
			//start=start<end?start:end;
			//end=end>start?end:temp;

			if(Math.abs(start-end)>240) {
				int multi=program.equals("blastn")?1:3;
				int adjstart=end>start?start-rb.getStartQuery().get(i)*multi:start+rb.getStartQuery().get(i)*multi;
				int adjend=end>start?end+(querylength-rb.getEndQuery().get(i))*multi:end-(querylength-rb.getEndQuery().get(i))*multi;
				//System.out.println(adjstart-adjend);
				if(adjend>seqlength) {
					adjend=seqlength;
				}
				if(adjstart<0) {
					adjstart=1;
				}
				blastIntervals.add(new Info(adjstart,adjend,rb.getDatabase().get(i)+"---"+rb.getEvalue().get(i)));
			}
		}
		return blastIntervals;
	}

}
