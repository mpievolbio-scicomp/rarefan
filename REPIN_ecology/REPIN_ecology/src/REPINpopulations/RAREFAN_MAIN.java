package REPINpopulations;
import java.io.*;
import java.util.*;
import frequencies.REPINProperties;
import identifyRAYTs.BlastRAYTs;
import util.*;
//the idea is to determine the REPIN populations that are present in a number of focal strains
//(runnning word freqs, identify groups etc) at first I will supply the sequence seeds (focal sequences)
//for each sequence group we determine the frequency in each of the given strains
//this data will then be displayed on a tree using R
import util.phylogenetics.RunTreePrograms;

public class RAREFAN_MAIN {

    //requires mcl, andi, clustDist and BLAST+
    String focalSeeds[];
    ArrayList<File> genomes;
    File inFolder;
    int numMuts=1;
    double minFrac=0.01;

    //distance from repin to rayt, if within vicinity then repin cluster is associated with that rayt
    String legacyBlastPerlLocation;
    File queryRAYT;
    File genomeFolder;
    String e;
    boolean analyseREPIN;
    File outFolder;
    HashMap<String/*genomes*/,HashMap<String/*focal seed*/,Integer/*pop size*/>> results=new HashMap<String,HashMap<String,Integer>>();
    public static HashSet<String> fastaExtensions=new HashSet<String>(Arrays.asList("fas","fasta","fna","fastn","fn"));
    int MCLThreads=1;

	  // Entry point.
    public static void main(String args[]) {

        // Handle wrong number of arguments.
        if(args.length<11 || args.length>12) {
            System.out.println("Usage: java -jar REPIN_ecology.jar IN_DIR OUT_DIR REFERENCE_STRAIN NMER_OCCURENCE MIN_NMER_LENGTH QUERY_RAYT TREEFILE E_VALUE_CUTOFF ANALYZE_REPINS MCL_THREADS DISTANCE_GROUP_SEEDS [PATH_TO_LEGACY_BLAST.PL]");

            System.exit(1);
        }

        File inFolder=new File(args[0]);
        File outFolder=new File(args[1]);
        String focalSeedGenome=args[2];
        int minRepFreq=Integer.parseInt(args[3]);
        int wordlength=Integer.parseInt(args[4]);
        File queryRAYT=new File(args[5]);
        File treeFile=new File(args[6]);
        String evalue=args[7];
        boolean analyseREPIN=args[8].equalsIgnoreCase("true");
        int MCLThreads=Integer.parseInt(args[9]);
        int distanceGroupSeeds=Integer.parseInt(args[10]);

        File out=new File(outFolder+"/results.txt");
        RAREFAN_MAIN dpf;
        String program="tblastn";

        // legacy_blast path not given.
        String legacyBlastPerlLocation="";

        // legacy_blast path given.
        if(args.length==12) {
            legacyBlastPerlLocation=args[11];
        }
        dpf=new RAREFAN_MAIN(inFolder, outFolder,focalSeedGenome,minRepFreq,wordlength,queryRAYT,program,treeFile,legacyBlastPerlLocation,evalue,analyseREPIN,MCLThreads,distanceGroupSeeds);
        dpf.print(out);
    }
	
    // Workhorse function.
    public RAREFAN_MAIN(File inFolder,File outFolder,String focalSeedGenome,int minRepFreq,int wordlength,File queryRAYT,String program,File treeFile,String legacyBlastPerlLocation,String evalue,boolean analyseREPIN,int MCLThreads,int distanceGroupSeeds){
        this.inFolder=inFolder;
        this.outFolder=outFolder;
        this.MCLThreads=MCLThreads;
        outFolder.mkdirs();
        genomes=getFiles();
        this.legacyBlastPerlLocation=legacyBlastPerlLocation;
        this.queryRAYT=queryRAYT;
        this.focalSeeds=getFocalSeeds(focalSeedGenome,minRepFreq,wordlength,distanceGroupSeeds);
        this.genomeFolder=inFolder;
        this.analyseREPIN=analyseREPIN;
        e=evalue;
        calculateResults();
        BlastRAYTs.runProgram(inFolder, queryRAYT, outFolder, e, program, getREPtype(), "yafM_relatives.fna",analyseREPIN);
        // 		treeFile=new File(outFolder+"/"+treeFile);
        // 		if(!treeFile.exists()) {
        // 			generateTree(treeFile);
        // 		}
    }

    private void generateTree(File treeFile) {
      
        System.out.println("Generating Tree.");

        String filenames=generateFileNameString();
        String treeID=treeFile.getName().split("\\.")[0];
        File distFile=new File(outFolder+"/"+treeID+".dist");
        System.out.println("Running andi.");
        RunTreePrograms.runProgram("andi "+filenames, "", outFolder,distFile);
        System.out.println("Running clustDist.");
        RunTreePrograms.runProgram("clustDist "+distFile, "", outFolder, treeFile);
    }

    private String generateFileNameString() {
        StringBuffer sb=new StringBuffer();
        File[] files=inFolder.listFiles();
        for(int i=0;i<files.length;i++) {
            if(hasCorrectExtension(files[i])) {
                sb.append(" "+files[i]);
            }
        }
        return sb.toString();
    }

    private String[] getREPtype() {
        ArrayList<String> list=new ArrayList<String>();
        for(int i=0;i<focalSeeds.length;i++) {
            list.add(i+"");
        }
        return list.toArray(new String[0]);
    }

    private String[] getFocalSeeds(String genome,int minRepFreq,int wl,int distanceGroupSeeds) {
        File fsg=new File(inFolder+"/"+genome);
        DetermineFocalSeeds dfs=new DetermineFocalSeeds(fsg,outFolder,minRepFreq,wl,distanceGroupSeeds);
        return dfs.getFocalSeeds();
    }

    public void print(File out) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(out));
            String[] genomes=results.keySet().toArray(new String[0]);
            for(int i=0;i<genomes.length;i++) {
                String[] seeds=results.get(genomes[i]).keySet().toArray(new String[0]);
                for(int j=0;j<seeds.length;j++) {
                    bw.write(genomes[i].replace("_", "\t")+"\t"+seeds[j]+"\t"+results.get(genomes[i]).get(seeds[j])+"\n");
                }
            }
            bw.close();

        }catch(IOException e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }
    public static String getGenomeID(File in) {
        return in.getName().split("\\.")[0];
    }
    private void calculateResults() {
        System.out.println("Calculating Results.");

        REPIN_RAYT_prox rrp=new REPIN_RAYT_prox(this.outFolder,focalSeeds.length);
        ArrayList<String> genomeIDs=new ArrayList<String>();
        // TODO: Can we parallelize this loop?
        for(int i=0;i<genomes.size();i++) {
            String onlyGenome=getGenomeID(genomes.get(i));
            // parallelize?
            genomeIDs.add(onlyGenome);
            ArrayList<REPINGenomePositions> rgp=new ArrayList<REPINGenomePositions>();
            ArrayList<Info> raytPos=writeRAYTLocation(genomes.get(i));

            for(int j=0;j<focalSeeds.length;j++) {
                String genomeID=onlyGenome+"_"+j;

                results.put(genomeID, new HashMap<String,Integer>());

                File outFolder=new File(this.outFolder+"/"+genomeID+"/");
                outFolder.mkdir();
                int wl=focalSeeds[j].length();
        
                REPINProperties rp=new REPINProperties(outFolder,genomeID,genomes.get(i),wl,numMuts,minFrac,null,focalSeeds[j],false,analyseREPIN,MCLThreads);
                System.out.println("Write REPINs as artemis files for "+genomeID+"...");

                writeREPINArtemis(new File(outFolder+"/"+genomeID+"_largestCluster.ss"),j);
                writeREPINArtemis(new File(outFolder+"/"+genomeID+".ss"),j);
                File cluster;
                int k=0;
                while((cluster=new File(outFolder+"/"+genomeID+"_"+k+".ss")).exists()){
                    writeREPINArtemis(cluster,k);
                    k++;
                }
                System.out.println("Write RAYT locations "+genomeID+"...");

                int popsize=rp.getPopSize();
                results.get(genomeID).put(focalSeeds[j],popsize);
                rgp.add(new REPINGenomePositions(rp.getREPINPositions()));
            }
            System.out.println("REPIN RAYT proximity calculation for "+onlyGenome+"...");

            rrp.addRAYT(raytPos, genomes.get(i), rgp);
            System.out.println("rgp: "+rgp.size()+"\n rrp: "+rrp.allRAYTs.size());

        }
        rrp.write(new File(outFolder+"/repin_rayt_association.txt"));
        rrp.writeREPINType(new File(outFolder+"/repin_rayt_association_byREPIN.txt"),genomeIDs,focalSeeds.length);
    }


    private ArrayList<Info> writeRAYTLocation(File genome) {
        String genomeID=getGenomeID(genome);
        ArrayList<Info> RAYTLocations;
        if(legacyBlastPerlLocation!="") {
            RAYTLocations=BlastRAYTs.blastQuery(genome, queryRAYT, outFolder, e, "tblastn",legacyBlastPerlLocation);
        }else {
            RAYTLocations=BlastRAYTs.blastQuery(genome, queryRAYT, outFolder, e, "tblastn");
        }
        WriteArtemis.write(RAYTLocations, new File(outFolder+"/rayt_"+genomeID+".tab"));
        return RAYTLocations;
    }


    private void writeREPINArtemis(File in,int group) {
        if(in.exists()) {
            ArrayList<Fasta> fas=Fasta.readFasta(in);
            ArrayList<Info> pos=new ArrayList<Info>();
            for(int i=0;i<fas.size();i++) {
                String ident=fas.get(i).getIdent();
                pos.addAll(getPos(ident,"Gr_"+group));
            }
            String genomeID=getGenomeID(in);
            WriteArtemis.write(pos, new File(in.getParent()+"/"+genomeID+".tab"));

        }
    }

    private ArrayList<Info> getPos(String ident,String inf){
        String split[]=ident.split("\\s+");
        ArrayList<Info> pos=new ArrayList<Info>();
        for(int i=1;i<split.length;i++) {
            String[] split2=split[i].split("_");
            int start=Integer.parseInt(split2[1]);
            int end=Integer.parseInt(split2[2]);
            pos.add(new Info(start,end,inf+" "+split[0]));
        }
        return pos;
    }
	
    public static  boolean hasCorrectExtension(File f) {
        String[] split=f.getAbsolutePath().split("\\.");
        String ext=split[split.length-1];
        if(fastaExtensions.contains(ext)) {
            return true;
        }
        return false;
    }
	
    private ArrayList<File> getFiles() {
        ArrayList<File> genomes=new ArrayList<File>();
        File[] all=inFolder.listFiles();
        for(int i=0;i<all.length;i++) {
			  
            if(hasCorrectExtension(all[i])) {
                genomes.add(all[i].getAbsoluteFile());
            }
        }
        return genomes;
    }

}
