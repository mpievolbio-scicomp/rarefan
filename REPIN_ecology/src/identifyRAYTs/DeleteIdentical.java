package identifyRAYTs;

import java.io.*;
import java.util.*;

import util.*;

public class DeleteIdentical {
	public static void main(String args[]) {
		File in=new File(args[0]);
		//File out=new File(in.getParentFile()+"/AsAdded_"+in.getName());
		File outUnique=new File(in.getParentFile()+"/Unique_"+in.getName());

		//addAs(in,out);
		deleteIdentical(in,outUnique);
	}
	
	public static void deleteIdentical(File in,File out) {
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		HashMap<String,ArrayList<String>> seqs=new HashMap<String, ArrayList<String>>();
		for(int i=0;i<fas.size();i++) {
			String id=getGenomeID(fas.get(i).getIdent());
			if(!seqs.containsKey(id)) {
				seqs.put(id, new ArrayList<String>());
			}
			seqs.get(id).add(fas.get(i).getSequence());
		}
		ArrayList<Fasta> fasOut=deleteIdentical(seqs);
		Fasta.write(fasOut, out);
	}
	
	public static ArrayList<Fasta> deleteIdentical(HashMap<String,ArrayList<String>> seqs){
		String[] ids=seqs.keySet().toArray(new String[0]);
		ArrayList<Fasta> deleted=new ArrayList<Fasta>();
		int count=0;
		for(int i=0;i<ids.length;i++) {
			ArrayList<String> list=seqs.get(ids[i]);
			HashSet<String> unique=new HashSet<String>();
			for(int j=0;j<list.size();j++) {
				if(!unique.contains(list.get(j))) {
					unique.add(list.get(j));
				}
			}
			if(unique.size()>1) {
				count++;
			}
			add(deleted,unique.toArray(new String[0]),ids[i]);
		}
		System.out.println(count+" genomes with more than one RAYT still present.");
		return deleted;
	}
	
	public static void add(ArrayList<Fasta> deleted,String[] unique,String id) {
		for(int i=0;i<unique.length;i++) {
			deleted.add(new Fasta(id+"_"+i,unique[i]));
		}
	}
	
	public static String getGenomeID(String fasID) {
		String split[]=fasID.split("\\s+");
		String[] split2=split[0].split("_");
		String id=split2[0];
		for(int i=1;i<split2.length-2;i++) {
			id=id+"_"+split2[i];
		}
		return id;
	}
	
	public static void addAs(File in,File out) {
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		for(int i=0;i<fas.size();i++) {
			if(!fas.get(i).getSequence().startsWith("A")) {
				fas.get(i).setSequence("A"+fas.get(i).getSequence().substring(1));
			}
		}
		Fasta.write(fas, out);
	}
	
}
