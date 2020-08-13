package identifyRAYTs;

import java.io.*;
import java.util.*;

import util.*;

public class Homopolymeric {
	public static void main(String args[]) {
		File in=new File(args[0]);
		File rayt=new File(args[1]);
		File out=new File(in.getParent()+"/homopoly.txt");
		findHomoPoly(in,rayt, out);
	}
	
	public static void findHomoPoly(File in,File rayt,File out) {
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		ArrayList<Fasta> raytfas=Fasta.readFasta(rayt);

		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write(makeCounts(fas).toString());
			bw.write(makeCounts(raytfas).toString());
			bw.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static  StringBuffer makeCounts(ArrayList<Fasta> fas) {
		String seq=fas.get(0).getSequence();
		String id=fas.get(0).getIdent();

		StringBuffer out=new StringBuffer();
		char[] c= {'T','G'};
		for(int j=0;j<c.length;j++) {
			for(int i=1;i<8;i++) {
				String query=getQuery(i,c[j]);
				int count=count(seq,query)+count(DNAmanipulations.reverse(seq),query);
				out.append(id.split("\\s+")[0]+"\t"+query+"\t"+count+"\n");
			}
		}
		return out;
	}
	
	
	private static String getQuery(int length,char c) {
		String query="";
		for(int i=0;i<length;i++) {
			query=query+c;
		}
		return query;
	}
	
	public static int count(String seq,String query) {
		int i=seq.indexOf(query);
		int count=0;
		while(i!=-1) {
			i=seq.indexOf(query, i+1);
				count++;
		}
		return count;
	}
	
}
