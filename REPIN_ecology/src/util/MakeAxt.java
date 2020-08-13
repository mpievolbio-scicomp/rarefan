package util;

import java.io.*;
import java.util.*;

public class MakeAxt {
	public static void main(String args[]) {
		File in=new File(args[0]);
		
		makeAxtFolder(in);
	}
	
	public static void makeAxtFolder(File folder) {
		File[] files=folder.listFiles();
		for(int i=0;i<files.length;i++) {
			if(files[i].getName().endsWith("fas")||files[i].getName().endsWith("fasta")) {
				makeAxt(files[i]);
			}
		}
	}
	
	public static void makeAxt(File in) {
		ArrayList<Fasta> fas=Fasta.readFasta(in);
		File out=new File(in.getParent()+"/"+in.getName().split("\\.")[0]+".axt");
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<fas.size();i++) {
				for(int j=i+1;j<fas.size();j++) {
					bw.write(fas.get(i).getIdent().split("\\s+")[0]+" "+fas.get(j).getIdent().split("\\s+")[0]+"\n");
					bw.write(fas.get(i).getSequence().replace("-", "")+"\n");
					bw.write(fas.get(j).getSequence().replace("-","")+"\n");
					bw.write("\n");
				}
			}
			bw.close();
		}catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
