package blastTools;

import java.io.*;
import java.util.*;

import util.Info;
import util.InfoTree;

public class GetBestOverlappingHits {
	public static void main(String args[]){
		File queryFolder=new File(args[0]);
		File dbFolder=new File(args[1]);
		String program=args[2];
		double eValue=Double.parseDouble(args[3]);
		File outFolder=new File(args[4]);
		String BLASTBINPATH=args[5];
		File[] dbFiles=dbFolder.listFiles();
		File[] queryFiles=queryFolder.listFiles();
		for(int i=0;i<dbFiles.length;i++){
			if(!dbFiles[i].getAbsolutePath().endsWith(".fasta"))continue;
			for(int j=0;j<queryFiles.length;j++){
				if(!queryFiles[j].getAbsolutePath().endsWith(".fasta"))continue;
				String query=queryFiles[j].getName().split("\\.")[0];
				String db=dbFiles[i].getName().split("\\.")[0];
				File temp=new File(outFolder+"/temp.txt");
				boolean DNA=program.endsWith("blastn");
				PerformBlast.blast(BLASTBINPATH+"/blastall",BLASTBINPATH+"/formatdb",program, eValue, temp, queryFiles[j], dbFiles[i],true, false,DNA,false);
				File output=new File(outFolder+"/"+db+"_"+query+".txt");
				ReadBlast rb=new ReadBlast(temp);
				InfoTree it=new InfoTree();
				for(int k=0;k<rb.getStartDB().size();k++){
					int start=Math.min(rb.getStartDB().get(k),rb.getEndDB().get(k));
					int end=Math.max(rb.getStartDB().get(k),rb.getEndDB().get(k));
					Info currentInf=new Info(start,end,rb.getEvalue().get(k)+"");
					
					if(!it.checkOverlap(currentInf)){
						it.insert(currentInf);
						//System.out.println(currentInf);
						
					}else{
						ArrayList<Info> al=new ArrayList<Info>();
						it.search(currentInf, al);
						if(al.size()==1){
							Info resident=al.get(0);
							//System.out.println(resident);
							double res=Double.parseDouble(resident.getInfo());
							double current=Double.parseDouble(currentInf.getInfo());
							if(res>current){
								it.remove(resident);
								
								it.insert(currentInf);
							}
						}else {
							
							for(int l=0;l<al.size();l++){
								it.remove(al.get(l));
							}
							it.insert(currentInf);
						}
						
					}
					
				}
				ArrayList<Info> al=it.parseTree();
				try{
					BufferedWriter bw=new BufferedWriter(new FileWriter(output));
					for(int k=0;k<al.size();k++){
						bw.write(al.get(k)+"\n");
					}
					bw.close();
				}catch(IOException e){
					e.printStackTrace();
				}
			}
		}
	}
}
