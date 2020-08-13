package ecoliREP;

import java.io.File;
import java.util.ArrayList;

import util.Fasta;
import util.Info;
import util.InfoTree;
import util.ReadGenbank;

public class AnnotateDeletions {
	public static void main(String args[]){
		File inFolder=new File(args[0]);
		File config=new File(args[1]);
		File outFolder=new File(args[2]);
		String reference=args[3];
		
		ArrayList<Infos> infos=Infos.readInfos(config);
		writeRecombinations(inFolder,infos,outFolder,reference);
	}
	
	public static void writeRecombinations(File inFolder,ArrayList<Infos> infos,File outFolder,String reference){

			for(int i=0;i<infos.size();i++){
				
				String name1=infos.get(i).name;
				if(!name1.equalsIgnoreCase(reference)){
					continue;
				}
				File currentFolder=new File(outFolder+"/"+name1+"/");
				if(!currentFolder.exists())currentFolder.mkdirs();
				ReadGenbank rgbk=new ReadGenbank(infos.get(i).genbank);
				for(int j=0;j<infos.size();j++){
					if(i==j)continue;
					String name2=infos.get(j).name;
					System.out.println(name1+"\t"+name2);
					File inMatches=new File(inFolder+"/"+name1+"/"+name2+".out");
					File out=new File(currentFolder+"/"+name2+".out");
					Matches matches=new Matches(inMatches);
					ArrayList<Info> deletions=getDeletions(matches,rgbk.getInfoTree("CDS"));
					Info.write(deletions, out);
					writeSequences(deletions,new File(currentFolder+"/"+name2+".fas"),infos.get(i).fasta);
				}
			}

	}
	
	public static void writeSequences(ArrayList<Info> deletions,File out,File genome){
		String Genome=Fasta.readFasta(genome).get(0).getSequence();
		ArrayList<Fasta> sequences=new ArrayList<Fasta>(); 
		for(int i=0;i<deletions.size();i++){
			int start=deletions.get(i).getStart();
			int end=deletions.get(i).getEnd();
			sequences.add(new Fasta("deletion_"+i,Genome.substring(start,end)));
		}
		Fasta.write(sequences, out);
	}
	
	public static ArrayList<Info> getDeletions(Matches matches,InfoTree gbk){
		ArrayList<Info> recs=new ArrayList<Info>();
		recs.addAll(matches.findDeletions());
		recs=addGeneInfo(recs,gbk);
		return recs;
	}
	
	public static ArrayList<Info> addGeneInfo(ArrayList<Info> recs,InfoTree gbkGenes){
		for(int i=0;i<recs.size();i++){
			ArrayList<Info> overlaps=new ArrayList<Info>();
			gbkGenes.search(recs.get(i), overlaps);
			recs.get(i).info+="_genes:";
			for(int j=0;j<overlaps.size();j++){
				recs.get(i).info+=overlaps.get(j).getGeneOrLocusInfo()+";";
				
			}
		}
		return recs;
	}
	

}
