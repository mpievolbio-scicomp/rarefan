package overrepresented;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import util.DNAmanipulations;
import util.Fasta;
import util.Find;
import util.Histogram;
import util.Info;
import util.ReadFasta;
import util.WriteGenomeAnnotation;

public class GroupSeedSequences {
	public static void main(String args[]){
		File genome=new File(args[0]);
		File searchStringFile=new File(args[1]);
		File outputFolder=new File(args[2]);
		int flanking=Integer.parseInt(args[3]);
		groupSeedSequences(genome,searchStringFile,outputFolder,flanking);

		
	}
	
	public static String[] groupSeedSequences(File genome,File searchStringFile,File outputFolder,int flanking) {
		HashMap<String,StringBuilder> genomeFasta=ReadFasta.readFasta(genome);
		//readFile, needs input from SelectOverrepresentedWords.java
		//pull out most abundant+flanking sequence
		//concatenated the sequences separated by |
		//write everything in a file which doesnt match+frequency
		//write everything in a file which does match+frequency
		//iterate process
		int i=0;
		ArrayList<String> seedSeqs=new ArrayList<String>();
		String id=searchStringFile.getName().split("\\.")[0];
		ArrayList<Info> wordPositions=new ArrayList<Info>();
		while(searchStringFile.length()>1){

			String word=getWord(searchStringFile);
			System.out.println("Most common sequence in group "+i+" in the reference genome: "+word);
			wordPositions.addAll(getPositions(Fasta.readFasta(genome).get(0).getSequence(),word,"Group_"+i));
			seedSeqs.add(word);
			Iterator<Entry<String,StringBuilder>> it=genomeFasta.entrySet().iterator();
			StringBuilder wordSequence=new StringBuilder();
			while(it.hasNext()){
				Entry<String,StringBuilder> e=it.next();
				wordSequence.append(getSequences(word,e,flanking));

			}
			//System.out.println(wordSequence);
			searchStringFile=writeFiles(searchStringFile,wordSequence.toString(),outputFolder,i,id);
			i++;
		}
		File artemisOut=new File(outputFolder+"/"+id+"_words.tab");
		File gffOut=new File(outputFolder+"/"+id+"_words.gff3");
		WriteGenomeAnnotation.writeTab(wordPositions, artemisOut);
		WriteGenomeAnnotation.writeGff(wordPositions,gffOut,id,null);

		return seedSeqs.toArray(new String[0]);
	}
	
	private static ArrayList<Info> getPositions(String genome,String search,String id) {
		ArrayList<Info> posInfo=new ArrayList<Info>();
		ArrayList<Integer> pos=Find.getPositions(search, genome);
		pos.addAll(Find.getPositions(DNAmanipulations.reverse(search), genome));
		for(int i=0;i<pos.size();i++) {
			posInfo.add(new Info(pos.get(i),pos.get(i)+search.length(),id));
		}
		return posInfo;
	}
	
	private static File writeFiles(File searchStringFile,String wordSequence,File path,int group,String id){
		File newSearchStringFile=new File(path+"/"+"GroupNotFound"+group+".out");
		newSearchStringFile.deleteOnExit();

		try{
			File found=new File(path+"/"+"Group_"+id+"_"+group+".out");
			BufferedReader br=new BufferedReader(new FileReader(searchStringFile));
			BufferedWriter bwNotFound=new BufferedWriter(new FileWriter(newSearchStringFile));
			HashMap<String,Integer> foundHash=new HashMap<String, Integer>();
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				String word=split[0];
				String revWord=DNAmanipulations.reverse(word);
				int freq=Integer.parseInt(split[1]);
				if(wordSequence.contains(word) || wordSequence.contains(revWord) ){
					if(foundHash.containsKey(revWord)){
						foundHash.put(revWord, foundHash.get(revWord)+freq);
					}else{
						foundHash.put(word, freq);
					}
				}else{

					bwNotFound.write(line+"\n");
				}
			}
			br.close();
			bwNotFound.close();
			Histogram.write(foundHash,found);
			writeFasta(foundHash,new File(found+".fas"));
		}catch(IOException e){
			System.err.println(e.toString());
		}
		
		return newSearchStringFile;
	}
	
	private static  void writeFasta(HashMap<String,Integer> map,File out) {
		try {
			String[] keys=map.keySet().toArray(new String[0]);
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<keys.length;i++) {
				bw.write(">seq"+i+" "+map.get(keys[i])+"\n"+keys[i]+"\n");
			}
			bw.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private static String getWord(File in){
		String word="";
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			int max=0;
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				int freq=Integer.parseInt(split[1]);
				if(max<freq) {
					word=split[0];
					max=freq;
				}
			}
			br.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
		return word;
	}
	
	
	private static String getSequences(String word,Entry<String,StringBuilder> e,int flanking){
		StringBuilder sb=buildString(word,e,flanking,true);
		sb.append(buildString(word,e,flanking,false));
		return sb.toString();
		
	}
	
	private static StringBuilder buildString(String word,Entry<String,StringBuilder> e,int flanking,boolean reverse){
		int i=0;
		StringBuilder result=new StringBuilder();
		word=word.toUpperCase();
		String sequence=reverse?DNAmanipulations.reverse(e.getValue().toString().toUpperCase()):e.getValue().toString().toUpperCase();
		while((i=sequence.indexOf(word, i))!=-1){
			String value=sequence.length()>i+word.length()+flanking && i-flanking>0?sequence.substring(i-flanking,i+word.length()+flanking):sequence.substring(i,i+word.length());
			result.append("|"+value);
			i++;
		}
		return result;
	}
}
