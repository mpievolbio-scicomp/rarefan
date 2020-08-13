package sequenceDistribution;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import util.DNAmanipulations;
import util.ReadFasta;

public class GetSequence {
	HashMap<String,StringBuilder> fasta=new HashMap<String, StringBuilder>();
	public GetSequence(File genome){
		fasta=ReadFasta.readFasta(genome);
	}
	
	public String getSequence(String key,int start,int end){
		boolean rev=start>end;
		String sequence=rev?DNAmanipulations.reverse(fasta.get(key).substring(end,start)):fasta.get(key).substring(start,end);
		return sequence;
	}
	
	public void writeSequenceDistance(String SequenceKey[],Integer[] StartPositions,int numBasesDistance,File out){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			for(int i=0;i<StartPositions.length-1;i++){
				if(StartPositions[i+1]-StartPositions[i]==numBasesDistance){
					String Genome=fasta.get(">"+SequenceKey[i]).toString();
					bw.write(">"+(StartPositions[i]-14)+"-"+(StartPositions[i+1]+28)+"_"+SequenceKey[i]+"\n");
					bw.write(Genome.substring(StartPositions[i]-14,StartPositions[i+1]+28)+"\n");
				}
			}
			
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	public void write(String key,int start,int end,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write(">"+key+"_"+start+"_"+end+"\n");
			bw.write(getSequence(key,start,end)+"\n");
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
}
