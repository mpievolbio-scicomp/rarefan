package seedAnalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map.Entry;

import util.WordFrequency;

public class OverrepresentedWordPlot {
	public static void main(String[] args){
		File words=new File(args[0]);
		File out=new File(args[1]);
		int stop=Integer.parseInt(args[2]);
		TreeMap<WordFrequency,Boolean> tm=sort(words);
		print(tm,out,stop);
	}
	public static void print(TreeMap<WordFrequency,Boolean> tm,File out,int stop){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			
			Iterator<Entry<WordFrequency,Boolean>> it=tm.entrySet().iterator();
			int i=0;
			while(it.hasNext()){
				Entry<WordFrequency,Boolean> e=it.next();
				bw.write(i+"\t"+e.getKey()+"\n");
				i++;
				if(stop<i){
					break;
				}
			}
			
			
			bw.close();
		}catch(IOException e){
			System.err.println(e);
		}
	}
	private static TreeMap<WordFrequency,Boolean> sort(File words){
		TreeMap<WordFrequency, Boolean> tm=new TreeMap<WordFrequency, Boolean>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(words));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				int freq=Integer.parseInt(split[1]);
				String word=split[0];
				tm.put(new WordFrequency(freq,word), true);
			}
			
		}catch(IOException e){
			System.err.println(e.toString());
		}
		return tm;
	}
	
}
