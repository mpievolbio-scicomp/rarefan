package statistics;

import java.io.*;
import java.util.*;

import statistics.Stats;

public class AnalyseWordFreqFlank {
	public static void main(String args[]){
		
		File wf=new File(args[0]);
		
		File[] folders=wf.listFiles();
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(wf+"/wfrStats.txt"));
			for(int i=0;i<folders.length;i++){
				if(folders[i].isDirectory()){

					File in=new File(folders[i]+"/maxWordFreq.txt");
					if(in.exists()){
						ArrayList<Double> freq=readWF(in);
						Stats stats=new Stats(freq);
						String name=folders[i].getName();
						String split[]=folders[i].getName().split("_");
						if(split.length>1){
							name=split[1];
						}
						bw.write(name+"\t"+stats.getAverage()+"\t"+stats.getStandardDeviation()+"\t"+stats.getStandardError()+"\n");
					}
				}
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	public static ArrayList<Double> readWF(File wf){
		ArrayList<Double> freq=new ArrayList<Double>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(wf));
			String line="";
			while((line=br.readLine())!=null){
				if(line.matches("\\D+.+")){
					String[] split=line.split("\\s+");
					double f=Double.parseDouble(split[1]);
					if(f>0)freq.add(f);
				}
			}
			
		}catch(IOException e){
			e.printStackTrace();
		}
		
		return freq;
	}
}
