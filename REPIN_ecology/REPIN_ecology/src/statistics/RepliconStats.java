package statistics;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import statistics.Stats;

public class RepliconStats {
	public static void main(String[] args){
		File inFolder=new File(args[0]);
		String fileName=args[1];
		try{
			File[] folders=inFolder.listFiles();
			BufferedWriter bw=new BufferedWriter(new FileWriter(new File(inFolder+"/repliconStats.txt")));
			for(int i=0;i<folders.length;i++){
				if(folders[i].isDirectory()){
					File in=new File(folders[i]+"/"+fileName);
					if(in.exists()){
						ArrayList<Double> list=readRepliconData(in);
						Stats stats=new Stats(list);
						String split[]=folders[i].getName().split("_");
						String name=folders[i].getName();
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
	public static ArrayList<Double> readRepliconData(File in){	
		HashMap<String,Integer> replicons=readReplicons(in);
		//bw.write("Total: "+hm.entrySet().size()+"\n");
		Iterator<Entry<String,Integer>> it=replicons.entrySet().iterator();
		ArrayList<Double> list=new ArrayList<Double>();
		while(it.hasNext()){
			Entry<String,Integer> e=it.next();
			list.add(0.0+e.getValue());
			//System.out.println(e.getValue());
		}
		return list;
	}
	public static HashMap<String,Integer> readReplicons(File in){
		HashMap<String,Integer> reps=new HashMap<String, Integer>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				if(line.startsWith(">")){
					String[] split=line.split("_");
					//System.out.println(split[2]);
					if(reps.containsKey(split[2])){
						reps.put(split[2],reps.get(split[2])+1);
					}else{
						reps.put(split[2], 1);
					}
				}
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return reps;
	}
}
