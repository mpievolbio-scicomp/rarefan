package ecoliREP;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;


import ecoliREP.SampleExSpace.ExtraSpaces;

public class ConstantREPs {
	public static void main(String args[]){
		File infoFile=new File(args[0]);
		File folder=new File(args[1]);
		File out=new File(args[2]);
		ArrayList<Infos> infos=Infos.readInfos(infoFile);
		HashMap<String,Integer> constants=new HashMap<String, Integer>();
		for(int i=1;i<infos.size();i++){
			readMatchFile(new File(folder+"/"+infos.get(0).name+"/"+infos.get(i).name+".out"),constants);
		}
		printConstants(out,constants);
	}
	
	public static void printConstants(File out,HashMap<String,Integer> constants){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			Iterator<Entry<String, Integer>> it=constants.entrySet().iterator();
			while(it.hasNext()){
				Entry<String,Integer> e=it.next();
				if(e.getValue()>=35)
				bw.write(e.getKey()+"\t"+e.getValue()+"\n");
			}
			bw.close();
			
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static void readMatchFile(File in,HashMap<String,Integer> constants){
		ArrayList<ExtraSpaces> al=new ArrayList<ExtraSpaces>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			String left="";
			while((line=br.readLine())!=null){
				String split[]=line.split("\\s+");
				if(split[0].contains("left")){
					left=line;
				}else{
					al.add(new ExtraSpaces(left,line));
				}
			}
			for(int i=0;i<al.size();i++){
				if(checkMatch(al.get(i))){
					left=al.get(i).left.split("\\s+")[0];
					String right=al.get(i).right.split("\\s+")[0];
					String key=left+"\t"+right;
					if(!constants.containsKey(key)){
						constants.put(key,1);
					}else{
						constants.put(key, constants.get(key)+1);
					}
				}
					
			}
			
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static boolean checkMatch(ExtraSpaces exSpace){
		String left=exSpace.left.split("\\s+")[1];
		String right=exSpace.right.split("\\s+")[1];
		if(left.contains("missing")||right.contains("missing")){
			return false;
		}else{
			if((left.contains("left") && right.contains("left"))||(left.contains("right") && right.contains("right"))){
				return false;
			}
		}
		return true;
	}
	
}	
