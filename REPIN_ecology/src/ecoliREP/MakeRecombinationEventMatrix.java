package ecoliREP;

import java.io.*;
import java.util.*;

public class MakeRecombinationEventMatrix {
	public static void main(String args[]){
		File in=new File(args[0]);
		File inFolder=new File(args[1]);
		ArrayList<Infos> config=Infos.readInfos(in);
		HashMap<String,HashMap<String,Integer>> hm=new HashMap<String, HashMap<String,Integer>>();
		ArrayList<String> list=new ArrayList<String>();
		for(int i=0;i<config.size();i++){
			list.add(config.get(i).name);
			String name1=config.get(i).name;
			for(int j=0;j<config.size();j++){
				if(i==j)continue;
				String name2=config.get(j).name;
				File current=new File(inFolder+"/"+name1+"/"+name2+".out");
				int recs=readFile(current);
				if(hm.containsKey(name1)){
					hm.get(name1).put(name2,recs);
				}else{
					HashMap<String,Integer> temp=new HashMap<String, Integer>();
					temp.put(name2,recs);
					hm.put(name1,temp);
				}
			}
		}
		CompareMatrices.writeMatrix(new File(inFolder+"/recMAtrix.out"), hm, list);
	}
	private static int readFile(File in){
		int lines=0;
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			while(br.readLine()!=null){
				lines++;
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return lines;
	}
}
