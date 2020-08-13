package frequencies;

import java.io.*;
import java.util.*;

public class CalculateConstantSites {
	public static void main(String args[]){
		File inFolder=new File(args[0]);
		File out=new File(inFolder+"/numConstantSites.txt");
		CalculateConstantSites ccs=new CalculateConstantSites(inFolder, 25);
		ccs.printNumConstantSites(out);
	}
	File inFolder;
	HashMap<String,HashSet<Integer>> constantSites=new HashMap<String,HashSet<Integer>>();
	
	public CalculateConstantSites(File inFolder,int seqLength){
		this.inFolder=inFolder;
		ArrayList<File> files=getFiles(inFolder);
		initHash(files,seqLength);
		readMutFreqs(files);
	}
	
	
	private ArrayList<File> getFiles(File inFolder){
		File[] all=inFolder.listFiles();
		ArrayList<File> files=new ArrayList<File>();
		for(int i=0;i<all.length;i++){
			if(all[i].getName().endsWith(".mf")){
				files.add(all[i]);
			}
		}
		return files;
	}
	
	private int getPos(String mut){
		String temp=mut.substring(1, mut.length()-1);
		return Integer.parseInt(temp);
	}
	
	private void readMutFreqs(ArrayList<File> files){
		for(int i=0;i<files.size();i++){
			readMutFreqs(files.get(i));
		}
	}
	
	private void readMutFreqs(File in){
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			String name=in.getName().split("\\.")[0];
			while((line=br.readLine())!=null){
				String mut=line.split("\\s+")[0];
				int pos=getPos(mut);
				constantSites.get(name).remove(pos);
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private void initHash(ArrayList<File> files,int length){
		for(int i=0;i<files.size();i++){
			String name=files.get(i).getName().split("\\.")[0];
			initHash(name,length);
		}
	}
	
	private void initHash(String key,int length){
		HashSet<Integer> constantSitesTemp=new HashSet<Integer>();
		for(int i=0;i<length;i++){
			constantSitesTemp.add(i);
		}
		constantSites.put(key, constantSitesTemp);
	}
	
	public void printNumConstantSites(File out){
		String[] keys=constantSites.keySet().toArray(new String[0]);
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<keys.length;i++){
				bw.write(keys[i]+"\t"+constantSites.get(keys[i]).size()+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	
}
