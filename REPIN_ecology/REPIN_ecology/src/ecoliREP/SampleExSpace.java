package ecoliREP;

import java.io.*;
import java.util.*;
import ecoliREP.CompareMatrices.MatEntry;
import util.SortArrayList;


//this program is supposed to sample a certain number of extragenic spaces repeatedly(100 times) and determine
//the number of invasion events, the output is a list of organisms that were compared to each other and for each pair
//a certain number of (100) invasion event values
//the inputs should be: 
//1. a file that contains the all paths to the extragenic space data (ExtraComparison.in)
//2. a path to the match files in order to find the number of extragenic spaces that are to be sampled for each comparison (REP match folder: /home/frederic/auckland/outputFiles/REPAnalysisEcoli/matches/Genomes/)
//the number of repetitions (100)
//output file
public class SampleExSpace {
	ArrayList<ExtraSpaces> spaces=new ArrayList<ExtraSpaces>();
	ArrayList<PairList> results=new ArrayList<PairList>();
	ArrayList<MatEntry> comparison=new ArrayList<MatEntry>();
	HashMap<String,HashMap<String,Double>> comparisonMatrix=new HashMap<String,HashMap<String,Double>>();
	ArrayList<String> names=new ArrayList<String>();
	File out;
	
	public ArrayList<ExtraSpaces> getSpaces(){
		return spaces;
	}
	
	public static void main(String[] args){
		File master=new File(args[0]);
		File matchFolderReps=new File(args[1]);
		int repetitions=Integer.parseInt(args[2]);
		File matrixFile=new File(args[3]);
		File out=new File(args[4]);
		SampleExSpace ses=new SampleExSpace(master,matchFolderReps,repetitions,out); 
		ArrayList<String> list=CompareMatrices.readList(matrixFile);
		ses.compare(CompareMatrices.readMatrix(matrixFile, list));
		ses.writeCompare();
		ses.writeCompareMatrix(list);
	}
	
	public double getPosition(ArrayList<Integer> list,double value){
		for(int i=0;i<list.size();i++){
			if(list.get(i)>value){
				return i;
			}
		}
		return list.size();
	}
	
	public void compare(HashMap<String,HashMap<String,Double>> REPs){
		for(int i=0;i<results.size();i++){
			String name1=results.get(i).strain1;
			String name2=results.get(i).strain2;
			ArrayList<Integer> list=results.get(i).list;
			double REPInvasion=REPs.get(name1).get(name2);
			double position=getPosition(list,REPInvasion);
			comparison.add(new MatEntry(name1,name2,position));
			if(comparisonMatrix.containsKey(name1)){
				comparisonMatrix.get(name1).put(name2, position);
				
			}else{
				HashMap<String,Double> temp=new HashMap<String, Double>();
				temp.put(name2,position);
				comparisonMatrix.put(name1, temp);
			}
		}
	}
	
	public void writeCompare(){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out+"/comparison.out"));
			for(int i=0;i<comparison.size();i++){
				bw.write(comparison.get(i)+"\n");
			}
			
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	public void writeCompareMatrix(ArrayList<String> list){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out+"/comparisonMatrix.out"));
			for(int i=0;i<list.size();i++){
				bw.write("\t"+list.get(i));
			}
			bw.write("\n");
			for(int i=0;i<list.size();i++){
				bw.write(list.get(i)+"\t");
				for(int j=0;j<list.size();j++){
					if(comparisonMatrix.containsKey(list.get(i))){
						if(comparisonMatrix.get(list.get(i)).containsKey(list.get(j))){
							bw.write(comparisonMatrix.get(list.get(i)).get(list.get(j))+"\t");
						}else{
							bw.write("NaN\t");
						}
					}else{
						bw.write("NaN\t");
					}
				}
				bw.write("\n");
			}
			
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	public SampleExSpace(){
		
	}
	
	public SampleExSpace(File master,File matchFolderREPs,int repetitions,File Out){
		out=Out;
		ArrayList<Infos> infos=Infos.readInfos(master);
		for(int i=0;i<infos.size();i++){
			System.out.println(infos.get(i).name);
			for(int j=0;j<infos.size();j++){
				if (j==i)continue;
				String strain1=infos.get(i).name;
				String strain2=infos.get(j).name;
				File matchFile=new File(infos.get(i).matchFolder+"/"+strain2+".out");
				loadSpaces(matchFile);
				int number=getFrequencies(new File(matchFolderREPs+"/"+strain1+"/"+strain2+".out"));
				ArrayList<Integer> list=new ArrayList<Integer>();
				for(int k=0;k<repetitions;k++){
					ArrayList<ExtraSpaces> sam=sample(number);
					list.add(getInvasions(sam));
				}
				list=SortArrayList.sort(list);
				results.add(new PairList(list,strain1,strain2));
				
			}
		}
		write();
	}
	
	public int getFrequencies(File in){
		int i=0;
		
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			while(br.readLine()!=null){
				i++;
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return i/2;
	}
	
	public void write(){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out+"/samples.out"));
			for(int i=0;i<results.size();i++){
				bw.write(results.get(i).strain1+"\t"+results.get(i).strain2+"\t");
				for(int j=0;j<results.get(i).list.size();j++){
					bw.write(results.get(i).list.get(j)+"\t");
				}
				bw.write("\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	
	public void loadSpaces(File in){
		try{
			spaces=new ArrayList<ExtraSpaces>();
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			int i=0;
			String l="";
			while((line=br.readLine())!=null){
				if(line.startsWith("REPleft")){
					l=line;
				}else if(line.startsWith("REPright")&&i>0){
					spaces.add(new ExtraSpaces(l,line));
				}
				i++;
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public ArrayList<ExtraSpaces> sample(int number){
		ArrayList<ExtraSpaces> sample=new ArrayList<ExtraSpaces>();
		HashMap<Integer,Boolean> exists=new HashMap<Integer, Boolean>();
		exists.put(-1, true);
		while(sample.size()<number && number<spaces.size()){
			int rand=-1;
			while(exists.containsKey(rand)){
				rand=(int)(Math.random()*spaces.size());
			}
			sample.add(spaces.get(rand));
			exists.put(rand, true);		
		}
		return sample;
	}
	
	public int getInvasions(ArrayList<ExtraSpaces> sample){
		int invasions=0;
		for(int i=0;i<sample.size();i++){
			if(isInvasion(sample.get(i))){
				invasions++;
			}
		}		
		return invasions;
	}
	
	public boolean isInvasion(ExtraSpaces sample){
		boolean a=sample.left.contains("missing");
		boolean b=sample.right.contains("missing");
		return (a||b)&&!(a&&b);
	}
	

	
	class PairList{
		ArrayList<Integer> list=new ArrayList<Integer>();
		String strain1="";
		String strain2="";
		public PairList(ArrayList<Integer> l,String s1,String s2){
			list=l;
			strain1=s1;
			strain2=s2;
		}
	}
	
	public static class ExtraSpaces{
		public ExtraSpaces(String l,String r){
			left=l;
			right=r;
		}
		String left;
		String right;
	}
	
}
