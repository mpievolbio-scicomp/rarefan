package ecoliREP;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

import util.SortArrayList;

public class CompareMatrices {
	public static class MatEntry implements Comparable<MatEntry>{
		Double Value;
		String Name1;
		String Name2;
		public MatEntry(String name1,String name2,Double value){
			Name1=name1;
			Name2=name2;
			Value=value;
		}
		public String toString(){
			return Name1+"\t"+Name2+"\t"+Value;
		}
		@Override
		public int compareTo(MatEntry o) {
			// TODO Auto-generated method stub
			
			return Value.compareTo(o.Value);
		}
		
	}
	public static void main(String[] args){
		File mat1=new File(args[0]);
		File mat2=new File(args[1]);
		File frequency=new File(args[2]);
		File out=new File(args[3]);

		File outList=new File(args[4]);
		File outListMat1=new File(args[5]);
		File outListMat2=new File(args[6]);

		ArrayList<String> list=readList(mat1);
		HashMap<String,HashMap<String,Double>> ratio=getRatio(readMatrix(mat1, list),readMatrix(mat2, list),list);
		writeMatrix(out,ratio,list);
		
		writeList(SortArrayList.sort((makeList(ratio))),outList,getHomologue(frequency));//RandomizeList.randomise(getHomologue(frequency)));
		writeList(SortArrayList.sort(makeList(readMatrix(mat1, list))),outListMat1,getHomologue(frequency));
		
		writeList(SortArrayList.sort(makeList(readMatrix(mat2, list))),outListMat2,getHomologue(frequency));
	}
	
	public static HashMap<String,String> getHomologue(File homo){
		HashMap<String,String> hm=new HashMap<String, String>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(homo));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				hm.put(split[0],split[1]+"\t"+split[2]);
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		
		return hm;
	}
	
	public static void writeList(ArrayList<MatEntry> list,File out,HashMap<String,String> homologue){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<list.size();i++){
				bw.write(list.get(i).Value+"\t"+list.get(i).Name1+"\t"+homologue.get(list.get(i).Name1)+"\t"+list.get(i).Name2+"\t"+homologue.get(list.get(i).Name2)+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	// puts out a list of a matrix that gets sorted subsequently
	
	public static ArrayList<MatEntry> makeList(HashMap<String, HashMap<String,Double>> matrix){
		Iterator<Entry<String,HashMap<String,Double>>> it=matrix.entrySet().iterator();
		ArrayList<MatEntry> list=new ArrayList<MatEntry>();
		while(it.hasNext()){
			Entry<String,HashMap<String,Double>> e=it.next();
			String name1=e.getKey();
			Iterator<Entry<String,Double>> it2=e.getValue().entrySet().iterator();
			while(it2.hasNext()){
				Entry<String,Double> e2=it2.next();
				String name2=e2.getKey();
				Double val=e2.getValue();
				list.add(new MatEntry(name1,name2,val));
			}
		}
		return list;
	}
	
	public static void writeMatrix(File out,HashMap Matrix,ArrayList<String> list){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write("\t");
			for(int i=0;i<list.size();i++){
				bw.write(list.get(i)+"\t");
			}
			bw.write("\n");
			for(int i=0;i<list.size();i++){
				bw.write(list.get(i)+"\t");
				for(int j=0;j<list.size();j++){
					if(j==i){
						bw.write("0\t");
						continue;
					}
					if(!((HashMap<String,HashMap<String,Object>>)Matrix).get(list.get(i)).containsKey(list.get(j))){
						bw.write("0\t");
					}else{
						Object value=((HashMap<String,HashMap<String,Object>>)Matrix).get(list.get(i)).get(list.get(j));
						
						bw.write(value+"\t");
					}
				}
				bw.write("\n");
			}
			
			bw.close();
			
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	public static HashMap<String,HashMap<String,Double>> readMatrix(File matrixFile,ArrayList<String> list){
		HashMap<String,HashMap<String,Double>> matrix=new HashMap<String, HashMap<String,Double>>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(matrixFile));
			String line="";
			int i=0;
			while((line=br.readLine())!=null){
				i++;
				if(i==1)continue;
				String[] split=line.split("\\s+");
				for(int j=1;j<split.length;j++){
					
					if(matrix.containsKey(split[0])){
						matrix.get(split[0]).put(list.get(j-1), Double.parseDouble(split[j]));
					}else{
						HashMap<String,Double> hm=new HashMap<String, Double>();
						hm.put(list.get(j-1), Double.parseDouble(split[j]));
						matrix.put(split[0], hm);
					}
				}
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		
		return matrix;
	}
	
	public static ArrayList<String> readList(File Matrix){
		ArrayList<String> list=new ArrayList<String>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(Matrix));
			String line="";
			line=br.readLine();
			String[] split=line.split("\\s+");
			br.close();
			for(int i=1;i<split.length;i++){
				list.add(split[i]);
			}
		}catch(IOException e){
			e.printStackTrace();
		}
		
		return list;
		
	}
	
	public static HashMap<String,HashMap<String,Double>> getRatio(HashMap<String,HashMap<String,Double>> matrix1,HashMap<String,HashMap<String,Double>> matrix2,ArrayList<String> list){
		HashMap<String,HashMap<String,Double>> ratio=new HashMap<String, HashMap<String,Double>>();
		for(int i=0;i<list.size();i++){
			String name1=list.get(i);
			for(int j=0;j<list.size();j++){
				String name2=list.get(j);
				double Ratio=matrix2.get(name1).get(name2)!=0?matrix1.get(name1).get(name2)/matrix2.get(name1).get(name2):Double.NaN;
				if(ratio.containsKey(name1)){
					ratio.get(name1).put(name2,Ratio);
				}else{
					HashMap<String,Double> hm=new HashMap<String, Double>();
					hm.put(name2,Ratio);
					ratio.put(name1, hm);
				}
				
			}
		}
		
		
		return ratio;
	}
	
	public static double getMin(HashMap<String,HashMap<String,Double>> matrix,ArrayList<String> list){
		double min=Double.MAX_VALUE;
		for(int i=0;i<list.size();i++){
			for(int j=0;j<list.size();j++){
				if(i==j)continue;
				if(matrix.get(list.get(i)).get(list.get(j))!=Double.NaN && matrix.get(list.get(i)).get(list.get(j))<min){
					min=matrix.get(list.get(i)).get(list.get(j));
				}
			}
		}
		return min;
	}
	
	public static double getMax(HashMap<String,HashMap<String,Double>> matrix,ArrayList<String> list){
		double max=Double.MIN_VALUE;
		for(int i=0;i<list.size();i++){
			for(int j=0;j<list.size();j++){
				if(i==j)continue;
				if(matrix.get(list.get(i)).get(list.get(j))!=Double.NaN &&matrix.get(list.get(i)).get(list.get(j))>max){
					max=matrix.get(list.get(i)).get(list.get(j));
				}
			}
		}
		return max;
	}
	public static double getAverage(HashMap<String,HashMap<String,Double>> matrix,ArrayList<String> list){
		double avg=0;
		for(int i=0;i<list.size();i++){
			for(int j=0;j<list.size();j++){
				if(i==j)continue;
				if(matrix.get(list.get(i)).get(list.get(j))!=Double.NaN)avg+=matrix.get(list.get(i)).get(list.get(j));
			}
		}
		return avg/list.size();
	}
}
