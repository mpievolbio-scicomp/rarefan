package ecoliREP;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

public class REPdifferenceHistogram {
	public static void main(String args[]){
		File in=new File(args[0]);
		int limit=Integer.parseInt(args[1]);
		File out=new File(args[2]);
		ArrayList<ResList> resList=readResList(in);
		ArrayList<Histogram> hist=createHistogram(resList);
		write(hist,out,limit);
	}
	
	public static ArrayList<Histogram> createHistogram(ArrayList<ResList> resList){
		ArrayList<Histogram> al=new ArrayList<Histogram>();
		String name="";
	
		for(int i=0;i<resList.size();i++){
			if(resList.get(i).s1.equals(name)){
				al.get(al.size()-1).put(resList.get(i).repNum);
			}else{
				Histogram hist=new Histogram(resList.get(i).s1);
				hist.put(resList.get(i).repNum);
				al.add(hist);
				name=resList.get(i).s1;
			}
		}
		return al;
	}
	
	public static void write(ArrayList<Histogram> al,File out,int limit){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<al.size();i++){
				Iterator<Entry<Integer,Integer>> it=al.get(i).frequency.entrySet().iterator();
				while(it.hasNext()){
					Entry<Integer,Integer> e=it.next();
					if(e.getValue()<=limit){
						bw.write(al.get(i).name+"\t"+e.getKey()+"\t"+e.getValue()+"\n");
					}
				}
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static ArrayList<ResList> readResList(File in){
		ArrayList<ResList> al=new ArrayList<ResList>();
			try{
				BufferedReader br=new BufferedReader(new FileReader(in));
				String line="";
				int i=0;
				while((line=br.readLine())!=null){
					if(i==0){
						i++;
						continue;
					}
					String[] split=line.split("\t");
					int start=Integer.parseInt(split[3]);
					int end=Integer.parseInt(split[4]);
					String annotation=split[5];
					al.add(new ResList(split[0],split[1],Integer.parseInt(split[2]),start,end,annotation));
				}
				br.close();
			}catch(IOException e){
				e.printStackTrace();
			}
		return al;
	}
	
	public static class Histogram{
		String name;
		TreeMap<Integer,Integer> frequency;
		public Histogram(String Name){
			name=Name;
			frequency=new TreeMap<Integer,Integer>();
		}
		public void put(int number){
			if(frequency.containsKey(number)){
				frequency.put(number,frequency.get(number) +1);
			}else{
				frequency.put(number, 1);
			}
		}
	}
	
	public static class ResList{
		String s1;
		String s2;
		int repNum;
		int start;
		int end;
		String annotation;
		public ResList(String strain1,String strain2,int RepNum,int Start,int End,String Annotation){
			s1=strain1;
			s2=strain2;
			repNum=RepNum;
			start=Start;
			end=End;
			annotation=Annotation;
		}
		
	}
	
}
