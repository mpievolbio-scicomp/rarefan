package seedAnalysis;
//uses output of PrintWordFrequency
//class to plot the word length against different frequencies like max number of word occurrences and avg number
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
public class CreateWordPlot {
	public static void main(String args[]){
		File in=new File(args[0]);
		File out=new File(args[1]);
		getStats(in,new File(out+"/max.txt"),new File(out+"/avg.txt"));
		
	}
	
	private static void write(String text,File out){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			bw.write(text);
			bw.close();
			
		}catch(IOException e){
			System.err.println(e.toString());
		}
		
	}
	
	private static void getStats(File in,File max,File avg){
		HashMap<Integer,Integer> Max=new HashMap<Integer, Integer>();
		HashMap<Integer,Double> Avg=new HashMap<Integer, Double>();
		HashMap<Integer,String> MaxWord=new HashMap<Integer, String>();


		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			int sum=0;
			int count=0;
			int oldlength=0;
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				String key=split[0];
				int number=Integer.parseInt(split[1]);

				if(Max.containsKey(key.length())){
					if(Max.get(key.length())<number){
						Max.put(key.length(),number);
						MaxWord.put(key.length(), key);
					}
				}else{
					if(oldlength>0){
						Avg.put(oldlength, (sum*1.0)/count);
					}

					oldlength=key.length();
					sum=0;
					count=0;
					Max.put(key.length(),number);
					MaxWord.put(key.length(),key);
				}
				count++;
				sum+=number;
			}
			if(oldlength>0){
				Avg.put(oldlength, (sum*1.0)/count);
			}
		}catch(IOException e){
			System.err.println(e.toString());
		}
		write(printHashWord(Max,MaxWord),max);
		write(printHash(Avg),avg);
	}
	private static String printHashWord(HashMap<Integer,?> max,HashMap<Integer,String> word){
		StringBuilder plot=new StringBuilder();
		
		Iterator<?> it=max.entrySet().iterator(); 
		while(it.hasNext()){
			Entry<?,?> e =(Entry<?, ?>) it.next();
			plot.append(e.getKey()+"\t"+e.getValue()+"\t"+word.get(e.getKey())+"\n");
		}
		
		return plot.toString();
		
	}
	private static String printHash(HashMap<Integer,?> max){
		StringBuilder plot=new StringBuilder();
		
		Iterator<?> it=max.entrySet().iterator(); 
		while(it.hasNext()){
			Entry<?,?> e =(Entry<?, ?>) it.next();
			plot.append(e.getKey()+"\t"+e.getValue()+"\n");
		}
		
		return plot.toString();
		
	}
}
