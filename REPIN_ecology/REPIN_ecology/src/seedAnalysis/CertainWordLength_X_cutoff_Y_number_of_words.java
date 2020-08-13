package seedAnalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class CertainWordLength_X_cutoff_Y_number_of_words {
	public static void main(String[] args){
		File words=new File(args[0]);
		int wordlength=Integer.parseInt(args[1]);
		File out=new File(args[2]);
		write(createPlot(words,wordlength),out);
	}
	
	public static void write(ArrayList<Double> plot,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));

			for(int i=2;i<plot.size();i++){
				bw.write(i+"\t"+plot.get(i)+"\n");
			}
			bw.close();	
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
	private static ArrayList<Double> createPlot(File in,int wl){
		ArrayList<Double> plot=new ArrayList<Double>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				if(wl<split[0].length())break;
				if(wl>split[0].length())continue;
				int freq=Integer.parseInt(split[1]);
				for(int i=2;i<=freq;i++){
					if(plot.size()>i){
						plot.set(i,plot.get(i)+1);
					}else{
						plot.add(1.0);
					}
				}
			}
			
		}catch(IOException e){
			System.err.println(e.toString());
		}
		
		
		return plot;
	}
	
}
