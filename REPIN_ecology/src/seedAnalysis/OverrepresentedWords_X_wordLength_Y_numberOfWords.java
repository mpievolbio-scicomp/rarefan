package seedAnalysis;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
public class OverrepresentedWords_X_wordLength_Y_numberOfWords {
	public static void main(String[] args){
		File words=new File(args[0]);
		//avg file
		File overrepresented=new File(args[1]);
		ArrayList<Double> over=readOverrepresented(overrepresented);
		File out=new File(args[2]);
		CertainWordLength_X_cutoff_Y_number_of_words.write(createPlot(words,over),out);
	}
	
	private static ArrayList<Double> readOverrepresented(File avg){
		ArrayList<Double> al=new ArrayList<Double>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(avg));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				int wl=Integer.parseInt(split[0]);
				double Avg=Double.parseDouble(split[1]);
				if(wl>=al.size()){
					for(int i=0;i<=wl;i++){
						al.add(0.0);
					}
					al.set(wl,Avg);
				}else{
					al.set(wl,Avg);
				}
			}
		}catch(IOException e){
			System.err.println(e.toString());
		}
		return al;
	}
	
	private static ArrayList<Double> createPlot(File in,ArrayList<Double> over){
		ArrayList<Double> plot=new ArrayList<Double>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				int freq=Integer.parseInt(split[1]);
				int wl=split[0].length();
				if(freq>over.get(wl)){
					if(plot.size()>wl){
						plot.set(wl,plot.get(wl)+1);
					}else{
						for(int i=plot.size();i<=wl;i++)
							plot.add(0.0);
						plot.set(wl,1+0.0);
					}
				}
			}

			
		}catch(IOException e){
			System.err.println(e.toString());
		}
		
		
		return plot;
	}
	
}
