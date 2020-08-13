package overrepresented;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class SelectOverRepresentedWords {
	public static void main(String args[]){
		File wf=new File(args[0]);
		int cutOff=Integer.parseInt(args[1]);
		File out=new File(args[2]);
		int wordLength=Integer.parseInt(args[3]);
		readAndWrite(wf, out, cutOff,wordLength);
		
	}
	
	public static void readAndWrite(File in,File out,int cutOff,int wordLength){
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				
				int frequency=Integer.parseInt(split[1]);
				if(split[0].length()==wordLength && frequency>=cutOff){
					bw.write(line+"\n");
				}
				if(split[0].length()>wordLength){
					break;
				}
			}
			bw.close();
			br.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
}
