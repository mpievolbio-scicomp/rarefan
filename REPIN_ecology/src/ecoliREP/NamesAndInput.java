package ecoliREP;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class NamesAndInput {
		public HashMap<String,Input> input;
		public ArrayList<String> names;
		public NamesAndInput(HashMap<String,Input> in,ArrayList<String> Names){
			input=in;
			names=Names;
		}
		public  NamesAndInput(File in){
			input=new HashMap<String, Input>();
			names=new ArrayList<String>();
			try{
				BufferedReader br=new BufferedReader(new FileReader(in));
				String line="";
				while((line=br.readLine())!=null){
					if(line.startsWith("#"))continue;
					String[] split=line.split("\\s+");
					File genome=new File(split[1]);
					File blastoutREP=null;
					String name=split[0];
					names.add(name);
					String homologue=split[4];
					File genbank=new File(split[2]);
					input.put(name,new Input(blastoutREP,genome,genbank,homologue,name));
				}
			}catch(IOException e){
				e.printStackTrace();
			}
			
			
		}
}
