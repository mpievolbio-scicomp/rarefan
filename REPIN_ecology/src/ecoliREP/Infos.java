package ecoliREP;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;



public  class Infos{
	public Infos(String Name,String Fasta,String Genbank,String MatchFolder){
		name=Name;
		fasta=new File(Fasta);
		genbank=new File(Genbank);
		matchFolder=new File(MatchFolder);
	}
	String name;
	File fasta;
	File genbank;
	File matchFolder;
	public static ArrayList<Infos> readInfos(File infos){
		ArrayList<Infos> al=new ArrayList<Infos>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(infos));
			String line="";
			while((line=br.readLine())!=null){
				if(line.startsWith("#"))continue;
				String[] split=line.split("\\s+");
				Infos info=new Infos(split[0],split[1],split[2],split[3]);
				al.add(info);
			}
			
		}catch(IOException e){
			e.printStackTrace();
		}
		return al;
	}
}

