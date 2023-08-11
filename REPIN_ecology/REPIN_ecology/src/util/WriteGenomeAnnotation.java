package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class WriteGenomeAnnotation {
	public static void writeTab(ArrayList<Info> info,File art){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(art));
			for(int i=0;i<info.size();i++){
				String ide=info.get(i).info.length()>5?info.get(i).info.substring(0,5):info.get(i).info;
				StringBuilder sb=new StringBuilder(ide);
				for(int j=0;j<5-ide.length();j++){
					sb.append(" ");
				}
				ide=sb.toString();
				int start=info.get(i).getStart();
				int end=info.get(i).getEnd();
				String inf=info.get(i).info;
				if(start<end)
					bw.write("FT   "+ide+"             "+start+".."+end+"\n");
				else bw.write("FT   "+ide+"             "+"complement("+start+".."+end+")"+"\n");
				bw.write("FT                   /note=\""+inf+"\"\n");
			}
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
	private static String deleteLastUnderscore(String seqId) {
		String[] split=seqId.split("_");
		StringBuffer newseq=new StringBuffer();
		newseq.append(split[0]);
		for(int i=1;i<split.length-1;i++) {
			newseq.append("_"+split[i]);
		}
		return newseq.toString();
	}
	
	public static void writeGff(ArrayList<Info> info,File art,String seqId,ArrayList<Boolean> orientation){ 
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(art));
			seqId=deleteLastUnderscore(seqId);
			for(int i=0;i<info.size();i++){
				StringBuilder sb=new StringBuilder();
				String source="RAREFAN";
				String type="terminal_inverted_repeat_element";
				int start=info.get(i).start;
				int end=info.get(i).end;
				int score=0;
				String strand="+";
				if(orientation!=null) {
					if(orientation.get(i)==false) {
						strand="-";
					}
				}
				String phase="0";
				String attributes="none";
				sb.append(seqId+"\t");
				sb.append(source+"\t");
				sb.append(type+"\t");
				sb.append(start+"\t"+end+"\t");
				sb.append(score+"\t");
				sb.append(strand+"\t");
				sb.append(phase+"\t");
				sb.append(attributes+"\n");
				bw.write(sb.toString());
			}
			bw.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
	public static void writeTab(ArrayList<Integer> start,ArrayList<Integer> end,ArrayList<String> id,File art){
			ArrayList<Info> info=new ArrayList<Info>();
			for(int i=0;i<start.size();i++){
				info.add(new Info(start.get(i),end.get(i),id.get(i)));
			}
			writeTab(info,art);
		
	}
	
	
	

}
