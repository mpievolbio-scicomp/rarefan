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
	
	

	private static int getMaxEnd(ArrayList<Info> info) {
		int max=Integer.MIN_VALUE;
		for(int i=0;i<info.size();i++) {
			int curr=Math.max(info.get(i).getEnd(),info.get(i).getStart());
			max=curr>max?curr:max;
		}
		return max;
	}

	
	public static void writeGff(ArrayList<Info> info,File art,String seqId,ArrayList<Boolean> orientation,boolean deleteUnderscore,String name){ 
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(art));
			seqId=deleteUnderscore?deleteLastUnderscore(seqId):seqId;
			
			for(int i=0;i<info.size();i++){
				if(i==0) {
				   int endAnn=getMaxEnd(info);
				   bw.write("##gff-version 3.1.26\n");
				   bw.write("##sequence-region\t"+seqId+"\t"+1+"\t"+endAnn+"\n");
				}
				StringBuilder sb=new StringBuilder();
				String source="RAREFAN";
				String type="terminal_inverted_repeat_element";
				int start=info.get(i).start;
				int end=info.get(i).end;
				int score=0;
				String strand="+";
				if(start>end) {
					strand="-";
					int temp=start;
					start=end;
					end=temp;
				}
				if(orientation!=null) {
					if(orientation.get(i)==false) {
						strand="-";
					}
				}
				String phase="0";
				String attributes="Name="+name+"_"+i;
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
