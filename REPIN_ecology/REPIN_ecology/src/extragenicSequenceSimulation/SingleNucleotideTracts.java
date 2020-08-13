package extragenicSequenceSimulation;

import java.io.File;

public class SingleNucleotideTracts {
	public static void main(String[] args){
		File in=new File(args[0]);
		File outFolder=new File(in.getParent()+"/out/");
		outFolder.mkdir();
		char[] chars=new char[]{'A','C'};
		for(int j=0;j<chars.length;j++){
			for(int i=5;i<11;i++){
				File out=new File(outFolder+"/tract_"+i+"_"+chars[j]+".txt");
				File outFas=new File(out+".fas");
				String word=getWord(i,chars[j]);
				MutatedSequenceOccurrences.main(new String[]{in+"","null",0+"",out+"",1+"",outFas+"",0+"","full",word});
			}
		}

	}
	private static String getWord(int length,char c){
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<length;i++){
			sb.append(c);
		}
		return sb.toString();
	}
}
