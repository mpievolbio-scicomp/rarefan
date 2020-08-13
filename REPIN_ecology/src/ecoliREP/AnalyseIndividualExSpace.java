package ecoliREP;

import java.io.*;
import java.util.*;


import util.Info;


//given a strain name, a specific extragenic space number in the particular strain, 
//a folder that contains all exspace fasta files, a file containing all strain names that are to be compared with
//the particular strain (REPComparison.in)

//number of REPs in Space


public class AnalyseIndividualExSpace {
	public static void main(String args[]){
		File REPseq=new File(args[0]);
		File input=new File(args[1]);//like REPComparison.in
		File matchFolder=new File(args[2]);
		File outFolder=new File(args[3]);
		BlastREPs bReps=new BlastREPs(outFolder,true);
		bReps.createBlast(input,"blastn", 0.001, REPseq,false,false);
		ArrayList<String> names=bReps.getNames();
		HashMap<String,String> totalSpaceMissing=new HashMap<String, String>();
		for(int i=0;i<names.size();i++){
			System.out.println(names.get(i));
			
			totalSpaceMissing.put(names.get(i), countREPs(matchFolder,names.get(i),bReps.getNames(),bReps,outFolder));
		}
		write(totalSpaceMissing,outFolder,names,bReps);
	}
	
	public static void write(HashMap<String,String> tsm,File out,ArrayList<String> names,BlastREPs brep){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out+"/spacesMissing.out"));
			for(int i=0;i<names.size();i++){
				String ident=names.get(i)+"\t"+brep.inputInfo.get(names.get(i)).homologue;
				bw.write(ident+"\t"+tsm.get(names.get(i))+"\n");
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public static String countREPs(File matchFolder,String name,ArrayList<String> names,BlastREPs brep, File out){
		int totalSpacesMissing=0;
		int maxNumber=1;
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out+"/"+name+"_REPcomparison.out"));
			BufferedWriter bw2=new BufferedWriter(new FileWriter(out+"/"+name+"_missings.out"));
			bw2.write("REP-Number\tNumStrainsWhereSpaceMissing\tNumberOfREPsInSpace\n");
			String strainInfoSelf=name+"_"+brep.inputInfo.get(name).homologue;
			int number=0;
			bw.write("REP-Num\tstrainFrom\t#REPs\tstrainTo\t#REPs\tinterS\tinterE\n");
			while(number<maxNumber){
				int missing=0;
				int selfREPs=0;
				boolean intervalSet=false;
				for(int i=0;i<names.size();i++){
					if(names.get(i).equals(name))continue;
					File matchFile=new File(matchFolder+"/"+name+"/"+names.get(i)+".out");
					Matches matches=new Matches(matchFile);
					maxNumber=matches.size()/2;
					Info interval=new Info(0,0,"");
					String strainInfoOther=names.get(i)+"_"+brep.inputInfo.get(names.get(i)).homologue;
					if(!intervalSet){
						selfREPs=brep.countREPs(name, getIntervalSelf(number,matches));
						intervalSet=true;
					}

					if(getIntervalOther(number,matches,interval)){
						bw.write(number+"\t"+strainInfoSelf+"\t"+selfREPs+"\t"+strainInfoOther+"\t"+brep.countREPs(names.get(i), interval)+"\t"+interval.getStart()+"\t"+interval.getEnd()+"\n");
					}else{
						bw.write(number+"\t"+strainInfoSelf+"\tmissing\t"+strainInfoOther+"\tmissing\n");
						missing++;
					}
				}
				if(missing>0){
					totalSpacesMissing++;
				}
				bw2.write(number+"\t"+missing+"\t"+selfREPs+"\n");
				
				number++;
			}
			bw.close();
			bw2.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return totalSpacesMissing+"\t"+maxNumber+"\t"+((1.0*totalSpacesMissing)/maxNumber);
	} 
	
	public static boolean getIntervalOther(int number,Matches mf,Info interval){
		boolean first=true;
		int start=0;
		int end=0;
		for(int i=0;i<mf.size();i++){
			if( i/2==number){
				if(number!=mf.getMatch1(i).number ){
					System.err.println("Something wrong with the matchfile. \nLine number does not correspond to REP number:");
					System.err.println("Line number: "+i+"\t"+mf.getLine(i));
					System.exit(-1);
				}
				if(!mf.getMissing(i)){
					if(first){
						start=mf.getMatch2(i).end;
						first=false;
					}else{
						end=mf.getMatch2(i).start;

					}
				}else{
					return false;
				}
			}
		} 
		interval.info="";
		interval.setStart(start);
		interval.setEnd(end);
		return true;
	}
	public static Info getIntervalSelf(int number,Matches mf){
		boolean first=true;
		int start=0;
		int end=0;
		Info interval;
		for(int i=0;i<mf.size();i++){
			if( i/2==number){
				if(number!=mf.getMatch1(i).number ){
					System.err.println("Something wrong with the matchfile. \nLine number does not correspond to REP number:");
					System.err.println("Line number: "+i+"\t"+mf.getLine(i));
					System.exit(-1);
				}
					if(first){
						start=mf.getMatch1(i).end;
						first=false;
					}else{
						end=mf.getMatch1(i).start;

					}

			}
		} 
		interval=new Info(start,end,"");
		return interval;
	}

	

	
}
