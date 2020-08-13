package ecoliREP;

import java.io.*;
import java.util.*;

import util.*;
//Datatype for reading match files
public class Matches {
	ArrayList<Match> matches;
	
	public static void main(String args[]){
		File in=new File(args[0]);
		File out=new File(args[1]);
		Matches m=new Matches(in);
		Info.write(m.findInversions(), out);
	}
	
	public Matches(File in){
		readMatchFile(in);
	}
	
	
	
	public int size(){
		return matches.size();
	}
	
	public String getLine(int index){
		return matches.get(index).toString();
	}
	
	public MatchPart getMatch1(int index){
		return matches.get(index).m1;
	}
	public boolean getMissing(int index){
		return matches.get(index).missing;
	}
	public MatchPart getMatch2(int index){
		return matches.get(index).m2;
	}
	
	public ArrayList<Info> findDuplications(){
		ArrayList<Info> duplicates=new ArrayList<Info>();
		HashMap<String,Boolean> dup=new HashMap<String, Boolean>();
		boolean dupRegion=false;
		int start=-1;
		int end=-1;
		for(int i=0;i<matches.size();i++){
			if(!getMissing(i) &&!dup.containsKey(matches.get(i).m2.toString())){			
				dup.put(matches.get(i).m2.toString(), true);
				if(dupRegion){
					duplicates.add(new Info(start,end,"duplicate"));
				}
				dupRegion=false;
			}else if(!getMissing(i)){
				if(dupRegion)end=matches.get(i).m1.end;
				else {
					start=matches.get(i).m1.start;
					end=matches.get(i).m1.end;
				}
				dupRegion=true;
			}else {
				if(dupRegion)duplicates.add(new Info(start,end,"duplicate"));
				dupRegion=false;
			}
		}
		if(dupRegion){
			duplicates.add(new Info(start,end,"duplicate"));
		}
		return duplicates;
	}
	
	public ArrayList<Info> findInversions(){
		ArrayList<Info> inversions=new ArrayList<Info>();
		int start=-1;
		int end=-1;
		boolean inv=false;
		int signegsum=0;
		int sigpossum=0;
		boolean invshort=false;
		for(int i=1;i<matches.size();i++){
			//at least 3 consecutive inversions to be counted as an inversion
			if(!getMissing(i) && !getMissing(i-1)){
				
				int diff=matches.get(i).m2.number-matches.get(i-1).m2.number;
				double sig=Math.signum(diff);
				if(sig<0 && diff==-1){
					signegsum++;
					if(inv==false){
						
						if(signegsum>=3){
							sigpossum=0;
							inv=true;
						}
						if(!invshort){
							invshort=true;
							start=matches.get(i-1).m1.start;
						}
					}
					end=matches.get(i).m1.end;

				}else if(sig>0){
					invshort=false;
					sigpossum++;
					if(sigpossum>=3 && inv){
						signegsum=0;
						inv=false;
						
						inversions.add(new Info(start,end,"inversion"));
					}
				}
			}
		}

		return inversions;
	}
	
	public ArrayList<Info> findDeletions(){
		ArrayList<Info> deletions=new ArrayList<Info>();
		int start=-1;
		int end=-1;
		boolean missing=false;
		for(int i=0;i<matches.size();i++){
			if(getMissing(i)){
				if(!missing){
					start=matches.get(i).m1.start;
					missing=true;
				}
				end=matches.get(i).m1.end;
			}else{
				if(missing){
					deletions.add(new Info(start,end,"deletion"));
					missing=false;
				}
				
			}
				
		}
		
		return deletions;
	}
	
	private void readMatchFile(File in){
		ArrayList<Match> al=new ArrayList<Match>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+|_");
				if(split.length>=9)al.add(new Match(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),Integer.parseInt(split[3]),split[4],Integer.parseInt(split[5]),Integer.parseInt(split[6]),Integer.parseInt(split[7]),Double.parseDouble(split[8])));
				else al.add(new Match(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),Integer.parseInt(split[3])));
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		matches=al;
	}
	public static class MatchPart{
		String name;
		int start;
		int end;
		int number;
		double eValue;
		public MatchPart(String Name,int Start,int End,int Number){
			
			name=Name;
			start=Start;
			end=End;
			number=Number;

		}
		public String toString(){
				
				return name+"_"+start+"_"+end+"_"+number;
			
		}
		
	}
	public static class Match{
		private boolean missing;
		MatchPart m1=null;
		MatchPart m2=null;

		
		
		public Match(String Name1,int Start1,int End1,int Number1,String Name2,int Start2,int End2,int Number2,double EValue){
			m1=new MatchPart(Name1,Start1,End1,Number1);
			m2=new MatchPart(Name2,Start2,End2,Number2);
			m2.eValue=EValue;
			missing=false;
		}

		public Match(String Name1,int Start1,int End1,int Number1){
			m1=new MatchPart(Name1,Start1,End1,Number1);
			missing=true;
		}
		
		public String toString(){
			if(!missing){
				return m1.name+"_"+m1.start+"_"+m1.end+"_"+m1.number+"\t"+m2.name+"_"+m2.start+"_"+m2.end+"_"+m2.number+"\t"+m2.eValue;
			}else{
				return m1.name+"_"+m1.start+"_"+m1.end+"_"+m1.number+"\tmissing";
			}
		}
		

		
	}
}
