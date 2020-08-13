package util;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class ConsenseOutfile {
	
	public static void main(String args[]){
		File inFolder=new File(args[0]);
		File out=new File(inFolder+"/bootstrap.txt");
		File[] list=inFolder.listFiles();
		for(int i=0;i<list.length;i++){
			if(list[i].getName().startsWith("percent")&&list[i].isDirectory()){
				File list2[]=list[i].listFiles();
				int percent=Integer.parseInt(list[i].getName().split("percent")[1]);
				for(int j=0;j<list2.length;j++){
					if(list2[j].getName().equals("outfile")){
						ConsenseOutfile co=new ConsenseOutfile(list2[j]);
						RefragmentAlignment.printBootstrapResults(-1,percent,0,co.getBranchPointsExcluded(),co.getBranchPointsIncluded(),out);

					}
				}
			}
		}

	}
	public ConsenseOutfile(File in){
		readFile(in);
	}
	int branchIn=0;
	int branchOut=0;
	
	public int getBranchPointsIncluded(){
		return branchIn;
	}
	public int getBranchPointsExcluded(){
		return branchOut;
	}
	private void readFile(File in){
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			boolean branchesIncluded=false;
			boolean branchesNotIncluded=false;
			 
			while((line=br.readLine())!=null){
				if(line.matches("^Sets included in the consensus tree.*")){
					branchesIncluded=true;
				}else if(line.matches("^Sets NOT included in consensus tree:.*")){
					branchesNotIncluded=true;
					branchesIncluded=false;
				}else if(branchesIncluded){
					int number=getNumberTrees(line);
					branchIn+=number;
					//System.out.println(line+" "+number);
					//if(number>0)System.out.println(line+" "+number);

					
				}else if(branchesNotIncluded){
					int number=getNumberTrees(line);
					branchOut+=number;
					//if(number>0)System.out.println(line+" "+number);
				}
			}
			
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
	public int getNumberTrees(String line){
		int number=0;
		if(line.matches("^[[\\.|\\*]+\\s]+\\d+.*")){
			Matcher m=Pattern.compile("^[[\\.|\\*]+\\s]+(\\d+).*").matcher(line);
			if(m.find()){
				number=Integer.parseInt(m.group(1));
			}else{
				System.err.println("There seems to be a problem with the following line:");
				System.err.println(line);
			}
		}
		return number;
	}
	
}
