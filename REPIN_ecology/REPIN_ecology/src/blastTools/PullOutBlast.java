package blastTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import sequenceDistribution.GetSequence;

public class PullOutBlast {
	public static void main(String args[]){
		File blast=new File(args[0]);
		File genome=new File(args[1]);
		double evalue=Double.parseDouble(args[2]);
		String directory=args[3];
		GetSequence gs=new GetSequence(genome);
		ReadBlast rb=new ReadBlast(blast);
		ArrayList<String> query=rb.getQuery();
		ArrayList<String> database=rb.getDatabase();
		ArrayList<Integer> start=rb.getStartDB();
		ArrayList<Integer> end=rb.getEndDB();
		ArrayList<Double> evalueList=rb.getEvalue();

		String q="";
		if(query.size()>0){
			File f=new File(directory+"/"+query.get(0)+".fasta");
			try{
				BufferedWriter bw=new BufferedWriter(new FileWriter(f));
				for(int i=0;i<query.size();i++){
					if(!q.equals(query.get(i))){
						q=query.get(i);
						f=new File(directory+"/"+q+".fasta");
						bw.close();
						bw=new BufferedWriter(new FileWriter(f));
					}
					if(evalueList.get(i)<=evalue){
						bw.write(">"+start.get(i)+"_"+end.get(i)+"\n");
						
						bw.write(gs.getSequence(">"+database.get(i), start.get(i),end.get(i))+"\n");
					}
				}
				bw.close();
			}catch(IOException e){
				System.err.println(e.toString());
			}
		}
	}
}
