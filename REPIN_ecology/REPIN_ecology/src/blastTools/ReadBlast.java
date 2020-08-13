package blastTools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


public class ReadBlast {
	ArrayList<Integer> startS=new ArrayList<Integer>();
	ArrayList<Integer> endS=new ArrayList<Integer>();
	ArrayList<String> query=new ArrayList<String>();
	ArrayList<String> database=new ArrayList<String>();
	ArrayList<Double> evalue=new ArrayList<Double>();
	ArrayList<Integer> length=new ArrayList<Integer>();
	ArrayList<Integer> startQ=new ArrayList<Integer>();
	ArrayList<Integer> endQ=new ArrayList<Integer>();
	ArrayList<Integer> identities=new ArrayList<Integer>();
	ArrayList<Double> score=new ArrayList<Double>();
	ArrayList<Double> percentIdentity=new ArrayList<Double>();
	int pos=0;
	int size=0;
	
	public ReadBlast nextQuery(){
		ReadBlast temp=new ReadBlast();
		String current=query.get(pos);
		while(pos<size&&current.equals(query.get(pos))){
			temp.score.add(this.score.get(pos));
			temp.startS.add(this.startS.get(pos));
			temp.endS.add(this.endS.get(pos));
			temp.query.add(this.query.get(pos));
			temp.database.add(this.database.get(pos));
			temp.evalue.add(this.evalue.get(pos));
			temp.length.add(this.length.get(pos));
			temp.startQ.add(this.startQ.get(pos));
			temp.endQ.add(this.endQ.get(pos));
			temp.identities.add(this.identities.get(pos));
			temp.percentIdentity.add(this.percentIdentity.get(pos));

			pos++;
		}
		return temp;
	}
	public ReadBlast(double score,int startS,int endS,String query,String database,double evalue,int length,int startQ,int endQ,int identities,double percentIdentity){
		pos=0;
		size=1;
		this.score.add(score);
		this.startS.add(startS);
		this.endS.add(endS);
		this.query.add(query);
		this.database.add(database);
		this.evalue.add(evalue);
		this.length.add(length);
		this.startQ.add(startQ);
		this.endQ.add(endQ);
		this.identities.add(identities);
		this.percentIdentity.add(percentIdentity);

	}

	public ReadBlast nextQueryFirstEntry(){
		ReadBlast temp=new ReadBlast();
		String current=query.get(pos);
		temp.score.add(this.score.get(pos));
		temp.startS.add(this.startS.get(pos));
		temp.endS.add(this.endS.get(pos));
		temp.query.add(this.query.get(pos));
		temp.database.add(this.database.get(pos));
		temp.evalue.add(this.evalue.get(pos));
		temp.length.add(this.length.get(pos));
		temp.startQ.add(this.startQ.get(pos));
		temp.endQ.add(this.endQ.get(pos));
		temp.identities.add(this.identities.get(pos));
		temp.percentIdentity.add(this.percentIdentity.get(pos));
		while(pos<size&&current.equals(query.get(pos))){

			pos++;
		}
		return temp;
	}
	
	public boolean hasNext(){
		return pos<size;
	}
	
	public void reset(){
		pos=0;
	}
	
	public HashMap<String,ReadBlast> getQueryHash(){
		HashMap<String,ReadBlast> hm=new HashMap<String, ReadBlast>();
		while(this.hasNext()){
			
			hm.put(query.get(pos), this.nextQuery());
		}
		this.reset();
		return hm;
	}
	public HashMap<String,ReadBlast> getQueryHashFirstEntry(){
		HashMap<String,ReadBlast> hm=new HashMap<String, ReadBlast>();
		while(this.hasNext()){
			
			hm.put(query.get(pos), this.nextQueryFirstEntry());
		}
		this.reset();
		return hm;
	}
	public ReadBlast(){
		pos=0;
		size=0;
	}
	
	private void getLines(File f){
		try{	
			BufferedReader br = new BufferedReader(new FileReader(f));
			String line="";
			while((line=br.readLine())!=null){
				String[] split=line.split("\\s+");
				score.add(Double.parseDouble(split[11]));
				startS.add(Integer.parseInt(split[8]));
				endS.add(Integer.parseInt(split[9]));	
				startQ.add(Integer.parseInt(split[6]));
				endQ.add(Integer.parseInt(split[7]));
				query.add(split[0]);
				evalue.add(Double.parseDouble(split[10]));
				database.add(split[1]);
				length.add(Integer.parseInt(split[3]));
				Double d=(Double.parseDouble(split[2])/100)*length.get(length.size()-1);
				identities.add(d.intValue());
				percentIdentity.add(Double.parseDouble(split[2]));
				size++;
			}
			br.close();
		}catch(IOException e){
			System.err.println(e.toString());
		}
	}
	
	public ReadBlast(File f){
		getLines(f);
	}
	public ArrayList<String> getQuery(){
		return query;
	}
	
	public ArrayList<Double> getEvalue(){
		return evalue;
	}
	
	public ArrayList<String> getDatabase(){
		return database;
	}
	public ArrayList<Integer> getIdentities(){
		return identities;
	}
	public ArrayList<Integer> getLength(){
		return length;
	}
	public ArrayList<Integer> getStartQuery(){
		return startQ;
	}
	
	public ArrayList<Integer> getEndQuery(){
		return endQ;
	}
	public ArrayList<Integer> getStartDB(){
		return startS;
	}
	public ArrayList<Double> getScore(){
		return score;
	}
	public ArrayList<Double> getPercentIdentity(){
		return percentIdentity;
	}
	public ArrayList<Integer> getEndDB(){
		return endS;
	}
	public String get(int i){
		return query.get(i)+"\t"+database.get(i)+"\t"+identities.get(i)+"\t"+startQ.get(i)+"\t"+endQ.get(i)+"\t"+startS.get(i)+"\t"+endS.get(i)+"\t"+evalue.get(i)+"\t"+score.get(i);
	}
}
