package blastTools;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
public class ReadHmmer {

		ArrayList<Integer> startDB=new ArrayList<Integer>();
		ArrayList<Integer> endDB=new ArrayList<Integer>();
		ArrayList<String> query=new ArrayList<String>();
		ArrayList<String> queryAccession=new ArrayList<String>();

		ArrayList<String> target=new ArrayList<String>();
		ArrayList<String> targetAccession=new ArrayList<String>();
		ArrayList<String> targetDescription=new ArrayList<String>();

		ArrayList<Double> evalue=new ArrayList<Double>();
		ArrayList<Integer> targetLength=new ArrayList<Integer>();
		ArrayList<Integer> queryLength=new ArrayList<Integer>();

		ArrayList<Integer> startQ=new ArrayList<Integer>();
		ArrayList<Integer> endQ=new ArrayList<Integer>();
		ArrayList<Double> score=new ArrayList<Double>();
		ArrayList<Double> bias=new ArrayList<Double>();

		private void getLinesDomain(File f){
			try{	
				BufferedReader br = new BufferedReader(new FileReader(f));
				String line="";
				while((line=br.readLine())!=null){
					if(!line.startsWith("#")){
						String[] split=line.split("\\s+");
						score.add(Double.parseDouble(split[7]));
						bias.add(Double.parseDouble(split[8]));

						startDB.add(Integer.parseInt(split[15]));
						endDB.add(Integer.parseInt(split[16]));	
						startQ.add(Integer.parseInt(split[17]));
						endQ.add(Integer.parseInt(split[18]));
						query.add(split[3]);
						evalue.add(Double.parseDouble(split[6]));
						target.add(split[0]);
						targetAccession.add(split[1]);
						targetDescription.add(split[22]);

						queryAccession.add(split[4]);

						targetLength.add(Integer.parseInt(split[2]));
						queryLength.add(Integer.parseInt(split[5]));

					}
				}
				br.close();
			}catch(IOException e){
				System.err.println(e.toString());
			}
		}
		private void getLinesTable(File f){
			try{	
				BufferedReader br = new BufferedReader(new FileReader(f));
				String line="";
				while((line=br.readLine())!=null){
					if(!line.startsWith("#")){
						String[] split=line.split("\\s+");
						score.add(Double.parseDouble(split[5]));
						bias.add(Double.parseDouble(split[6]));


						query.add(split[2]);
						evalue.add(Double.parseDouble(split[4]));
						target.add(split[0]);
						targetAccession.add(split[1]);

						queryAccession.add(split[3]);


					}
				}
				br.close();
			}catch(IOException e){
				System.err.println(e.toString());
			}
		}
		public ReadHmmer(File f,String type){
			if(type.equals("tblout")){
				getLinesTable(f);
			}else if(type.equals("domtblout")){
				getLinesDomain(f);
			}
		}
		public ArrayList<String> getQuery(){
			return query;
		}
		public ArrayList<String> getQueryAccession(){
			return queryAccession;
		}
		public ArrayList<Double> getEvalue(){
			return evalue;
		}
		
		public ArrayList<String> getTarget(){
			return target;
		}
		public ArrayList<String> getTargetAccession(){
			return targetAccession;
		}
		public ArrayList<String> getTargetDescription(){
			return targetDescription;
		}
		public ArrayList<Integer> getTargetLength(){
			return targetLength;
		}
		
		public ArrayList<Integer> getQueryLength(){
			return queryLength;
		}
		public ArrayList<Integer> getStartQuery(){
			return startQ;
		}
		
		public ArrayList<Integer> getEndQuery(){
			return endQ;
		}
		public ArrayList<Integer> getStartDB(){
			return startDB;
		}
		public ArrayList<Double> getScore(){
			return score;
		}
		public ArrayList<Integer> getEndDB(){
			return endDB;
		}
		public String get(int i){
			return query.get(i)+"\t"+target.get(i)+"\t"+startQ.get(i)+"\t"+endQ.get(i)+"\t"+startDB.get(i)+"\t"+endDB.get(i)+"\t"+evalue.get(i)+"\t"+score.get(i);
		}
}


