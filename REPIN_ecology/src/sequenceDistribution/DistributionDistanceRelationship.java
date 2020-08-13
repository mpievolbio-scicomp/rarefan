package sequenceDistribution;
//obsolete!!!!
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map.Entry;

import util.Histogram;
import util.PositionAndName;
import util.ReadSoap;

public class DistributionDistanceRelationship {
	public static void main(String args[]){
		File soap=new File(args[0]);
		File out=new File(args[1]);
		int upperBound=Integer.parseInt(args[2]);
		boolean distance = Boolean.parseBoolean(args[3]);
		File cluster=new File(args[4]);
		ReadSoap rs=new ReadSoap(soap);
		Histogram<String> h=new Histogram<String>(rs.getQuery().toArray(new String[0]));
		HashMap<String,String> trans=groupTranslation(h.sortByValues());
		ArrayList<PositionAndName> pan=toArrayList(rs.getStart(),rs.getQuery(),trans);
		Histogram<String> h2=new Histogram<String>(addEntries(DistanceAndSorting.sortPAN(pan),upperBound,distance,cluster));
		h2.writeSortedByValues(out);

		
	}
	
	private static HashMap<String,String> groupTranslation(TreeMap<Double,ArrayList<String>> highscore){
		HashMap<String,String> trans=new HashMap<String, String>();
		Iterator<Entry<Double,ArrayList<String>>> it=highscore.descendingMap().entrySet().iterator();
		int i=0;
		while(it.hasNext()){
			Entry<Double,ArrayList<String>> e=it.next();
			i++;
			for(int j=0;j<e.getValue().size();j++){
				trans.put(e.getValue().get(j), "Group "+i);
			}
		}
		
		return trans;
	}
	
	private static TreeMap<String,Double> addEntries(PositionAndName[] sortedPan,int upperBound,boolean distance,File cluster){
		TreeMap<String, Double> h=new TreeMap<String, Double>();
		int j=0;
		ArrayList<Integer> elements=new ArrayList<Integer>();
		for(int i=0;i<sortedPan.length-1;i++){
			int dist=sortedPan[i+1].position-sortedPan[i].position;
			String name1=distance?sortedPan[i+1].name+"\t"+sortedPan[i].name+"\t"+dist:sortedPan[i+1].name+"\t"+sortedPan[i].name;
			String name2=distance?sortedPan[i].name+"\t"+sortedPan[i+1].name+"\t"+dist:sortedPan[i].name+"\t"+sortedPan[i+1].name;
			j++;
			if(dist>200){
				elements.add(j);
				if(j>1){
					int prev=0;
					for(int k=-j+1;k<=0;k++){
						System.out.println(sortedPan[i+k].name+" "+sortedPan[i+k].position+" "+(sortedPan[i+k].position-prev));
						prev=sortedPan[i+k].position;
					}
					System.out.println("________________________");
				}
				
				j=0;
				
				
			}
			
			
			if(dist<=upperBound){
				if(h.containsKey(name1)){
					h.put(name1,h.get(name1)+1);
				}else if(h.containsKey(name2)){
					h.put(name2,h.get(name2)+1);

				}else{
					h.put(name2,1.0);
				}
			}
		}
		Histogram<Integer> packet=new Histogram<Integer>(elements.toArray( new Integer[0]));
		packet.write(cluster,"std");
		return h;
	}
	private static ArrayList<PositionAndName> toArrayList(ArrayList<Integer> start,ArrayList<String> query,HashMap<String,String> trans){
		ArrayList<PositionAndName> pan=new ArrayList<PositionAndName>();
		for(int i=0;i<start.size();i++){
			pan.add(new PositionAndName(start.get(i),trans.get(query.get(i))));
			
		}
		return pan;
	}
}
