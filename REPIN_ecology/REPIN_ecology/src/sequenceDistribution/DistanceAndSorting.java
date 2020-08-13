package sequenceDistribution;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map.Entry;

import util.PositionAndName;

public class DistanceAndSorting implements Comparator<PositionAndName>{
	public int compare(PositionAndName x,PositionAndName y){
		return x.position==y.position && x.name.equals(y.name)?0:x.position==y.position && !x.name.equals(y.name)?1:x.position>y.position?1:-1;
		
	}
	public static Integer[] sort(ArrayList<Integer> list){
		TreeMap<Integer,Integer> sort=new TreeMap<Integer,Integer>();
		for(int i=0;i<list.size();i++){
			if(!sort.containsKey(list.get(i))){
				sort.put(list.get(i),1);
			}else{
				sort.put(list.get(i),sort.get(list.get(i))+1);
			}
			
		}
		//Integer[] array=sort.keySet().toArray(new Integer[0]);
		
		return toArray(sort);
	}
	
	public static PositionAndName[] sortPAN(ArrayList<PositionAndName> list){
		TreeMap<PositionAndName,Integer> sort=new TreeMap<PositionAndName,Integer>(new DistanceAndSorting());
		for(int i=0;i<list.size();i++){
			if(!sort.containsKey(list.get(i))){
				sort.put(list.get(i),1);
			}else{
				sort.put(list.get(i),sort.get(list.get(i))+1);
			}
			
		}
		//Integer[] array=sort.keySet().toArray(new Integer[0]);
		
		return toArrayPAN(sort);
	}
	
	private static Integer[] toArray(TreeMap<Integer,Integer> tm){
		Iterator<Entry<Integer,Integer>> it=tm.entrySet().iterator();
		ArrayList<Integer> a=new ArrayList<Integer>();
		while(it.hasNext()){
			Entry<Integer,Integer> e=it.next();
			for(int i=0;i<e.getValue();i++){
				a.add(e.getKey());
			}
		}
		return a.toArray(new Integer[0]);
	}
	
	private static PositionAndName[] toArrayPAN(TreeMap<PositionAndName,Integer> tm){
		Iterator<Entry<PositionAndName,Integer>> it=tm.entrySet().iterator();
		ArrayList<PositionAndName> a=new ArrayList<PositionAndName>();
		while(it.hasNext()){
			Entry<PositionAndName,Integer> e=it.next();
			for(int i=0;i<e.getValue();i++){
				a.add(e.getKey());
			}
		}
		return a.toArray(new PositionAndName[0]);
	}
	public static Integer[] calcDistance(Integer[] start){
		
		ArrayList<Integer> dist=new ArrayList<Integer>();
		for(int i=0;i<start.length-1;i++){
			dist.add(start[i+1]-start[i]);
			if(dist.get(dist.size()-1)<=1){
				System.out.println(start[i+1]+" "+start[i]);
			}
		}
		return dist.toArray(new Integer[0]);
	}
	public static Integer[] calcDistance(PositionAndName[] start){
		ArrayList<Integer> dist=new ArrayList<Integer>();
		for(int i=0;i<start.length-1;i++){
			dist.add(start[i+1].position-start[i].position);
			if(dist.get(dist.size()-1)<=200){
			}
		}
		return dist.toArray(new Integer[0]);
	}
}
