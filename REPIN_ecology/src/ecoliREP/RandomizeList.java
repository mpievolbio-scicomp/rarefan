package ecoliREP;

import java.util.*;
import java.util.Map.Entry;

public class RandomizeList {
	public static HashMap<String,String> randomise(HashMap<String,String> list){
		HashMap<String,String> newList=new HashMap<String, String>();
		Iterator<Entry<String,String>> it=list.entrySet().iterator();
		ArrayList<String> names=new ArrayList<String>();
		ArrayList<String> values=new ArrayList<String>();
		while(it.hasNext()){
			Entry<String,String> e=it.next();
			names.add(e.getKey());
			values.add(e.getValue());
		}
		
		while(names.size()>0){
			double rand=Math.random();
			int pos=(int)(names.size()*rand);
			if(pos==names.size())pos=pos-1;
			String name=names.get(pos);
			String val=values.get(0);
			names.remove(pos);
			values.remove(0);
			newList.put(name,val);
		}
		
		return newList;
	}
}
