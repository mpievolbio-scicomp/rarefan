package util;

import java.io.*;
import java.util.*;

public class MapOperations {
	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue( Map<K, V> map ){
		List<Map.Entry<K, V>> list =new LinkedList<Map.Entry<K, V>>( map.entrySet() );
		Collections.sort( list, new Comparator<Map.Entry<K, V>>(){
			public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 ){
				return -(o1.getValue()).compareTo( o2.getValue() );
			}
		} );

		Map<K, V> result = new LinkedHashMap<K, V>();
		for (Map.Entry<K, V> entry : list)
		{
			result.put( entry.getKey(), entry.getValue() );
		}
		return result;
	}
	
	public static <K> int getSumEntries(Map<K,Integer> map){
		Iterator<Integer> it=map.values().iterator();
		int sum=0;
		while(it.hasNext()){
			sum+=it.next();
		}
		return sum;
	}
	
	public static <K, V extends Comparable<? super V>> void write(Map<K,V> map,File out){
		try{
			Iterator<K> it=map.keySet().iterator();
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			while(it.hasNext()){
				K k=it.next();
				bw.write(k+"\t"+map.get(k)+"\n");
			}
			bw.close();
			
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static <K> void add(Map<K,Integer> map,K key,Integer number){
		if(!map.containsKey(key)){
			map.put(key, number);
		}else{
			map.put(key,map.get(key)+number);
		}
	}

}
