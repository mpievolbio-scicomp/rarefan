package util;
import java.io.*;
import java.util.*;
public class FindConnectedComponent<T> {
	
	public static void main(String args[]){
		File in=new File(args[0]);
		String node=args[1];
		HashMap<String,ArrayList<String>> graph=readGraph(in);
		FindConnectedComponent<String> fcc=new FindConnectedComponent<String>(graph);
		fcc.printComponents(fcc.getComponent(node));
		
	}
	
	
	HashMap<T,ArrayList<T>> graph=new HashMap<T, ArrayList<T>>();
	/**removes a random component from the graph
	 * 
	 * @return
	 * returns all elements that are part of the removed connected component
	 */
	
	public ArrayList<T> removeConnectedComponent(){
		ArrayList<T> list=new ArrayList<T>();
		if(graph.size()>0){
			Iterator<T> it=graph.keySet().iterator();
			T node=it.next();
			list=getComponent(node);
			removeComponents(list);
		}
		return list;
	}
	
	public static HashMap<String,ArrayList<String>> readGraph(File in){
		HashMap<String,ArrayList<String>> hm=new HashMap<String, ArrayList<String>>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line;
			while((line=br.readLine())!=null){
				String split[]=line.split("\\s+");
				if(!hm.containsKey(split[0])){
					hm.put(split[0],new ArrayList<String>());
				}
				hm.get(split[0]).add(split[1]);
			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return hm;
	}
	
	public FindConnectedComponent(HashMap<T,ArrayList<T>> graph){
		this.graph=graph;
	}
	
	public int size(){
		return graph.size();
	}
	
	public ArrayList<T> getComponent(T node){
		HashMap<T,Boolean> visited=new HashMap<T, Boolean>();
		ArrayList<T> components=new ArrayList<T>();
		recursive(node,visited,components);
		return components;
	}
	
	public void removeComponents(ArrayList<T> nodes){
		for(int i=0;i<nodes.size();i++){
			graph.remove(nodes.get(i));
		}
	}
	
	private void printComponents(ArrayList<T> comp){
		for(int i=0;i<comp.size();i++){
			System.out.println(comp.get(i));
		}
	}
	
	private void recursive(T node,HashMap<T,Boolean> visited,ArrayList<T> components){
		if(graph.containsKey(node)){
			visited.put(node,true);
			components.add(node);
			ArrayList<T> children=graph.get(node);
			for(int i=0;i<children.size();i++){
				T child=children.get(i);
				if(!visited.containsKey(child)){
					recursive(child,visited,components);
				}
			}
		}
	}
	
	
}
