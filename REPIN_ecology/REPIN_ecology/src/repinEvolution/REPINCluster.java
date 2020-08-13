package repinEvolution;

import java.util.*;

import util.*;

public class REPINCluster {
	TreeMap<Integer,ArrayList<REPIN>> cluster;
	String refID;
	int proximity;
	
	public static void main(String[] args) {
		REPINCluster rc=new REPINCluster("test", 100);
		rc.addREPIN("test",1,new REPIN(1, 2, 0, new Fasta("test","ATG"), "test",0,0));
	}
	
	public int size() {
		return cluster.size();
	}
	
	public REPINCluster(String referenceID,int proximity) {
		this.proximity=proximity;
		cluster=new TreeMap<Integer, ArrayList<REPIN>>();
		refID=referenceID;
	}
	
	public void addREPIN(String referenceID,int refPos,REPIN repin) {
		if(referenceID.equals(refID)) {
			addREPIN(repin,refPos);
		}else {
			System.err.println("Did not add REPINs, reference IDs did not match.");
		}
	}
	
	private void addREPIN(REPIN repin,int pos) {
		if(cluster.containsKey(pos)) {
			cluster.get(pos).add(repin);
		}else {
			int posHigher=cluster.higherKey(pos)!=null?cluster.higherKey(pos):-1;
			int posLower=cluster.lowerKey(pos)!=null?cluster.lowerKey(pos):-1;
			if(posHigher>-1 && posHigher-pos<proximity) {
				cluster.get(posHigher).add(repin);
			}else if(posLower>-1 && pos-posLower<proximity) {
				cluster.get(posLower).add(repin);
			}else {
				cluster.put(pos, new ArrayList<REPIN>());
				cluster.get(pos).add(repin);
			}
				
		}
	}
	
	public ArrayList<REPIN> getREPINCluster(int pos){
		return cluster.get(pos);
	}
	
	public Integer[] getStart() {
		return cluster.keySet().toArray(new Integer[0]);
	}
	
	
}
