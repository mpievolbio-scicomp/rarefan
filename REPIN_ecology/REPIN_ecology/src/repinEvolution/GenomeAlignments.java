package repinEvolution;

import java.util.*;

public class GenomeAlignments {
	HashMap<String,HashMap<String,AlignTwoGenomes>> alignments=new HashMap<String, HashMap<String,AlignTwoGenomes>>();
	public GenomeAlignments() {
		alignments=new HashMap<String, HashMap<String,AlignTwoGenomes>>();
	}
	
	public void put(String g1,String g2,AlignTwoGenomes atg) {
		if(!alignments.containsKey(g1)) {
			alignments.put(g1, new HashMap<String,AlignTwoGenomes>());
		}
		alignments.get(g1).put(g2, atg);
		if(!alignments.containsKey(g2)) {
			alignments.put(g2, new HashMap<String, AlignTwoGenomes>());
		}
		alignments.get(g2).put(g1, atg);
	}
	
	public AlignTwoGenomes get(String g1,String g2) {
		return alignments.get(g1).get(g2);
	}
	
	public String[] getAllKeys() {
		return alignments.keySet().toArray(new String[0]);
	}
	
}
