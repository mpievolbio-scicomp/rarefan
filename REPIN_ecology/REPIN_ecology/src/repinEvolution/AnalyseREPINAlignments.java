package repinEvolution;

import java.io.*;
import java.util.*;

import util.*;
import util.phylogenetics.*;

public class AnalyseREPINAlignments {
	HashMap<String,ArrayList<Fasta>> fasMap;
	HashMap<Double,Integer> propNotEqual;
	HashMap<String,ArrayList<Fasta>> alignments;
	public AnalyseREPINAlignments(HashMap<String,ArrayList<Fasta>> fasMap) {
		this.fasMap=fasMap;
		calcAlignments();
		calcPropNotEqual();
	}
	
	public HashMap<Double,Integer> getPropNotEqual(){
		return propNotEqual;
	}
	
	public void writePropNotEqual(File out) {
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			Double[] props=propNotEqual.keySet().toArray(new Double[0]);
			for(int i=0;i<props.length;i++) {
				bw.write(props[i]+"\t"+propNotEqual.get(props[i])+"\n");
			}
			bw.close();
		}catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private void calcPropNotEqual() {
		String[] names=alignments.keySet().toArray(new String[0]);
		propNotEqual=new HashMap<Double, Integer>();
		for(int i=0;i<names.length;i++) {
			ArrayList<Fasta> fas=alignments.get(names[i]);
			if(fas.size()>1) {
				double propDiv=getPropDiv(fas);
				if(!propNotEqual.containsKey(propDiv)) {
					propNotEqual.put(propDiv, 0);
				}
				propNotEqual.put(propDiv, propNotEqual.get(propDiv)+1);
			}
		}
	}
	
	private double getPropDiv(ArrayList<Fasta> fas) {
		Alignment al=new Alignment(fas);
		int differences=0;
		for(int i=0;i<al.getLength();i++) {
			String col=al.getColumn(i);
			char c=col.charAt(0);
			boolean isDiff=false;
			for(int j=0;j<col.length();j++) {
				if(col.charAt(j)!=c) {
					isDiff=true;
					j=col.length();
				}
			}
			if(isDiff) {
				differences++;
			}
		}
		return((1.0*differences)/al.getLength());
	}
	
	private void calcAlignments(){
		alignments=new HashMap<String, ArrayList<Fasta>>();
		String[] names=fasMap.keySet().toArray(new String[0]);
		for(int i=0;i<names.length;i++) {
			ArrayList<Fasta> fas=fasMap.get(names[i]);
			if(sameLength(fas)) {
				alignments.put(names[i],fas);
			}
		}
	}
	
	private boolean sameLength(ArrayList<Fasta> fas) {
		int length=fas.get(0).getSequence().length();
		for(int i=1;i<fas.size();i++) {
			if(fas.get(i).getSequence().length()!=length) {
				return false;
			}
		}
		return true;
	}
	
	
}
