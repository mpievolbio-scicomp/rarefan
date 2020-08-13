package sequenceDistribution;

import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;

import util.Histogram;
import util.ReadArtemis;

public class LengthDistributionOfExtragenicRegions {
	public static void main(String args[]){
		File artemis=new File(args[0]);
		File out = new File(args[1]);
		File proportionout = new File(args[2]);
		ReadArtemis rA=new ReadArtemis(artemis);
		Histogram<Integer> h=new Histogram<Integer>(convertBitArrayToLengthArray(rA.getBoolArray("CDS")));
		h.write(out,"std");
		h.set(Histogram.calculateProportionIntegral(h.getHistogram()));
		h.write(proportionout,"std");
	}
	

	private static Integer[] convertBitArrayToLengthArray(BitSet pos){
		ArrayList<Integer> elements=new ArrayList<Integer>();
		
		for(int i=0;i<pos.size();i++){
			if(pos.get(i)==false){
				int j=0;
				while(pos.get(i)==false && i<pos.size()){
					j++;
					i++;
				}
				elements.add(j);
			}
		}
		
		return  elements.toArray(new Integer[0]);
	}
	
}
