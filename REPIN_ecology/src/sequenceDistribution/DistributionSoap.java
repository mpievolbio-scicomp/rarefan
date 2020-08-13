package sequenceDistribution;

import java.io.File;
import java.util.ArrayList;

import util.Histogram;
import util.PositionAndName;
import util.ReadSoap;

public class DistributionSoap {
	public static void main(String args[]){
		File soap=new File(args[0]);
		File out=new File(args[1]);
		ReadSoap rs=new ReadSoap(soap);

		ArrayList<PositionAndName> pan=toArrayList(rs.getStart(),rs.getQuery());
		
		Integer[] dist=DistanceAndSorting.calcDistance(DistanceAndSorting.sortPAN(pan));
		Histogram<Integer> hist=new Histogram<Integer>(dist);
		hist.write( out,"std");

		
	}
	private static ArrayList<PositionAndName> toArrayList(ArrayList<Integer> start,ArrayList<String> query){
		ArrayList<PositionAndName> pan=new ArrayList<PositionAndName>();
		for(int i=0;i<start.size();i++){
			pan.add(new PositionAndName(start.get(i),query.get(i)));
			
		}
		return pan;
	}

	
}
