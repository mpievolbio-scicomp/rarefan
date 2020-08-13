package frequencies.util;

import java.io.File;
import java.util.*;

import statistics.Stats;

public class CalculateAverageNetworkProperties {
	NetworkProperties avg;
	NetworkProperties stddev;
	public CalculateAverageNetworkProperties(ArrayList<NetworkProperties> list){
		
		ArrayList<Stats> propStats=getPropertyStats(list);
		String regime=list.get(0).getRegime();
		String genomeID=list.get(0).genomeID;
		File outFolder=list.get(0).outFolder.getParentFile();
		avg=new NetworkProperties(regime,genomeID,outFolder);
		stddev=new NetworkProperties(regime,genomeID,outFolder);
		avg.setAll(getAvg(propStats));
		stddev.setAll(getStdDev(propStats));
		avg.regime=avg.regime+"_avg";
		stddev.regime=stddev.regime+"_stddev";
	}
	
	public NetworkProperties getAvg(){
		return avg;
	}
	
	public NetworkProperties getStdDev(){
		return stddev;
	}
	
	private ArrayList<Double> getAvg(ArrayList<Stats> propStats){
		ArrayList<Double> avg=new ArrayList<Double>();
		for(int i=0;i<propStats.size();i++){
			avg.add(propStats.get(i).getAverage());
		}
		return avg;
	}
	private ArrayList<Double> getStdDev(ArrayList<Stats> propStats){
		ArrayList<Double> avg=new ArrayList<Double>();
		for(int i=0;i<propStats.size();i++){
			avg.add(propStats.get(i).getStandardDeviation());
		}
		return avg;
	}

	private ArrayList<Stats> getPropertyStats(ArrayList<NetworkProperties> list){
		ArrayList<ArrayList<Double>> listOfProperties=convertPropertyList(list);
		ArrayList<Stats> propStats=new ArrayList<Stats>();
		for(int i=0;i<listOfProperties.size();i++){
			propStats.add(new Stats(listOfProperties.get(i)));
		}
		return propStats;
	}
	
	private ArrayList<ArrayList<Double>> convertPropertyList(ArrayList<NetworkProperties> list){
		ArrayList<ArrayList<Double>> converted=new ArrayList<ArrayList<Double>>();
		for(int i=0;i<list.size();i++){
			ArrayList<Double> properties=list.get(i).getAll();
			for(int j=0;j<properties.size();j++){
				if(converted.size()<=j){
					converted.add(new ArrayList<Double>());
				}
				converted.get(j).add(properties.get(j));
			}
		}
		
		return converted;
	}
	
	
}
