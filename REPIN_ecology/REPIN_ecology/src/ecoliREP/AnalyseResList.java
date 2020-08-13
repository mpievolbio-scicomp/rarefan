package ecoliREP;

import java.io.*;
import java.util.*;

import ecoliREP.REPdifferenceHistogram.*;

public class AnalyseResList {
	public static void  main(String args[]){
		File resListREPsFile =new File(args[0]);
		File repNumber=new File(args[1]);
		File ExSpaceNumber=new File(args[2]);
		File resListExSpaceFile=new File(args[3]);
		File out=new File(args[4]);
		
		ArrayList<String> list=CompareMatrices.readList(repNumber);
		ArrayList<ResList> resListREPs=REPdifferenceHistogram.readResList(resListREPsFile);
		HashMap<String,HashMap<String, Double>> matrix1ISREPs=makeMatrix(resListREPs,list);
		HashMap<String,HashMap<String, Double>> matrixREPNum=CompareMatrices.readMatrix(repNumber, list);
		HashMap<String,HashMap<String, Double>> repRatio=CompareMatrices.getRatio(matrix1ISREPs, matrixREPNum, list);
		
		ArrayList<ResList> resListExSpace=REPdifferenceHistogram.readResList(resListExSpaceFile);
		HashMap<String,HashMap<String, Double>> matrix1ISExSpaces=makeMatrix(resListExSpace,list);
		HashMap<String,HashMap<String, Double>> matrixExSpaceNum=CompareMatrices.readMatrix(ExSpaceNumber, list);
		HashMap<String,HashMap<String, Double>> ExSpaceRatio=CompareMatrices.getRatio(matrix1ISExSpaces, matrixExSpaceNum, list);
		
		HashMap<String,HashMap<String, Double>> RepExSpaceRatio=CompareMatrices.getRatio(repRatio, ExSpaceRatio, list);
		CompareMatrices.writeMatrix(out, RepExSpaceRatio, list);
		
	}
	
	public static HashMap<String,HashMap<String,Double>> makeMatrix(ArrayList<ResList> resList,ArrayList<String> list){
		HashMap<String,HashMap<String,Double>> matrix=new HashMap<String, HashMap<String,Double>>();
		for(int i=0;i<resList.size();i++){
			String s1=resList.get(i).s1;
			String s2=resList.get(i).s2;
			if(resList.get(i).annotation.contains("IS629")){
				if(matrix.containsKey(s1)){
					if(matrix.get(s1).containsKey(s2)){
						matrix.get(s1).put(s2, matrix.get(s1).get(s2)+1);
					}else{
						matrix.get(s1).put(s2, 1.0);
					}
				}else{
					HashMap<String,Double> temp=new HashMap<String, Double>();
					temp.put(s2,1.0);
					matrix.put(s1,temp);
				}
			}
		}
		for(int i=0;i<list.size();i++){
			for(int j=0;j<list.size();j++){
				String s1=list.get(i);
				String s2=list.get(j);
				if(matrix.containsKey(s1)){
					if(!matrix.get(s1).containsKey(s2)){
						matrix.get(s1).put(s2, 0.0);
					}
				}else{
					HashMap<String,Double> temp=new HashMap<String, Double>();
					temp.put(s2,0.0);
					matrix.put(s1,temp);
				}
			}
		}
		return matrix;
	}
	

	

	
}
