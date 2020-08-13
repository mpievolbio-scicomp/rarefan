package util;

import java.util.ArrayList;

public class List_Array {
	public static double[] toDouble(ArrayList<Double> list){
		int size=list.size();
		double[] d=new double[size];
		for(int i=0;i<size;i++){
			d[i]=list.get(i);
		}
		return d;
	}
	public static double[] toDouble(Double[] list){
		int size=list.length;
		double[] d=new double[size];
		for(int i=0;i<size;i++){
			d[i]=list[i];
		}
		return d;
	}
	public static int[] BooleanToInteger(ArrayList<Boolean> list){
		int size=list.size();
		int[] d=new int[size];
		for(int i=0;i<size;i++){
			d[i]=list.get(i)?1:0;
		}
		return d;
	}
	
	public static int[] toInt(Integer[] list){
		int size=list.length;
		int[] d=new int[size];
		for(int i=0;i<size;i++){
			d[i]=list[i];
		}
		return d;
	}
	
	public static int[] toInt(ArrayList<Integer> list){
		int size=list.size();
		int[] d=new int[size];
		for(int i=0;i<size;i++){
			d[i]=list.get(i);
		}
		return d;
	}
	public static ArrayList<Double> toDouble(double[] list){
		int size=list.length;
		ArrayList<Double> d=new ArrayList<Double>();
		for(int i=0;i<size;i++){
			d.add(list[i]);
		}
		return d;
	}
	public static ArrayList<String> toString(String[] list){
		int size=list.length;
		ArrayList<String> d=new ArrayList<String>();
		for(int i=0;i<size;i++){
			d.add(list[i]);
		}
		return d;
	}
	
	public static String[] toString(ArrayList<String> list){
		int size=list.size();
		String[] d=new String[size];
		for(int i=0;i<size;i++){
			d[i]=list.get(i);
		}
		return d;
	}
}
