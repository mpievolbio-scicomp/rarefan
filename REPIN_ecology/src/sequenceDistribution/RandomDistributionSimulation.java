package sequenceDistribution;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;
public class RandomDistributionSimulation {
	//TODO percentage in intragenic regions!!!
	BitSet genomeInter;
	BitSet genomeIntra;
	ArrayList<Integer> posList;
	
	public RandomDistributionSimulation(BitSet sequence,int numberOfElements,int length,double percentIntra){
		
		genomeInter=(BitSet)sequence.clone();
		genomeIntra=(BitSet)sequence.clone();
		posList=setElementsRandomly(numberOfElements,length,percentIntra);
	}
	
	public ArrayList<Integer> getPosList(){
		return posList;
	}
	
	private ArrayList<Integer> setElementsRandomly(int nOE,int length,double percentIntra){
		ArrayList<Integer> posList=new ArrayList<Integer>();
		int i=nOE;
		Random r=new Random();
		while(i>0){
			
			int genomePos=r.nextInt(genomeInter.size()-length);
			boolean intra=r.nextDouble()<=percentIntra;
			if(isFree(genomePos,length,intra)){
				posList.add(genomePos);
				setOccupied(genomePos,length,intra);
				i--;
			}
		}
		return posList;
	}
	

	private boolean isFree(int pos,int length,boolean intra){
		for(int i=pos;i<pos+length;i++){
			if(intra?!genomeIntra.get(i):genomeInter.get(i)){
				return false;
			}
		}
		return true;
	}
	
	private void setOccupied(int pos,int length,boolean intra){
		for(int i=pos;i<pos+length;i++){

			if(intra)genomeIntra.set(i,false);
			else genomeInter.set(i,true);

		}

	}
	
	public static ArrayList<Integer> createFreeSpaces(boolean genome[]){
		ArrayList<Integer> free=new ArrayList<Integer>();
		for(int i=0;i<genome.length;i++){
			if(!genome[i]){
				free.add(i);
			}
		}
		return free;
	}
	
	
}
