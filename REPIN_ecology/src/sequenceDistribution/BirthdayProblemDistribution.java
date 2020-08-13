package sequenceDistribution;

import java.util.Random;

public class BirthdayProblemDistribution {
	int[] hm;
	public  BirthdayProblemDistribution(int numberOfElements,int sites){
		setElementsRandomly(numberOfElements,sites);
	}
	private void setElementsRandomly(int nOE,int sites){
		hm=new int[sites];
		int i=nOE;
		Random r=new Random();
		while(i>0){
			int genomePos=r.nextInt(sites);
			hm[genomePos]++;
			i--;
		}
	}
	public int[] getDistribution(){
		return hm;
	}
}
