package simulation;

import java.io.*;
import java.util.*;

import util.*;


public class SimulateSequenceEvolution {
	String sequence="CAATGTGGGAGGGGGCTTGCCCCCGATGGCGGTGGTTCAGTCAACGATATTCCAACTGACCCACCGCCATCGGGAGCAAGCCCCCTCCCACCGTTT";
	double mutationRate=5e-11;
	double movingRate=0.0000390625;
    double movingMR=0.01;
	double replicationRate;
	int max=120;
	int mutAdj=1000;
	int movingAdj=1;
	double deathT=0.95;
	int fitlength=24;
	ArrayList<StringBuffer> population=new ArrayList<StringBuffer>();
	ArrayList<Long> generation=new ArrayList<Long>();
	ArrayList<Integer> parents=new ArrayList<Integer>();
	ArrayList<Fasta> fas=new ArrayList<Fasta>();
	ArrayList<Integer> moves=new ArrayList<Integer>();
	boolean directedMutation=false;
	boolean allReplicate=true;
	double recombinationRate=0.001;
	//int recombStart=38;
	//int recombLength=15;
	int recombStart=29;
	int recombLength=15+18;

	boolean recombine=true;
	double GC;
	long currgeneration;
	boolean fitness=false;
	
	
	
	public static void main(String args[]){
		File REPINalignment=new File(args[0]);
		File outFolder=new File(args[1]);
	
		File out=new File(outFolder+"/seqPopulationMove.fas");
		SimulateSequenceEvolution sse=new SimulateSequenceEvolution(REPINalignment,true);
		sse.writePopulation(out);
	}
	
	public void setRecombine(boolean recomb){
		recombine=recomb;
		recombinationRate=replicationRate;

	}
	
	public void setReplicationRate(double rate){
		replicationRate=rate;
		if(recombine)recombinationRate=replicationRate;
	}
	
	public void setAllReplicate(boolean allReplicate){
		this.allReplicate=allReplicate;
	}
	
	public void setMovingRate(double mr){
		this.movingRate=mr;
	}
	
	public void setMovingMutationRate(double mmr){
		this.movingMR=mmr;
	}
	
	public void setDirectedMutation(boolean dm){
		directedMutation=dm;
	}
	
	public SimulateSequenceEvolution(File alg,boolean fitness){
		this.fitness=fitness;
		ArrayList<Fasta> fas=Fasta.readFasta(alg);
		max=fas.size();
		sequence=fas.get(0).getSequence();
		init(1);
		
	}
	
	
	
	private void init(int num){
		for(int i=0;i<num;i++){
			population.add(new StringBuffer(sequence));
			generation.add(0L);
			parents.add(-1);
			moves.add(0);
		}
		GC=getGC();
	}
	
	private double getGC(){
		int count=0;
		for(int i=0;i<sequence.length();i++){
			if(sequence.charAt(i)=='G'||sequence.charAt(i)=='C'){
				count++;
			}
		}
		return count*(1.0)/sequence.length();
	}
	
	public void simulateEvolution(){
		for(long i=0;population.size()<max&&i<25000000000L&&population.size()>0;i++){
			int lengthold=population.size();
			if(i%mutAdj==0){
				mutate();
				replicate();
				move();
				recombineIntern();
			}
			currgeneration++;
			if(lengthold<population.size()){
				System.out.println(population.size()+" out of "+max+" replications. Generation: "+i);
			}
			
		}
	}
	
	
	//fitness decreases linearly for non-pairing basepairs for REPIN
	private double getFitness(String seq){
		if(!fitness){
			return 1;
		}else{
			double maxFit=(seq.length()/2)-2;
			double sum=0;
			for(int i=0;i<seq.length()/2-2;i++){
				char s=seq.charAt(i);
				char e=seq.charAt(seq.length()-1-i);
				sum+=DNAmanipulations.getComplement(s) ==e?1:0;
			}
			return sum/maxFit;
		}
	}
	
	//fitness decreases linearly for each mutation that occurs in the last 24 and first 24 base pairs
	private double getFitness2(String seq){
		if(!fitness){
			return 1;
		}else{
			double maxFit=fitlength*2;
			int seqlength=seq.length();
			int first=getIdentity(seq.substring(0,fitlength),sequence.substring(0,fitlength));
			int last=getIdentity(seq.substring(seqlength-fitlength,seqlength),sequence.substring(seqlength-fitlength,seqlength));
			
			return (first+last)/maxFit;
		}
	}
	
	private int getIdentity(String seq1,String seq2){
		int sum=0;
		for(int i=0;i<seq1.length();i++){
			char s=seq1.charAt(i);
			char e=seq2.charAt(i);
			sum+=s==e?1:0;
		}
		return sum;
	}
	
	private void recombineIntern(){
		for(int i=0;i<population.size();i++){
			double p=Math.random();
			if(p<=recombinationRate*mutAdj){
				recombineIntern(population.get(i));
			}
		}
	}
	
	private void recombineExtern(){
		for(int i=0;i<population.size();i++){
			double p=Math.random();
			if(p<=recombinationRate*mutAdj){
				recombineRandom(population.get(i));
			}
		}
	}
	
	private void recombineIntern(StringBuffer sb){
		int recombStart=(int)(Math.random()*(sb.length()/2));
		int recombEnd=sb.length()-recombStart;
		for(int i=recombStart;i<recombEnd;i++){
			char c=sb.charAt(i);
			int index=sb.length()-i-1;
			sb.setCharAt(index,DNAmanipulations.getComplement(c));
		}
	}
	
	private void recombineRandom(StringBuffer sb){
		for(int i=recombStart;i<recombStart+recombLength;i++){
			sb.setCharAt(i,DNAmanipulations.getRandomBase(0.4));
		}
	}
	
	private void replicate(){
		ArrayList<String> add=new ArrayList<String>();
		ArrayList<Integer> parent=new ArrayList<Integer>();
		for(int i=0;i<population.size();i++){
			double p=Math.random();
			double fit=getFitness2(population.get(i).toString());
			if(fit<deathT){
				remove(i);
			}else{
				
				if(p<=replicationRate*mutAdj*((fit-deathT)/(1-deathT))){
					add.add(population.get(i).toString());
					parent.add(i);
					moves.add(0);
				}
			}
			if(!allReplicate){
				break;
			}
		}
		addPopulation(add,parent);
	}
	
	private void move(){
		for(int i=0;i<population.size();i++){
			double p=Math.random();
			if(p<=movingRate*mutAdj){
				if(!directedMutation)mutate(population.get(i),mutAdj,movingMR,0);
				else mutate(population.get(i),mutAdj,movingMR,fitlength);
				moves.set(i,moves.get(i)+1);
			}
		}
	}
	
	private void remove(int index ){
		population.remove(index);
		moves.remove(index);
		generation.remove(index);
		parents.remove(index);
	}
	
	

	
	private void addPopulation(ArrayList<String> add,ArrayList<Integer> parent){
		for(int i=0;i<add.size();i++){
			population.add(new StringBuffer(add.get(i)));
			generation.add(currgeneration);
			parents.add(parent.get(i));
		}
	}
	
	private void mutate(){
		for(int i=0;i<population.size();i++){
			StringBuffer seq=population.get(i);
			mutate(seq,mutAdj,mutationRate,0);
		}
	}
	
	private void mutate(StringBuffer seq,int mutAdj,double mutationRate,int fitlength){
		for(int i=fitlength;i<seq.length()-fitlength;i++){
			double p=Math.random();
			if(p<mutationRate*mutAdj){
				char c=seq.charAt(i);
				seq.setCharAt(i, DNAmanipulations.getMutation(c, GC));
			}
		}
	}
	

	
	private ArrayList<Fasta> makeFasta(){
		String moveParam="movingRate:"+movingRate+";movingMutationRate:"+movingMR;
		String fitnessParam=fitness?"deathThreshold:"+deathT+";fitnessLength:"+fitlength:"";
		String parameters="RepRate:"+replicationRate+";"+moveParam+";"+fitnessParam;
		ArrayList<Fasta> fas=new ArrayList<Fasta>();
		for(int i=0;i<population.size();i++){
			fas.add(new Fasta(i+"_gen:"+generation.get(i)+"_parent:"+parents.get(i)+" "+parameters,population.get(i).toString()));
		}
		return fas;
	}
	
	public ArrayList<Fasta> getPopulation(){
		
		return makeFasta();
	}
	public void writePopulation(File out){
		Fasta.write(makeFasta(), out);
	}
}
