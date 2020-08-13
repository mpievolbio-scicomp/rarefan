package simulation;

import java.io.*;
import java.util.*;

import util.*;
/*TODO!
 * positions of mutations between all, between next neighbors, and between next neighbors in clusters and between clusters
 * 
 */
public class ReplicationRateInference {
	public static void main(String args[]){
		File REPINalg=new File(args[0]);
		//int numTrials=Integer.parseInt(args[1]);
		File outFolder=new File(args[1]);
		int id=Integer.parseInt(args[2]);
		double repAdj=Double.parseDouble(args[3]);
		//for(int id=1;id<100;id++){
		ReplicationRateInference rri=new ReplicationRateInference(REPINalg,outFolder,id);
		//rri.runWithFitness();
		rri.setAdjRepR(repAdj);
		//rri.setAllReplicate(false);
		rri.setRecombine(true);
		rri.runWithoutFitness();
		//rri.setMove();
		//rri.setAdjMoveMR(0.01);
		//rri.runWithoutFitness();
		//rri.setAdjMoveMR(1);
		//rri.setAdjRepR(10);
		//rri.runWithFitness();
		//rri.setAdjMoveMR(1000);
		//rri.setAdjRepR(1000);
		//rri.runWithFitnessDirected();
		//}
	}
	
	//int trials;
	File outFolder;
	File alignment;
	boolean overwrite=false;
	boolean move=false;
	int id;
	double adjMoveMR=1;
	double adjRepR=1;
	boolean allReplicate=true;
	boolean recombine=false;
	
	public void setRecombine(boolean recomb){
		recombine=recomb;
	}
	
	public void setAllReplicate(boolean set){
		allReplicate=set;
	}
	
	public void setAdjMoveMR(double adj){
		adjMoveMR=adj;
	}
	public void setAdjRepR(double adj){
		adjRepR=adj;
	}
	public void setMove(){
		move=true;
	}
	
	public void setOverwrite(){
		overwrite=true;
	}

	public void runWithoutFitness(){
		runSimulation(false,"NoFitness",false);
	}

	public void runWithFitness(){
		runSimulation(true,"Fitness",false);
	}

	public void runWithFitnessDirected(){
		runSimulation(true,"FitnessDirected",true);
	}
	
	private void runSimulation(boolean fitness,String suffix,boolean directed){
		String moveStr=move?"Move":"";


		int subID=id%100;
		//repRate needs to run from 5e-10 to 5e-9; 10 repeats
		int repID=subID%10;
		double repRate=(repID)*5*Math.pow(10, -10)*adjRepR;
		//for each of those the moveMR needs to run from 1e-2 to 1e-9; 10 repeats
		int moveID=subID/10+2;
		double moveMR=Math.pow(10, -moveID)*adjMoveMR;
		File tempFolder=new File(outFolder+"/"+suffix+moveStr+"/");
		tempFolder.mkdirs();
		File outFas=new File(tempFolder+"/pop_"+id+"_"+suffix+moveStr+".fas");
		File Log=new File(tempFolder+"/pop_"+id+"_"+suffix+moveStr+".log");

		CalculateTraits ct;
		ct=simulate(repRate,moveMR,alignment,outFas,fitness,directed);
		log(ct,repRate,moveMR,outFas,Log);


	}
	
	public ReplicationRateInference(File REPINalg,File outFolder,int id){
		this.id=id;
		
		//this.trials=trials;
		this.outFolder=outFolder;
		this.outFolder.mkdir();
		alignment=REPINalg;


	}

	private void log(CalculateTraits ct,double repRate,double moveMR,File outFas,File logFile){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(logFile));
			String moveMRLog=move?moveMR+"\t":"";
			bw.write(repRate+"\t"+moveMRLog+ct.getAvgPairwiseDist()+"\t"+ct.getAvgNearestNeighborDist()+"\t"+outFas+"\t"+ct.alignment.size()+"\n");
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private CalculateTraits simulate(double repRate,double moveMR,File alignment,File outFas,boolean fitness,boolean directed){
		SimulateSequenceEvolution sse=new SimulateSequenceEvolution(alignment,fitness);
		sse.setDirectedMutation(directed);
		sse.setRecombine(recombine);
		sse.setReplicationRate(repRate);
		sse.setAllReplicate(allReplicate);
		if(move){
			sse.setMovingMutationRate(moveMR);
		}else{
			sse.setMovingRate(0);
		}
		sse.simulateEvolution();
		ArrayList<Fasta> fas=sse.getPopulation();
		Fasta.write(fas, outFas);
		CalculateTraits ct=new CalculateTraits(fas);
		return ct;
	}


}
