package frequencies;

import java.io.*;
import java.util.*;


public class CalculateAmplificationRates {
	File inFolder;
	File scriptPath=new File("/usr/local/bin/Rscript");
	double genomeMutationRate=8.9e-11;//for E. coli from Wielgoss et al. "Mutation Rate Inferred From Synonymous Substitutions in a Long-Term Evolution Experiment With Escherichia coli"
	HashMap<Double/*mutationRate*/,HashMap<String/*species*/,/*fitness*/ArrayList<Double>>> similarityData=new HashMap<Double, HashMap<String/*species*/,ArrayList<Double>>>();

	HashMap<Double/*mutationRate*/,HashMap<String/*species*/,/*fitness*/ArrayList<Double>>> fitnessData=new HashMap<Double, HashMap<String/*species*/,ArrayList<Double>>>();
	HashMap<Double/*mutationRate*/,HashMap<String,ArrayList<Double>>/*amplificationRate*/> amplificationRates=new HashMap<Double, HashMap<String,ArrayList<Double>>>();
	HashMap<Double,HashMap<String,Double>> distances=new HashMap<Double, HashMap<String,Double>>();
	
	public CalculateAmplificationRates(File inFolder){
		this.inFolder=inFolder;
		readData(fitnessData,"quasispecies_fitness_");
		readData(similarityData,"minDiff_");
		calculateAmpRates();
		
	}
	
	private double calculateAmpRate(double fit,double refMutRate){
		double ampRate=fit/(refMutRate/genomeMutationRate);
		return ampRate;
	}
	
	
	private void calculateAmpRates(){
		Double[] rates=fitnessData.keySet().toArray(new Double[0]);
		for(int i=0;i<rates.length;i++){
			double refMutRate=rates[i];
			amplificationRates.put(rates[i],new HashMap<String, ArrayList<Double>>());
			String[] species=fitnessData.get(rates[i]).keySet().toArray(new String[0]);
			for(int j=0;j<species.length;j++){
				ArrayList<Double> amp=new ArrayList<Double>();
				ArrayList<Double> fit=fitnessData.get(rates[i]).get(species[j]);
				for(int k=0;k<fit.size();k++){
					amp.add(calculateAmpRate(fit.get(k),refMutRate));
				}
				
				amplificationRates.get(rates[i]).put(species[j],amp );
			}
			
		}
	}
	
	private void readData(HashMap<Double/*mutationRate*/,HashMap<String/*species*/,/*fitness*/ArrayList<Double>>> data,String startsWith){
		File list[]=inFolder.listFiles();
		for(int i=0;i<list.length;i++){
			String name=list[i].getName();
			if(name.startsWith(startsWith)){
				String split[]=name.substring(0, name.length()-4).split("_");
				double mr=Double.parseDouble(split[split.length-1]);
				data.put(mr, getData(list[i]));
			}
		}
	}
	
	private HashMap<String,ArrayList<Double>> getData(File in){
		HashMap<String,ArrayList<Double>> data=new HashMap<String, ArrayList<Double>>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(in));
			String line="";
			while((line=br.readLine())!=null){
				line=line.replace("}", "");
				line=line.replace("{", "");
				line=line.replace(",", "");
				String split[]=line.split("\\s+");
				String species=split[split.length-1];
				ArrayList<Double> fit=new ArrayList<Double>();
				for(int i=0;i<3;i++){
					fit.add(Double.parseDouble(split[i]));
				}
				data.put(species, fit);
			}
			
			br.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return data;
	}
	

	

	
	private  File  write(String name,HashMap<Double,HashMap<String,ArrayList<Double>>> results){
		File out=new File(inFolder+"/"+name+".txt");
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			Double[] rates=results.keySet().toArray(new Double[0]);
			bw.write("baseMutRate\tspecies");
			
			for(int i=0;i<rates.length;i++){
				String[] species=results.get(rates[i]).keySet().toArray(new String[0]);
				for(int j=0;j<species.length;j++){
					ArrayList<Double> data=results.get(rates[i]).get(species[j]);

					//annotation
					if(i==0&&j==0){

						for(int k=0;k<data.size();k++){
							bw.write("\t"+k);
						}
						bw.write("\n");

					}

					
					bw.write(rates[i]+"\t"+species[j]);

					
					
					//write Data
					for(int k=0;k<data.size();k++){
						
						bw.write("\t"+data.get(k));
					}
					bw.write("\n");
				}
			}
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		return out;
	}

//	private void plot(File in){
//		File out=new File(in+".R");
//		File pdf=new File(in+".pdf");
//		RCode rc=new RCode();
//		rc.addRCode("library(ggplot2)");
//		rc.addRCode("t<-read.table(\""+in+"\",header=TRUE)");
//		rc.addRCode("df=as.data.frame(t)");
//		rc.addRCode("ggplot(df,aes(x=baseMutRate,y=X0,group=species))+geom_line(aes(color=factor(species)))");
//		rc.addRCode("ggsave(\""+pdf+"\")");
//		R_functions.runRCode(rc, scriptPath);
//		R_functions.writeRCode(rc, out);
//	}
	
	public void writeAll(){
		File inf=write("inferredAmplificationRates",amplificationRates);
		//plot(inf);
		File minDist=write("minDistances",similarityData);
		//plot(minDist);
		File fitDat=write("fitnessData",fitnessData);
		//plot(fitDat);
	}
}
