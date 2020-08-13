package util;

import java.io.File;

public class SGEutils {
	/**
	 * Creates a job script that can be submitted via qsub on a sun grid engine.
	 * 
	 * @param command
	 * the command that is to be executed in the job file
	 * @param jobscriptFile
	 * The name and location of the job file that is to be executed.
	 * @param mem
	 * the maximum amount of memory to be used by the job
	 */
	
	String script="#!/bin/bash\n"+
			"#$ -S /bin/bash\n"+
			"#$ -P %project%\n"+
			"#$ -o %file%\n"+
			"#$ -q %queue%\n"+
			"#$ -l mem_total=%mem%G\n"+
			"source ~/.bashrc";
	
	public void setQueue(String queue){
		
	}
	
	public void setFile(File jobFile){
		
	}
	
	public void setCommand(File jobFile){
		
	}
	
	public void setProject(File jobFile){
		
	}
	
	public void makeJobScript(String command,File jobscriptFile,int mem){
		
	}
}
