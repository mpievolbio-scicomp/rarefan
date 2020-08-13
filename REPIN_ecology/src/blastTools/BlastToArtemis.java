package blastTools;

import java.io.File;

import util.WriteArtemis;

//creates an artemis output of a blast tabular file
//input blast tab
//input artemis output 
//output artemis

public class BlastToArtemis {
	public static void main(String args[]){
		File blast=new File(args[0]);
		File artout=new File(args[1]);
		ReadBlast rb=new ReadBlast(blast);
		WriteArtemis.write(rb.getStartDB(),rb.getEndDB(),rb.getQuery(),artout); 
	}

}
