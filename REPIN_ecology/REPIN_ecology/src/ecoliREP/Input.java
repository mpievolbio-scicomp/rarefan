package ecoliREP;

import java.io.File;




public class Input {
		public File blastoutREP;
		public String homologue;
		public File genome;
		public File genbank;
		public String name;
		public Input(File BlastoutREP,File Genome,File Genbank,String Homologue,String Name){
			blastoutREP=BlastoutREP;
			homologue=Homologue;
			genome=Genome;
			genbank=Genbank;
			name=Name;
		}
		
	
	
}
