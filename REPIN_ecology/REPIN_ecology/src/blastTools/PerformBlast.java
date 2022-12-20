package blastTools;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
class StreamGobbler extends Thread
{
    InputStream is;
    String type;
    
    StreamGobbler(InputStream is, String type)
    {
        this.is = is;
        this.type = type;
    }
    
    public void run()
    {
        try
        {
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line=null;
            while ( (line = br.readLine()) != null)
                System.out.println(type + ">" + line);    
            } catch (IOException ioe)
              {
                ioe.printStackTrace();  
              }
    }
}

public class PerformBlast {
	public static boolean blast(String blastpath,String builderpath,String program,double eValue,File output,File input,File database,boolean verbose,boolean createDB,boolean DNA,boolean Filter){
		try{
			File db=new File(database+".nhr");
			String dbType=DNA?"F":"T";
			if(!db.exists()||createDB){
				String blastcom=builderpath+" -p "+dbType+" -i "+database;

				Process p=Runtime.getRuntime().exec(blastcom);
		        System.out.println("test: "+blastcom);
	            // any error message?
	            StreamGobbler errorGobbler = new 
	                StreamGobbler(p.getErrorStream(), "ERROR");            
	            
	            // any output?
	            StreamGobbler outputGobbler = new 
	                StreamGobbler(p.getInputStream(), "OUTPUT");
	                
	            // kick them off
	            errorGobbler.start();
	            outputGobbler.start();

				if(p.waitFor()!=0){

					System.err.println(blastcom);
					InputStream i=p.getErrorStream();
					int c=0;

					while((c=i.read())!=-1){
						System.err.print((char)c);
					}
					System.err.println("Formatdb was not successful!");
					return false;
				}

			}
	        System.out.println("test: ");

			char filter=Filter?'T':'F';
			String blastcom=blastpath+" -p "+ program+" -m 8 -o "+output+" -i "+input+" -e "+eValue+" -d "+database+" -F "+filter+" -b 50000000";
			if(verbose)System.out.println(blastcom);
			Process p=Runtime.getRuntime().exec(blastcom);

			if(p.waitFor()!=0){
				System.err.println(blastcom);
				InputStream i=p.getErrorStream();
				int c=0;
				while((c=i.read())!=-1){
					System.err.print((char)c);
				}
				System.err.println("Blast was not successful!");
				return false;
			}
		}catch(InterruptedException e){
			e.printStackTrace();
			return false;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
		return true;
		
	}
	public static void deleteDatabases(File db){
		File d1=new File(db+".nhr");
		File d2=new File(db+".nin");
		File d3=new File(db+".nsq");
		if(d1.exists())d1.delete();
		if(d2.exists())d2.delete();
		if(d3.exists())d3.delete();

	}
}
