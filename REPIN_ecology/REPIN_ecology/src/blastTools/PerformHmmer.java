package blastTools;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

public class PerformHmmer {


		public static boolean hmmer(double eValue,File output,File input,File database,boolean verbose,String outtype){
			try{
				String hmmercom="hmmscan --"+outtype+" "+output+" -E "+eValue+" "+database+" "+input;
				if(verbose)System.out.println(hmmercom);
				Process p=Runtime.getRuntime().exec(hmmercom);
				if(p.waitFor()!=0){
					System.err.println(hmmercom);
					InputStream i=p.getErrorStream();
					int c=0;
					while((c=i.read())!=-1){
						System.err.print((char)c);
					}
					System.err.println("Hmmer was not successful!");
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
		
}
