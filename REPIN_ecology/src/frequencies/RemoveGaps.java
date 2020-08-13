package frequencies;

import java.io.File;
import java.util.ArrayList;

import util.Fasta;
import util.phylogenetics.Alignment;

public class RemoveGaps {
	public static void main(String args[]) {
		File alignment=new File(args[0]);
		File out=new File(args[1]);
		ArrayList<Fasta> fas=Fasta.readPhylip(alignment);
		changeIdents(fas);
		File fasOut=new File(alignment+".fas");

		Fasta.write(fas, fasOut);
		Alignment alg=new Alignment(Fasta.readFasta(fasOut));
		alg=removeGaps(alg);
		fas=alg.toFasta();
		Fasta.write(fas, out);
		Fasta.writePhylip(fas, new File(out+".phy"), 20);
	}
	
	private static void changeIdents(ArrayList<Fasta> fas) {
		for(int i=0;i<fas.size();i++) {
			fas.set(i, new Fasta(fas.get(i).getIdent().split("\\.")[0],fas.get(i).getSequence()));
		}
	}
	
	private static Alignment removeGaps(Alignment alg) {
		Alignment newAlg=new Alignment();
		newAlg.changeIdents(alg.getIdents());
		for(int i=0;i<alg.getLength();i++) {
			String col=alg.getColumn(i);
			if(!col.contains("-")) {
				newAlg.addColumn(col);
			}
		}
		return newAlg;
	}

}
