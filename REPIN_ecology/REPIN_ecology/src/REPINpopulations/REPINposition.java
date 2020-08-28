package REPINpopulations;

public class REPINposition {
	public int start;
	public int end;
	public int id;
	public REPINposition(int start,int end,int id/*sequence position in genome fasta file*/) {
		this.start=start;
		this.end=end;
		this.id=id;
	}
}
