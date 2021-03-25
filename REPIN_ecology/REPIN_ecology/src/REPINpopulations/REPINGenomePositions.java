package REPINpopulations;

import java.util.ArrayList;
import java.util.HashMap;

public class REPINGenomePositions {
	public HashMap<String,ArrayList<REPINposition>> repinGenomePositions=new HashMap<String, ArrayList<REPINposition>>();
	public REPINGenomePositions(HashMap<String,ArrayList<REPINposition>> repinGenomePositions) {
		this.repinGenomePositions=repinGenomePositions;
	}
	public HashMap<String,ArrayList<REPINposition>> get(){
		return repinGenomePositions;
	}
}
