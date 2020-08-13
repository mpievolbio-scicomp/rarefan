package repinEvolution;

import java.util.*;

import frequencies.*;

public class REPINpopulations {
		HashMap<String,REPINProperties> bacterium=new HashMap<String, REPINProperties>();
		public REPINpopulations() {
			bacterium=new HashMap<String, REPINProperties>();
		}
		public void addREPINPopulation(String REPIN,REPINProperties repinProps) {
			bacterium.put(REPIN, repinProps);
		}
		public REPINProperties getREPINProperties(String REPIN) {
			return bacterium.get(REPIN);
		}
		public String[] getAllREPINs() {
			return bacterium.keySet().toArray(new String[0]);
		}
}
