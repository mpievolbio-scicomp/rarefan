package repinEvolution;


public class REPINTransition {
	
	//input is a single cluster
	//calculate transition: genomeID1(ref) genomeID2 REPINtype1 REPINType2 isREPIN1 isREPIN2
	REPINSignature r1;
	REPINSignature r2;
	public REPINTransition(REPINSignature r1,REPINSignature r2,String reference) {
		
		if(r1.genomeID.equals(reference)) {
			this.r1=r1;
			this.r2=r2;
		}else {
			this.r1=r2;
			this.r2=r1;
		}
	}
	
	public boolean equals(Object rt) {
		return(rt.toString().equals(this.toString()));
	}
	public String toString() {
		return r1.toString()+"\t"+r2.toString();
	}
	public int compareTo(REPINTransition rt) {
		
		return this.toString().compareTo(rt.toString());
	
	}
	
	public int hashCode() {
		return (this.toString()).hashCode();
	}
	

	
}
