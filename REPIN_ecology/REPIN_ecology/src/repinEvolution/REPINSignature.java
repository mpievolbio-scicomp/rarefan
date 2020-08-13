package repinEvolution;

public class REPINSignature {
	String genomeID;
	int REPINType;
	int isREPIN;
	String all;
	public REPINSignature(String genomeID,int REPINtype,int isREPIN) {
		this.genomeID=genomeID;
		this.REPINType=REPINtype;
		this.isREPIN=isREPIN;
		setAll();
	}
	public REPINSignature(REPIN r) {
		this.genomeID=r.genomeID;
		this.REPINType=r.REPINtype;
		this.isREPIN=r.isREPIN;
		setAll();
	}
	
	private void setAll() {
		all=this.genomeID+"\t"+this.REPINType+"\t"+this.isREPIN;

	}
	public int compareTo(REPINSignature rs) {
		return all.compareTo(rs.all);
	}
	
	public String toString() {
		return all;
	}

	
	
	
}
