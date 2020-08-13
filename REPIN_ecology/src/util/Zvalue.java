package util;


public class Zvalue {
	
	public static void main(String[] args){
		System.out.println(getZvalue(new NP(19,0.68,"parallel"),new NP(1079,0.409,"others")));
	}
	
	public static class NP{
		double p;
		int n;
		String name;
		public NP(int n,double p,String name){
			this.n=n;
			this.p=p;
			this.name=name;
		}
		public String toString(){
			return name+" "+n+" "+p;
		}
	}

	public static double getZvalue(NP np1,NP np2){
		double pb=(np1.n*np1.p+np2.n*np2.p)/(np1.n+np2.n);
		double z=(np1.p-np2.p)/Math.pow(pb*(1.0-pb)*(1.0/np1.n+1.0/np2.n),0.5);
		return z;
	}
}
