package pairwiseAlignment;

public class PairwiseAlignment {
	String alA;
	String alB;
	public PairwiseAlignment(String a,String b){
		alA=a;
		alB=b;
	}
	public String toString(){
		return alA+"\n"+alB+"\n";
	}
	public String getSeq1(){
		return alA;
	}
	public String getSeq2(){
		return alB;
	}
	
	public void deleteGapPositions(){
		StringBuffer sbA=new StringBuffer();
		StringBuffer sbB=new StringBuffer();
		for(int i=0;i<alA.length();i++){
			if(alA.charAt(i)!='-'&&alB.charAt(i)!='-'){
				sbA.append(alA.charAt(i));
				sbB.append(alB.charAt(i));
			}
		}
		alA=sbA.toString();
		alB=sbB.toString();
	}
	
	public void deleteGapPositionsFirst(){
		StringBuffer sbA=new StringBuffer();
		StringBuffer sbB=new StringBuffer();
		for(int i=0;i<alA.length();i++){
			if(alA.charAt(i)!='-'){
				sbA.append(alA.charAt(i));
				sbB.append(alB.charAt(i));
			}
		}
		alA=sbA.toString();
		alB=sbB.toString();
	}
}
