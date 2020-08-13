package util.phylogenetics;

public class Score{
	double score;
	int windowLength;
	int numberSites;
	String sitePattern;
	public Score(double s,int win,int num,String pattern){
		score=s;
		windowLength=win;
		numberSites=num;
		sitePattern=pattern;
	}
	public String toString(){
		return score+"\t"+windowLength+"\t"+numberSites+"\t"+sitePattern;
	}
}
