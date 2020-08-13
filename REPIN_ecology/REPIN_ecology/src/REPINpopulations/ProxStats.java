package REPINpopulations;

public class ProxStats {
	String name;
	int numRAYTs;
	int numDifferentClusters;
	boolean containsNoCluster;
	int focalSeedGroup;
	int largestClusters;
	public ProxStats(String name,int focalSeedGroup,int numRAYTs,int numDifferentClusters,int numLC) {
		this.name=name;
		this.numRAYTs=numRAYTs;
		this.focalSeedGroup=focalSeedGroup;
		this.numDifferentClusters=numDifferentClusters;
		this.largestClusters=numLC;
	}
	public static String printHeading() {
		return "Name\tfocalSeedGroup\tnumRAYTs\tnumClusters\tnumLargestClusters\n";
	}
	public String print() {
		return name+"\t"+focalSeedGroup+"\t"+numRAYTs+"\t"+numDifferentClusters+"\t"+largestClusters+"\n";
	}
}
