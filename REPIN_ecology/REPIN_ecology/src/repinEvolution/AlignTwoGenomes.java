
	package repinEvolution;

	import java.io.*;
	import java.util.*;

	import util.*;
	import util.phylogenetics.RunTreePrograms;

	public class AlignTwoGenomes {
		
		
		HashMap<String,int[]> positionMap=new HashMap<String,int[]>();
		File alignmentOut;
//		HashMap<String,ArrayList> alignment=new HashMap<String, ArrayList>();
		//HashMap<String, ArrayList> alignmentC=new HashMap<String, ArrayList>();
		HashMap<String,String> genomeSequences=new HashMap<String, String>();
		String[] ids;
		File[] genomes;
		public static void main(String[] args) {
			File G1=new File("/Users/bertels/Documents/Arne/mutationFrequencies/chlororaphis/mauve/chlTAMOak81.fas");
			File G2=new File("/Users/bertels/Documents/Arne/mutationFrequencies/chlororaphis/mauve/chl50083.fas");
			String reference=getID(G1);
			AlignTwoGenomes atg=new AlignTwoGenomes(G1, G2);
		}
		
		public Integer getPositionQuery(int pos,String ref){
			return positionMap.get(ref)[pos];
		}
		
		public AlignTwoGenomes(File G1,File G2) {
			genomes= new File[]{G1,G2};
			ids=getIDs();
			HashMap<String,ArrayList> alignment=alignGenomes(G1,G2);

					setGenomeSequences();
			
			makePositionMap(alignment);
		}
		
		private void setGenomeSequences() {
			for(int i=0;i<ids.length;i++) {
				genomeSequences.put(ids[i], Fasta.readFasta(genomes[i]).get(0).getSequence());
			}
		}
		
		public String getGenomeSequence(String id) {
			return genomeSequences.get(id);
		}
		
		private void makePositionMap(HashMap<String,ArrayList> alignment) {
			for(int i=0;i<ids.length;i++) {
				ArrayList<Integer> ref=alignment.get(ids[i]);
				ArrayList<Integer> query=alignment.get(ids[(i+1)%2]);
				fillMap(ref,query,ids[i]);
			}
		}
		
		private void fillMap(ArrayList<Integer> ref,ArrayList<Integer> query,String id) {
			HashMap<Integer,Integer> map=new HashMap<Integer,Integer>();
			System.out.println(ref.size()+" "+query.size());
			
			for(int j=0;j<ref.size();j++) {
				if(!map.containsKey(ref.get(j))) {
					map.put(ref.get(j),query.get(j));
				}
			}
			int[] temp=convertMapToInt(map);
			positionMap.put(id,temp);
		}
		
		private int[] convertMapToInt(HashMap<Integer,Integer> map) {
			Integer[] keys=map.keySet().toArray(new Integer[0]);
			int max=getMax(keys);
			int[] result=new int[max+1];
			for(int i=0;i<max+1;i++) {
				if(map.containsKey(i)) {
					result[i]=map.get(i);
				}else {
					result[i]=-1;
				}
			}
			return result;
		}
		
		private int getMax(Integer[] array) {
			int max=-1;
			for(int i=0;i<array.length;i++) {
				if(max<array[i]) {
					max=array[i];
				}
			}
			return max;
		}
		
		
//		public HashMap<String, ArrayList> getAlignment(boolean character){
//			//if(character) {
//			//	return alignmentC;
//			//}else {
//				return alignment;
//			//}
//		}
		
		public String[] getIDs() {
			String[] ids=new String[genomes.length];
			for(int i=0;i<genomes.length;i++) {
				ids[i]=getID(genomes[i]);
			}
			return ids;
		}
		
		private HashMap<String, ArrayList> alignGenomes(File G1,File G2) {
			HashMap<String,ArrayList> alignment=new HashMap<String, ArrayList>();
			File runDir=G1.getParentFile();
			alignmentOut=new File(runDir+"/"+ids[0]+"_"+ids[1]+".xmfa");
			String call="/Applications/Mauve.app/Contents/MacOS/progressiveMauve --weight=1000000  --output="+alignmentOut+" "+G1+" "+G2;
			if(!alignmentOut.exists())RunTreePrograms.runProgram(call, "", runDir);
			alignment=readAlignment(alignmentOut);
			//print(alignment,100);
			return alignment;
		}
		
//		private void print(HashMap<String,ArrayList> alignment,int size) {
//			String[] ids=alignment.keySet().toArray(new String[0]);
//				for(int i=alignment.get(ids[0]).size()-1;i>alignment.get(ids[0]).size()-size;i--) {
//					for(int j=0;j<ids.length;j++) {
	//
//						System.out.print(" "+alignment.get(ids[j]).get(i));
//					}
//				System.out.println();
//			}
//		}
		
		//reads only the first two fasta entries
		//need to read this properly, maybe the output should be HAshMap<String,ArrayList<Integer>> then you will have a position mapping because the integer arrays should be the same length, 
		//maybe once we have them we can convert it to something else HashMap<String,HashMap<Integer,ArrayList<Integer>>> oder so
		private HashMap<String,ArrayList> readAlignment(File align) {
			HashMap<String,ArrayList> alignment=new HashMap<String, ArrayList>();
			try {
				BufferedReader br=new BufferedReader(new FileReader(align));
				String line;
				String currentSeqID="";
				int inc=0;
				int from=0;
				int k=0;
				System.out.println(align);
				while((line=br.readLine())!=null) {
					k++;
					if(k%10000==0) {
						long heapFreeSize = Runtime.getRuntime().freeMemory();
						long heapSize = Runtime.getRuntime().totalMemory(); 
						//System.out.println(heapFreeSize/1000000+" "+ heapSize/1000000 );
					}
					if(line.startsWith(">")) {
						String split[]=line.split("\\s+");
						currentSeqID=getID(new File(split[3]));
						if(!alignment.containsKey(currentSeqID)) {
							alignment.put(currentSeqID,new ArrayList<Integer>());
							//alignmentC.put(currentSeqID,new ArrayList<Character>());

						}
						String range[]=split[1].split(":|-");
						from=Integer.parseInt(range[1]);
						int to=Integer.parseInt(range[2]);
						boolean rev=false;
						if(from>to) {
							rev=true;
						}
						inc=rev?-1:1;
						
					}else if(line.startsWith("#")){
						
					}else {
						line=line.strip();
						for(int i=0;i<line.length();i++) {
							if(line.charAt(i)!='-' && line.charAt(i)!='=') {
								alignment.get(currentSeqID).add(from);
								//alignmentC.get(currentSeqID).add(line.charAt(i));

								from+=inc;
							}else if(line.charAt(i)=='=') {
								String otherSeqID=getOtherSeqID(currentSeqID,alignment);
								int currsize=alignment.get(currentSeqID).size();
								int othersize=alignment.get(otherSeqID).size();
								if(currsize>othersize) {
									int diff=currsize-othersize;
									for(int j=0;j<diff;j++) {
										alignment.get(otherSeqID).add(alignment.get(otherSeqID).get(othersize-1));
									}
								}
								currentSeqID="";
								inc=0;
								from=0;

							}else {
								alignment.get(currentSeqID).add(from);
								//alignmentC.get(currentSeqID).add(line.charAt(i));

							}
						}
					}
				}
				br.close();
			}catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
			return alignment;
		}
		
		private String getOtherSeqID(String curr,HashMap<String,ArrayList> alignment) {
			String[] keys=alignment.keySet().toArray(new String[0]);
			for(int i=0;i<keys.length;i++) {
				if(!curr.equals(keys[i])) {
					return keys[i];
				}
			}
			return curr;
		}
		
		public static String getID(File f) {
			String split[]=f.getName().split("\\.");
			String id=split[split.length-2];
			return id;
		}
		
	}

