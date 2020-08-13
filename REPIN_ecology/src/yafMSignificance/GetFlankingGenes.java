package yafMSignificance;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

import overrepresented.GetFrequencyBelowPvalue;

import blastTools.ReadBlast;

import ecoliREP.BlastREPs;
import ecoliREP.Input;
import extragenicSequenceSimulation.GenerateExtragenicSequences;

import util.*;
import util.IntervalTree.Flank;
//this program is supposed to output a fasta file that includes all genes that flank the occurrences of words in extragenic space specified
//a genbank-file, fasta file and a list of words (equal sizes) are needed as input
//input file has the following structure: name\tGenomeFasta\tGBKfile
public class GetFlankingGenes {
	public static void main(String args[]) {
		File genbank = new File(args[0]);
		HashMap<String, StringBuilder> genomeFasta = ReadFasta
				.readFasta(new File(args[1]));
		String genome = genomeFasta.values().toArray(new StringBuilder[0])[0]
				.toString();
		String name = args[2];
		File out = new File(args[3]);
		File inputFile = new File(args[4]);
		int length = Integer.parseInt(args[5]);
		File wordFreqFolder = new File(args[6]);
		double threshold = Double.parseDouble(args[7]);
		int minWordNumber = Integer.parseInt(args[8]);
		ArrayList<Info> intervals = getExSpacesAboveWordNumber(minWordNumber,
				genbank, genome, new File(wordFreqFolder + "/" + name + ".out"),length);
		ArrayList<Info> genes = getFlankingGenes(intervals, genbank, name);
		File flankingGenes = new File(out + "/" + name + "_flankedGenes.out");
		Fasta.write(Fasta.makeFasta(genes, genome, true), flankingGenes);
		BlastREPs br = new BlastREPs(out, true);
		br.createBlast(inputFile, "tblastn", 1e-30, flankingGenes,false,false);
		ArrayList<String> names = br.getNames();
		HashMap<String, Input> inputInfo = br.inputInfo;
		HashMap<String, ArrayList<String>> statsHM = makeHashMap(genes);

		for (int i = 0; i < names.size(); i++) {
			System.out.println(names.get(i));
			// find all homologues genes, intervals
			File wordFreqs = new File(wordFreqFolder + "/" + names.get(i)
					+ ".out");
			Input in = inputInfo.get(names.get(i));
			HashMap<String, StringBuilder> genomeFastaLocal = ReadFasta
					.readFasta(in.genome);
			String genomeLocal = genomeFastaLocal.values().toArray(
					new StringBuilder[0])[0].toString();
			InfoTree geneTree = makeGeneTree(in.genbank);
			ArrayList<Info> blastIntervals = blastToIntervals(in.blastoutREP);
			ArrayList<Info> geneOverlaps = getOverlappingGenes(geneTree,
					blastIntervals,true);
			System.out.println(geneOverlaps.size());
			// get all exSpaces flanking homologues genes
			GenerateExtragenicSequences ge = new GenerateExtragenicSequences(
					genomeLocal, in.genbank, true);
			ArrayList<Info> extragenic = ge.getIntervals();
			InfoTree extraTree = intervalsToTree(extragenic);
			adjustGeneMap(extraTree, geneOverlaps, genomeLocal, threshold,
					wordFreqs, length, statsHM, names.get(i));
			writeStatsMap(statsHM, new File(out + "/flankedGenes.out"),
					new File(out + "/flankedGenesDetails.out"));
		}
		writeStatsMap(statsHM, new File(out + "/flankedGenes.out"), new File(
				out + "/flankedGenesDetails.out"));
	}

	public static ArrayList<Info> getExSpacesAboveWordNumber(int number,File genbank,String genome,File wordFrequencies,int length){
		GenerateExtragenicSequences ge=new GenerateExtragenicSequences(genome,genbank,true);
		ArrayList<Info> exSpaces=ge.getIntervals();
		GetFrequencyBelowPvalue gfbp=new GetFrequencyBelowPvalue(wordFrequencies,length);
		ArrayList<Info> newExSpaces=new ArrayList<Info>();
		for(int i=0;i<exSpaces.size();i++){
			if(exSpaces.get(i).getEnd()-exSpaces.get(i).getStart()>length){
				String space=Fasta.getSequence(exSpaces.get(i), genome, false);
				if(checkSpace(space, number, gfbp, length)!=null){
					newExSpaces.add(exSpaces.get(i));
				}
			}
		}
		return newExSpaces;
	}

	public static void writeStatsMap(
			HashMap<String, ArrayList<String>> statsHM, File out, File details) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			BufferedWriter bw2 = new BufferedWriter(new FileWriter(details));
			Iterator<Entry<String, ArrayList<String>>> it = statsHM.entrySet()
					.iterator();
			while (it.hasNext()) {
				Entry<String, ArrayList<String>> e = it.next();
				bw.write(e.getKey() + "\t" + e.getValue().size() + "\n");
				for (int i = 0; i < e.getValue().size(); i++) {
					bw2.write(e.getKey() + "\t" + e.getValue().get(i) + "\n");
				}
			}
			bw.close();
			bw2.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// deletes genes that are only flanked by 16mers that are below the
	// threshold and increases the counter for those that are
	// flanked on at least one side by a 16mer that is above the threshold
	public static void adjustGeneMap(InfoTree extras,
			ArrayList<Info> geneOverlaps, String genome, double threshold,
			File wordFreqs, int length,
			HashMap<String, ArrayList<String>> statsHM, String name) {
		GetFrequencyBelowPvalue gfbp = new GetFrequencyBelowPvalue(wordFreqs,length);
		for (int i = 0; i < geneOverlaps.size(); i++) {
			Info interval = geneOverlaps.get(i);
			String statistics = null;
			if (statsHM.containsKey(interval.info)
					&& !overlaps(extras, interval, false)) {
				Flank<Info> f = extras.findFlank(interval);
				String seq1 = Fasta.getSequence(f.prev, genome, false);
				String seq2 = Fasta.getSequence(f.next, genome, false);
				if ((statistics = checkSpace(seq1, threshold, gfbp, length)) != null
						|| (statistics = checkSpace(seq2, threshold, gfbp,
								length)) != null) {
					statsHM.get(interval.info).add(
							name + "\t" + statistics + "\t"
									+ interval.getStart() + "\t"
									+ interval.getEnd());
				} else {
					//statsHM.remove(interval.info);
				}
			}
		}
	}

	public static HashMap<String, ArrayList<String>> makeHashMap(
			ArrayList<Info> genes) {
		HashMap<String, ArrayList<String>> hm = new HashMap<String, ArrayList<String>>();
		for (int i = 0; i < genes.size(); i++) {
			String newName = genes.get(i).info + "_" + genes.get(i).getStart()
					+ "_" + genes.get(i).getEnd();
			hm.put(newName, new ArrayList<String>());
		}
		return hm;
	}

	public static void increaseGeneCounterAboveThreshold(
			ArrayList<Info> genesAbove, HashMap<String, Integer> Genes) {
		for (int i = 0; i < genesAbove.size(); i++) {
			if (Genes.containsKey(genesAbove.get(i).info)) {
				Genes.put(genesAbove.get(i).info, Genes
						.get(genesAbove.get(i).info) + 1);
			}
		}
	}

	public static void deleteGenesBelowThreshold(ArrayList<Info> genesBelow,
			HashMap<String, Integer> Genes) {
		for (int i = 0; i < genesBelow.size(); i++) {
			if (Genes.containsKey(genesBelow.get(i).info)) {
				Genes.remove(genesBelow.get(i).info);
			}
		}
	}

	public static String checkSpace(String space, int threshold,
			GetFrequencyBelowPvalue gfbp, int length) {
		ArrayList<String> words = chopWords(space, length);
		String stats = null;
		String bestWord = "";
		int number = 0;
		double pVal = 1;
		for (int i = 0; i < words.size(); i++) {
			if (gfbp.getFreq(words.get(i)) > threshold) {
				if (number < gfbp.getFreq(words.get(i))) {
					number = gfbp.getFreq(words.get(i));
					pVal = gfbp.getpvalue(words.get(i));
					bestWord = words.get(i);
					stats = bestWord + "\t" + number + "\t" + pVal;
				}
			}
		}

		return stats;
	}

	public static String checkSpace(String space, double threshold,
			GetFrequencyBelowPvalue gfbp, int length) {
		ArrayList<String> words = chopWords(space, length);
		String stats = null;
		String bestWord = "";
		int number = 0;
		double pVal = 1;
		for (int i = 0; i < words.size(); i++) {
			if (gfbp.getpvalue(words.get(i)) < threshold) {
				if (number < gfbp.getFreq(words.get(i))) {
					number = gfbp.getFreq(words.get(i));
					pVal = gfbp.getpvalue(words.get(i));
					bestWord = words.get(i);
					stats = bestWord + "\t" + number + "\t" + pVal;
				}
			}
		}

		return stats;
	}

	public static ArrayList<String> chopWords(String space, int length) {
		ArrayList<String> words = new ArrayList<String>();
		for (int i = 0; i < space.length() - 16; i++) {
			words.add(space.substring(i, i + 16));
		}
		return words;
	}

	public static ArrayList<Info> blastToIntervals(File blast) {
		ArrayList<Info> intervals = new ArrayList<Info>();
		ReadBlast rb = new ReadBlast(blast);
		HashMap<String, Boolean> queries = new HashMap<String, Boolean>();
		for (int i = 0; i < rb.getDatabase().size(); i++) {
			if (!queries.containsKey(rb.getQuery().get(i))) {
				queries.put(rb.getQuery().get(i), true);
				int start = rb.getStartDB().get(i) > rb.getEndDB().get(i) ? rb
						.getEndDB().get(i) : rb.getStartDB().get(i);
				int end = rb.getStartDB().get(i) < rb.getEndDB().get(i) ? rb
						.getEndDB().get(i) : rb.getStartDB().get(i);
				intervals.add(new Info(start, end, rb.getQuery().get(i)));
			}
		}
		return intervals;
	}


	
	public static ArrayList<Info> getOverlappingGenes(InfoTree geneTree,
			ArrayList<Info> blastHits,boolean keepInfo) {
		ArrayList<Info> genes = new ArrayList<Info>();
		InfoTree foundIntervals = new InfoTree();
		for (int i = 0; i < blastHits.size(); i++) {
			ArrayList<Info> temp = new ArrayList<Info>();
			foundIntervals.search(blastHits.get(i), temp);
			// if interval already found skip
			if (temp.size() > 0) {
				continue;
			}
			ArrayList<Info> al = new ArrayList<Info>();
			geneTree.search(blastHits.get(i), al);
			if (al.size() != 1) {
				System.err.println("Found " + al.size()
						+ " overlapping CDS! Ignored! "
						+ blastHits.get(i).getStart() + " "
						+ blastHits.get(i).getEnd());
			} else {
				foundIntervals.insert(blastHits.get(i));

				// keep the gene information from the reference genome..i e blast query string
				if(keepInfo)al.get(0).info = blastHits.get(i).info;
				else {
					String complement=al.get(0).info.contains("complement")?"_complement":"";
					al.get(0).info=al.get(0).getGeneOrLocusInfo()+complement+"|"+blastHits.get(i).info;
				}
				genes.add(al.get(0));
			}
		}
		return genes;
	}

	public static InfoTree makeGeneTree(File genbank) {
		ReadGenbank rgb = new ReadGenbank(genbank);
		ArrayList<Info> genes = rgb.getIntervals("CDS");
		InfoTree geneTree = new InfoTree();
		for (int i = 0; i < genes.size(); i++) {
			geneTree.insert(genes.get(i));
		}
		return geneTree;
	}

	public static InfoTree intervalsToTree(ArrayList<Info> intervals) {
		InfoTree intervalTree = new InfoTree();
		for (int i = 0; i < intervals.size(); i++) {
			intervalTree.insert(intervals.get(i));
		}
		return intervalTree;
	}

	public static ArrayList<Info> getFlankingGenes(ArrayList<Info> intervals,
			File genbank, String name) {
		InfoTree geneTree = makeGeneTree(genbank);
		return getFlankingGenes(geneTree, intervals, name);
	}
	public static ArrayList<Info> getFlankingGenesBothSides(ArrayList<Info> intervals,
			File genbank, String name) {
		InfoTree geneTree = makeGeneTree(genbank);
		return getFlankingGenesBothSides(geneTree, intervals, name);
	}
	
	public static ArrayList<Info> getFlankingGenesLeftSide(ArrayList<Info> intervals,
			File genbank, String name) {
		InfoTree geneTree = makeGeneTree(genbank);
		return getFlankingGenesLeftSide(geneTree, intervals, name);
	}
	public static ArrayList<Info> getFlankingGenesRightSide(ArrayList<Info> intervals,
			File genbank, String name) {
		InfoTree geneTree = makeGeneTree(genbank);
		return getFlankingGenesRightSide(geneTree, intervals, name);
	}
	
	public static ArrayList<Info> getFlankingGenesRightSide(InfoTree geneTree,ArrayList<Info> intervals, String name) {
		ArrayList<Info> righties=new ArrayList<Info>();
		
		for (int i = 0; i < intervals.size(); i++) {
			Info right=geneTree.successor(intervals.get(i));
			String newName = right.getGeneOrLocusInfo();
			right.info = right.info.endsWith("complement") ? name
					+ "_" + newName + "_" + "complement" : name + "_"
					+ newName;
			righties.add(right);
		}
		
		return righties;
	}
	public static ArrayList<Info> getFlankingGenesLeftSide(InfoTree geneTree,ArrayList<Info> intervals, String name) {
		ArrayList<Info> lefties=new ArrayList<Info>();
		
		for (int i = 0; i < intervals.size(); i++) {
			Info left=geneTree.predecessor(intervals.get(i));
			String newName = left.getGeneOrLocusInfo();
			left.info = left.info.endsWith("complement") ? name
					+ "_" + newName + "_" + "complement" : name + "_"
					+ newName;
			lefties.add(left);
		}
		
		return lefties;
	}
	public static ArrayList<Info> getFlankingGenesBothSides(InfoTree geneTree,ArrayList<Info> intervals, String name) {
		InfoTree flankingTreePrev = new InfoTree();
		InfoTree flankingTreeNext = new InfoTree();
		for (int i = 0; i < intervals.size(); i++) {
			Info interval = intervals.get(i);
			if (!overlaps(geneTree, interval, false)) {
				Flank<Info> f = geneTree.findFlank(interval);
				if (!overlaps(flankingTreePrev, f.prev, false)) {
					if(!f.prev.info.contains(name)){
						String newName = f.prev.getGeneOrLocusInfo();
						f.prev.info = f.prev.info.endsWith("complement") ? name
								+ "_" + newName + "_" + "complement" : name + "_"
								+ newName;
					}
					flankingTreePrev.insert(f.prev);

				}
				if (!overlaps(flankingTreeNext, f.next, false)) {
					if(!f.next.info.contains(name)){
						String newName = f.next.getGeneOrLocusInfo();
						f.next.info = f.next.info.endsWith("complement") ? name
								+ "_" + newName + "_" + "complement" : name + "_"
								+ newName;
					}
					flankingTreeNext.insert(f.next);
				}

			}
		}
		return getIntersection(flankingTreeNext,flankingTreePrev);
	}
	
	public static ArrayList<Info> getIntersection(InfoTree next,InfoTree prev){
		ArrayList<Info> intersection=new ArrayList<Info>();
		ArrayList<Info> nextList=next.parseTree();
		for(int i=0;i<nextList.size();i++){
			
			if(prev.find(nextList.get(i))!=null){
				intersection.add(nextList.get(i));
			}
		}
		return intersection;
	}
	
	
	public static ArrayList<Info> getFlankingGenes(InfoTree geneTree,
			ArrayList<Info> intervals, String name) {
		InfoTree flankingTree = new InfoTree();
		for (int i = 0; i < intervals.size(); i++) {
			Info interval = intervals.get(i);
			if (!overlaps(geneTree, interval, false)) {
				Flank<Info> f = geneTree.findFlank(interval);
				if (!overlaps(flankingTree, f.prev, false)) {
					String newName = getGeneOrLocusTag(f.prev.info);
					f.prev.info = f.prev.info.endsWith("complement") ? name
							+ "_" + newName + "_" + "complement" : name + "_"
							+ newName;
					flankingTree.insert(f.prev);
				}
				if (!overlaps(flankingTree, f.next, false)) {
					String newName = getGeneOrLocusTag(f.next.info);
					f.next.info = f.next.info.endsWith("complement") ? name
							+ "_" + newName + "_" + "complement" : name + "_"
							+ newName;
					flankingTree.insert(f.next);
				}

			}
		}
		return flankingTree.parseTree();
	}

	public static String getGeneOrLocusTag(String info) {
		String[] split = info.split("\\s+");
		return split[getGeneOrLocusTagIndex(split)];
	}

	public static int getGeneOrLocusTagIndex(String[] split) {
		for (int i = 0; i < split.length; i++) {
			if (split[i].contains("gene=")) {
				return i + 1;
			}
		}
		for (int i = 0; i < split.length; i++) {
			if (split[i].contains("locus_tag=")) {
				return i + 1;
			}
		}
		return 0;
	}

	public static boolean overlaps(InfoTree tree, Info x, boolean show) {
		ArrayList<Info> al = new ArrayList<Info>();
		tree.search(x, al);
		if (show) {
			for (int i = 0; i < al.size(); i++) {
				System.out.println(al.get(i).info + " " + al.get(i).getStart()
						+ " " + al.get(i).getEnd());
			}
		}
		return al.size() > 0;
	}

	public static ArrayList<Info> find(String genome, ArrayList<String> words) {
		ArrayList<Info> infos = new ArrayList<Info>();
		for (int i = 0; i < words.size(); i++) {
			String word = words.get(i);
			infos.addAll(find(genome, word));

		}
		return infos;
	}
	
	public static ArrayList<Info> find(String genome,String word){
		ArrayList<Info> all=new ArrayList<Info>();
		String inverted = DNAmanipulations.reverse(word);
		all.addAll(findBase(genome.toUpperCase(),word.toUpperCase()));
		all.addAll(findBase(genome.toUpperCase(),inverted.toUpperCase()));
		return all;
	}
	
	private static ArrayList<Info> findBase(String genome, String word) {
		int j = 0;
		ArrayList<Info> infos = new ArrayList<Info>();
		while (j != -1) {
			j = genome.indexOf(word, j);
			if (j == -1) {
				return infos;
			}
			j++;
			infos.add(new Info(j, j + word.length(), ""));
		}
		return infos;
	}
	
	

	public static ArrayList<Info> posToInfo(ArrayList<SequencePositions> sq,
			int size) {
		ArrayList<Info> infos = new ArrayList<Info>();
		for (int i = 0; i < sq.size(); i++) {
			if (sq.get(i).position > 0) {
				infos.add(new Info(sq.get(i).position, sq.get(i).position
						+ size / 2, ""));
				System.out.println(new Info(sq.get(i).position,
						sq.get(i).position + size / 2, ""));
			} else {
				infos.add(new Info(sq.get(i).position, sq.get(i).position
						+ size / 2, "inverted"));
				System.out.println(new Info(sq.get(i).position,
						sq.get(i).position + size / 2, "inverted"));
			}

		}
		return infos;
	}

}
