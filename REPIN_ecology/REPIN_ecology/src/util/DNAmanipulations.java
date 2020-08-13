package util;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;



public class DNAmanipulations {
	
	
	
	public static String reverse(String seq){
		StringBuilder sb=new StringBuilder();
		seq=seq.toUpperCase();
		for(int i=seq.length()-1;i>-1;i--){
			sb.append(seq.charAt(i)=='A'?'T':seq.charAt(i)=='T'?'A':seq.charAt(i)=='C'?'G':seq.charAt(i)=='G'?'C':'N');
		}
		
		return sb.toString();
	}
	
	public static char getComplement(char c){
		return c=='A'?'T':c=='T'?'A':c=='C'?'G':c=='G'?'C':'N';
	}
	
	public static String generateRandomSequence(int length,double GC){
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<length;i++){
			sb.append(getRandomBase(GC));
		}
		return sb.toString();
	}


	
	public static char getRandomBase(double GC){
		double rand=Math.random();
		double r=Math.random();
		if(rand>GC){
			if(r>0.5){
				return 'A';
			}else{
				return 'T';
			}
		}else {
			if(r>0.5){
				return 'C';

			}else{
				return 'G';
			}
		}
	}
	
	public static char getRandomBase(double transit,double transver,char base){
		double rand=Math.random();
		double r=Math.random()*(transit+transver);
		if(r<transit){
			if(base=='A'){
				return 'G';
			}else if(base=='G'){
				return 'A';
			}else if(base=='T'){
				return 'C';
			}else if(base=='C'){
				return 'T';
			}else{
				return 0;
			}
		}else {
			if(base=='A'){
				if(rand>0.5)return 'T';
				else return 'C';

			}else if(base=='T'){
				if(rand>0.5){
					return 'A';
				}else{
					return 'G';
				}
			}else if(base=='C'){
				if(rand>0.5){
					return 'A';
				}else{
					return 'G';
				}
			}else if(base=='G'){
				if(rand>0.5){
					return 'C';
				}else{
					return 'T';
				}
			}else{
				return 0;
			}
		}
	}
	
	public static String randomizeSequenceMutateToSelf(String seq,double mutationProb,double GC){
		StringBuffer sb=new StringBuffer(seq);
		int seqlength=seq.length();
		
		int numberMutations=(int)(seq.length()*mutationProb);
		for(int i=0;i<numberMutations;i++){
			int randPosition=(int)(Math.random()*seqlength);
			sb.setCharAt(randPosition, getMutationSelfMutate(seq.charAt(randPosition),GC));
			
		}
		return sb.toString();
	}
	
	/**
	 * Introduces mutations in a given sequences. Possible mutations are repeats, single basepair insertions and substitutions.
	 * @param seq
	 * DNA sequence.
	 * @param mutationProb
	 * Proportion of the sequence that will be mutated. Eg. if the sequence is 100bp long and the mutationProb is 0.1 then 10
	 * mutations will be randomly introduced. Mutations can also occur twice at the same position.
	 * @param insertionProb
	 * Same as above but with random base introductions, gaps are introduced in the ancestor.
	 * @param repeatProb
	 * Copies a piece of DNA of a randomly selected length between minRepLength and maxRepLength and randomly introduces it into the genome,
	 *  gaps are introduced in the ancestor.
	 * @param minRepLength
	 * @param maxRepLength
	 * @param GC
	 * Determines the GC content of the base introductions and substitutions.
	 * @return
	 * Returns the randomized sequence as well as the original sequence (with GAPs).
	 */
	public static ArrayList<Fasta> randomizeSequence(String seq,double mutationProb,double insertionProb,double repeatProb,int minRepLength,int maxRepLength,double GC){
		int seqlength=seq.length();
		//define mutation types
		final int subs=0;
		final int ins=1;
		final int reps=2;
		ArrayList<Integer> mutations=new ArrayList<Integer>();
		//calculate number of substitutions
		mutations.add((int)(mutationProb*seqlength));
		//calculate number of insertions
		mutations.add((int)(insertionProb*seqlength));
		//calculate number of repeats
		mutations.add((int)(repeatProb*seqlength));
		int allMutations=mutations.get(subs)+mutations.get(ins)+mutations.get(reps);
		//initialize ancestor and offspring in mutable stringbuffer
		StringBuffer offspring=new StringBuffer(seq);
		StringBuffer ancestor=new StringBuffer(seq);
		for(int i=0;i<allMutations;i++){
			seqlength=offspring.length();
			//calculate the position of the event
			int randPosition=(int)(Math.random()*seqlength);

			//randomly choose one mutation type, the more mutations there are of one type the more likely this mutation gets
			int choice=choose(mutations);
			
			//if there are no such mutations left then try again, this should never be the case...
			if(mutations.get(choice)<=0){
				System.out.println(mutations.get(choice));
				throw new RuntimeException("THERE MAY BE SOMETHING WRONG WITH THE CODE.");
				
			}
			
			//reduce the final choice by one
			mutations.set(subs,mutations.get(choice)-1);
			
			//then start to insert mutations
			if(choice==subs){
				
				//this is simple, change the base at this position, nothing to do in the ancestor
				offspring.setCharAt(randPosition, getMutation(offspring.charAt(randPosition),GC));
				
			}else if(choice==ins){
				//determine base to be inserted
				char base=getRandomBase(GC);
				//insert base
				offspring.insert(randPosition, base);
				//insert gap in ancestor
				ancestor.insert(randPosition, '-');
			}else if(choice==reps){
				
				//determine length of repeat
				int length=(int)(Math.random()*(maxRepLength-minRepLength))+minRepLength;
				
				//get repeat sequence
				CharSequence repeat;
				if(randPosition+length<offspring.length())
					repeat=offspring.subSequence(randPosition, randPosition+length);
				else
					repeat=offspring.subSequence(randPosition, offspring.length());
				//detmine target position 
				int target=(int)(Math.random()*seqlength);
				
				//insert into target
				offspring.insert(target, repeat);
				
				//insert GAPs for now...not sure what else to do TODO
				StringBuffer gaps=new StringBuffer();
				for(int j=0;j<repeat.length();j++){
					gaps.append('-');
				}
				ancestor.insert(target,gaps);
				
			}else{
				RuntimeException e=new RuntimeException("The random number generator should only choose numbers from 0 to 2. So you see this then there is something seriously wrong. Hence I will terminate.");
				throw e;
			}
			
		}
		seq=randomizeSequence(seq, mutationProb, GC);
		
		
		ArrayList<Fasta> newAlignment=new ArrayList<Fasta>();
		newAlignment.add(new Fasta("ancestor",ancestor.toString()));
		newAlignment.add(new Fasta("offspring",offspring.toString()));
		
		return newAlignment;
	}
	
	/**
	 * Method to randomly determine an index of the given ArrayList. Where the probability of choosing an index should correlate to the value of the stored integer. 
	 * All integers should be larger than -1.
	 * @param mutations
	 * @return
	 */
	private static int choose(ArrayList<Integer> mutations){
		int sum=0;
		//first determine the number of elements in the list
		for(int i=0;i<mutations.size();i++){
			if(mutations.get(i)<0){
				throw new RuntimeException("Negative values are not allowed in this array! CHECK CODE!");
			}
			sum+=mutations.get(i);
		}
		//then randomly choose a number that is smaller than the total number of elements
		int rand=(int)(Math.random()*sum);
		int part=0;
		//determine the position of this number
		for(int i=0;i<mutations.size();i++){
			part+=mutations.get(i);
			if(rand<part){
				//and return the index
				return i;
			}
			
		}
		throw new RuntimeException("Something went wrong with the code! Please make sure that this case cannot occur!");
		
	}
	
	
	public static String randomizeSequence(String seq,double mutationProb,double GC){
		StringBuffer sb=new StringBuffer(seq);
		int seqlength=seq.length();
		//calculate number of mutations 
		int numberMutations=(int)(seq.length()*mutationProb);
		for(int i=0;i<numberMutations;i++){
			//and distribute them
			int randPosition=(int)(Math.random()*seqlength);
			sb.setCharAt(randPosition, getMutation(seq.charAt(randPosition),GC));
			
		}
		return sb.toString();
	}
	
	public static String randomizeSequence(String seq,double transit,double transver,double GC){
		StringBuffer sb=new StringBuffer(seq);
		int seqlength=seq.length();
		//calculate number of mutations 
		int numberTransitMutations=(int)(seq.length()*transit);
		int numberTransverMutations=(int)(seq.length()*transver);
		for(int i=0;i<numberTransitMutations+numberTransverMutations;i++){
			//and distribute them
			int randPosition=(int)(Math.random()*seqlength);
			sb.setCharAt(randPosition, getMutation(seq.charAt(randPosition),GC,transit,transver));
			
		}
		
		return sb.toString();
	}
	
	public static String randomizeSequenceSelfMutate(String seq,double mutationProb,double GC){
		StringBuffer sb=new StringBuffer(seq);
		int seqlength=seq.length();
		
		int numberMutations=(int)(seq.length()*mutationProb);
		for(int i=0;i<numberMutations;i++){
			int randPosition=(int)(Math.random()*seqlength);
			sb.setCharAt(randPosition, getMutationSelfMutate(seq.charAt(randPosition),GC));
			
		}
		return sb.toString();
	}
	public static char getMutationSelfMutate(char in,double GC){

			char base=getRandomBase(GC);

				return base;

	}
	
	public static char getMutation(char in,double GC){
		boolean mut=false;
		while(mut!=true){
			char base=getRandomBase(GC);
			if(in!=base){
				return base;
			}
		}
		return 0;
	}
	
	public static char getMutation(char in,double transit,double transver,double GC){
			char base=getRandomBase(transit,transver,in);
			return base;
			
	}
	
	public static BitSet mutate(BitSet seq,int pos,int mutation){
		BitSet mut=new BitSet();
		mut.or(seq);
		if(pos%2==1){
			System.err.println("Position is wrong, has to be dividable by two!");
		}else{
			if(mut.get(pos)!= ((mutation)%2==0) ||mut.get(pos+1)!= ((mutation/2)%2==0)  ){
				mut.set(pos,(mutation)%2==0);
				mut.set(pos+1,(mutation/2)%2==0);
				return mut;
			}
		}
		return new BitSet();
	}
	
	public static BitSet reverse(BitSet code){
		BitSet rev=new BitSet();
		int j=0;
		for(int i=code.length()-3;i>-1;i=i-2){
			if(code.get(i)==false && code.get(i+1)==false){
				//T
				rev.set(j*2,false);
				rev.set(j*2+1,true);
			}else if(code.get(i)==false && code.get(i+1)==true){
				//A
				rev.set(j*2,false);
				rev.set(j*2+1,false);
			}else if(code.get(i)==true && code.get(i+1)==false){
				//G
				rev.set(j*2,true);
				rev.set(j*2+1,true);
			}else if(code.get(i)==true && code.get(i+1)==true){
				//C
				rev.set(j*2,true);
				rev.set(j*2+1,false);
			}
			j++;
			
		}
		rev.set(j*2,true);
		return rev;
	}
	public static BitSet flip(BitSet code){
		BitSet rev=new BitSet();
		int j=0;
		for(int i=code.length()-3;i>-1;i=i-2){
				rev.set(j*2,code.get(i));
				rev.set(j*2+1,code.get(i+1));

			j++;
			
		}
		rev.set(j*2,true);
		return rev;
	}
	public static BitSet append(BitSet b1,BitSet b2){
		int b1size=b1.length();
		int b2size=b2.length();
		if(b1size==0)return b2;

		for(int i=b1size-1;i<b1size+b2size;i++){
			
			b1.set(i,b2.get(i-b1size+1));
		}
		return b1;
	}
	
	public static int next(BitSet reference,BitSet search,int from){
		int pos=-1;
		if(from%2!=0){
			System.err.println("\"From\" has to be dividable by 2!");
			return -1;
		}
		for(int i=from;i<reference.length()-search.length()+1;i+=2){
			if(DNAmanipulations.get(reference,i,i+search.length()-1).equals(search)){
				return i;
			}
		}
		return pos;
	}
	
	public static BitSet get(BitSet b,int start,int end){
		if(end>=b.length())end=b.length()-1;
		if(start>=b.length())return new BitSet();
		BitSet nB=b.get(start,end);
		nB.set(end-start);
		return nB;
	}
/*	public static BitSet complement(BitSet code){
		BitSet rev=new BitSet();
		int j=0;
		for(int i=0;i<code.length()-2;i=i+2){
			if(code.get(i)==false && code.get(i+1)==false){
				//T
				rev.set(j*2,false);
				rev.set(j*2+1,true);
			}else if(code.get(i)==false && code.get(i+1)==true){
				//A
				rev.set(j*2,false);
				rev.set(j*2+1,false);
			}else if(code.get(i)==true && code.get(i+1)==false){
				//G
				rev.set(j*2,true);
				rev.set(j*2+1,true);
			}else if(code.get(i)==true && code.get(i+1)==true){
				//C
				rev.set(j*2,true);
				rev.set(j*2+1,false);
			}
			j++;
			
		}
		rev.set(j*2,true);
		return rev;
	}
	*/
	public  static String translate(String DNA,HashMap<String,String> code){
		StringBuilder AA=new StringBuilder("");
		for(int i=0;i<DNA.length()-2;i+=3){
			AA.append(code.get(DNA.substring(i,i+3).toUpperCase()));
		}
		return AA.toString();
	}
	
	public static HashMap<String,String> code(){
		HashMap<String,String> code=new HashMap<String, String>();
		code.put("TTT", "F");
		code.put("TTC", "F");
		code.put("TTG", "L");
		code.put("TTA", "L");
		code.put("CTT", "L");
		code.put("CTC", "L");
		code.put("CTA", "L");
		code.put("CTG", "L");
		code.put("ATT", "I");
		code.put("ATC", "I");
		code.put("ATA", "I");
		code.put("ATG", "M");
		code.put("GTT", "V");
		code.put("GTC", "V");
		code.put("GTA", "V");
		code.put("GTG", "V");
		code.put("TCT", "S");
		code.put("TCC", "S");
		code.put("TCA", "S");
		code.put("TCG", "S");
		code.put("CCT", "P");
		code.put("CCC", "P");
		code.put("CCA", "P");
		code.put("CCG", "P");
		code.put("ACT", "T");
		code.put("ACC", "T");
		code.put("ACA", "T");
		code.put("ACG", "T");
		code.put("GCT", "A");
		code.put("GCC", "A");
		code.put("GCA", "A");
		code.put("GCG", "A");
		code.put("TAT", "Y");
		code.put("TAC", "Y");
		code.put("TAA", "*");
		code.put("TAG", "+");
		code.put("CAT", "H");
		code.put("CAC", "H");
		code.put("CAA", "Q");
		code.put("CAG", "Q");
		code.put("AAT", "N");
		code.put("AAC", "N");
		code.put("AAA", "K");
		code.put("AAG", "K");
		code.put("GAT", "D");
		code.put("GAC", "D");
		code.put("GAA", "E");
		code.put("GAG", "E");
		code.put("TGT", "C");
		code.put("TGC", "C");
		code.put("TGA", "#");
		code.put("TGG", "W");
		code.put("CGT", "R");
		code.put("CGC", "R");
		code.put("CGA", "R");
		code.put("CGG", "R");
		code.put("AGT", "S");
		code.put("AGC", "S");
		code.put("AGA", "R");
		code.put("AGG", "R");
		code.put("GGT", "G");
		code.put("GGC", "G");
		code.put("GGA", "G");
		code.put("GGG", "G");
		return code;
	}
	public  static ArrayList<BitSet> FastaToBitSet(ArrayList<Fasta> seqs){
		ArrayList<BitSet> newSeqs=new ArrayList<BitSet>();
		for(int i=0;i<seqs.size();i++){
			String sequence=seqs.get(i).getSequence().toUpperCase();
			if(sequence.contains("N")){
				System.err.println("Left out a "+sequence.length()+"bp region since it contained Ns.");
				continue;
			}
			newSeqs.add(DNAmanipulations.codeDNA(sequence));
		}
		return newSeqs;
	}
	public  static ArrayList<BitSet> toBitSet(ArrayList<String> seqs){
		ArrayList<BitSet> newSeqs=new ArrayList<BitSet>();
		for(int i=0;i<seqs.size();i++){
			if(seqs.get(i).toUpperCase().contains("N")){
				System.err.println("Left out a "+seqs.get(i).length()+"bp region since it contained Ns.");
				continue;
			}
			newSeqs.add(DNAmanipulations.codeDNA(seqs.get(i).toUpperCase()));
		}
		return newSeqs;
	}
	public static String decodeDNA(BitSet code){
		
		StringBuilder DNA=new StringBuilder("");
		for(int i=0;i<code.length()-1;i+=2){
			if(code.get(i)==false && code.get(i+1)==false){
				DNA.append('A');
			}else if(code.get(i)==false && code.get(i+1)==true){
				DNA.append('T');

			}else if(code.get(i)==true && code.get(i+1)==false){
				DNA.append('C');

			}else if(code.get(i)==true && code.get(i+1)==true){
				DNA.append('G');

			}
		}
		return DNA.toString();
	}
	

	public static String decodeDNA(long code,int length){
		StringBuilder DNA=new StringBuilder("");
		for(int i=0;i<length;i+=1){
			if(code%4==0){
				DNA.append('A');
			}else if(code%4==1){
				DNA.append('T');

			}else if(code%4==2){
				DNA.append('C');

			}else if(code%4==3){
				DNA.append('G');

			}
			code/=4;
		}
		return DNA.toString();
	}
	
	//for quick coding
	
	public static long codeDNALong(String dna){
		long code=0;
		long A=0;
		long T=1;
		long C=2;
		long G=3;
		for (int i=0;i<dna.length();i++){
			long newLetter=0;
			if(dna.charAt(i)=='A'){
				newLetter=A << i*2l;
				
			}else if(dna.charAt(i)=='T'){
				newLetter=T << i*2l;

			}else if(dna.charAt(i)=='C'){
				newLetter=C << i*2l;

			}else if(dna.charAt(i)=='G'){
				newLetter=G << i*2l;

			}else{
				System.err.println("Wrong letter in DNA sequence: "+dna.charAt(i));
				return -1;
			}
			code+=newLetter;
		}
		return code;
	}
	//for long coding...slow
//	public static BitSet codeDNA(String dna){
//		BitSet code=new BitSet(dna.length()*2);
//		for (int i=0;i<dna.length();i++){
//			if(dna.charAt(i)=='A'){
//				code.set(i*2,false);
//				code.set(i*2+1,false);
//			}else if(dna.charAt(i)=='T'){
//				code.set(i*2,false);
//				code.set(i*2+1,true);
//			}else if(dna.charAt(i)=='C'){
//				code.set(i*2,true);
//				code.set(i*2+1,false);
//			}else if(dna.charAt(i)=='G'){
//				code.set(i*2,true);
//				code.set(i*2+1,true);
//			}else{
//				System.err.println("Wrong letter in DNA sequence: "+dna.charAt(i));
//				return null;
//			}
//		}
//		return code;
//	}
	public static BitSet codeDNA(String dna){
		BitSet code=new BitSet(dna.length()*2);
		dna=dna.toUpperCase();
		for (int i=0;i<dna.length();i++){
			if(dna.charAt(i)=='A'){
				code.set(i*2,false);
				code.set(i*2+1,false);
			}else if(dna.charAt(i)=='T'){
				code.set(i*2,false);
				code.set(i*2+1,true);
			}else if(dna.charAt(i)=='C'){
				code.set(i*2,true);
				code.set(i*2+1,false);
			}else if(dna.charAt(i)=='G'){
				code.set(i*2,true);
				code.set(i*2+1,true);
			}else{
				System.err.println("Wrong letter in DNA sequence: "+dna.charAt(i));
				return null;
			}
		}
		code.set((dna.length())*2,true);
		return code;
	}
}
