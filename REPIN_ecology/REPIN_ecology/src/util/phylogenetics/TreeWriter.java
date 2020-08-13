package util.phylogenetics;
import java.io.*;






	/**
	 * Static utility class to write a tree to a file.  Does not create an object, no constructor.
	 * Newick files are created, but the functions here could be used to create a Nexus
	 * file writer.
	 * @author jslack
	 *
	 */
	public class TreeWriter  {

		/** End of tree character. */
	    private final static char lineTerminator = ';';
	    /** End of tree character. */
	    private final static char treeTerminator = lineTerminator;
	    /** Open bracket to wrap children of a subtree root. */
	    private final static char openBracket = '(';
	    /** Close bracket to wrap children of a subtree root. */
	    private final static char closeBracket = ')';
	    /** Place between children of a subtree root. */
	    private final static char childSeparator = ',';
	    /** Quote object to wrap node names. */
	    private final static char doubleQuote = '"';


		/**
		 * Writes Tree to file by name of fileName.
		 * @param t Tree to be written to file.
		 * @param fileName output file name.
		 */
	    static public void writeTree(Tree t, String fileName,boolean append){
	    	writeTree(t, fileName,true,append);
	    }
		static public void writeTree(Tree t, String fileName,boolean dq,boolean append)
		{
			try
			{
				FileWriter treeFileWriter = new FileWriter(fileName,append);
				BufferedWriter writer = new BufferedWriter(treeFileWriter);
				nodeWriter(t.getRoot(), writer,dq);
				writer.write(treeTerminator);
				writer.close();
			}
			catch (IOException ioe)
			{
				System.err.println("Could not find file: " + fileName);
			}
		}
		
		/**
		 * Write a single node to the writer object, after calling the subtree writer on each child.
		 * All nodes are quoted with {@link #doubleQuote}.
		 * @param currNode Node to write.
		 * @param writer Writer object used to send the node and its subtree to file.
		 * @throws IOException
		 */
		static private void nodeWriter(TreeNode currNode, BufferedWriter writer) throws IOException
		{
			nodeWriter(currNode, writer,true);
		}
		static private void nodeWriter(TreeNode currNode, BufferedWriter writer,boolean dq) throws IOException
		{
			String quote=dq?"\"":"";
			if (currNode.numberChildren() > 0)
			{
				subtreeWriter(currNode, writer,dq);
			}
			
			if (currNode.getName().length() > 0)
				writer.write(quote + currNode.getName() + quote+":"+currNode.getWeight());	
			else 
				writer.write(":"+currNode.getWeight());	
		}
		
		/**
		 * Wrap a series of nodes in brackets, then call the node writer on each child of the input subtree root.
		 * @param subtreeRoot Root of the subtree to write.  This node itself is not written by this function.
		 * @param writer Writer object used to send the subtree to file.
		 * @throws IOException
		 */
		static private void subtreeWriter(TreeNode subtreeRoot, BufferedWriter writer,boolean dq)
		throws IOException
		{
			writer.write(openBracket);
			for (int i = 0; i < subtreeRoot.numberChildren() - 1; i++)
			{
				nodeWriter(subtreeRoot.getChild(i), writer,dq);
				writer.write(childSeparator);
			}
			nodeWriter(subtreeRoot.getChild(subtreeRoot.numberChildren()-1), writer,dq);
			writer.write(closeBracket);
		}

		
	
		
		


}
