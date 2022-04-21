suppressMessages(library(Biostrings))
suppressMessages(library(ape))
suppressMessages(library(muscle))
suppressMessages(library(stringr))
suppressMessages(library(logging))
suppressMessages(library(argparse))

parser <- ArgumentParser(description='Run muscle alignment on provided fasta data.')

# Input data dir
parser$add_argument('-i', '--input',
                    metavar='INPUT',
                    type="character",
                    dest='raytseqFile',
                    help='The fasta data to align.'
)


args = parser$parse_args(commandArgs(TRUE))

raytseqFile = args$raytseqFile
raytAlnFile = 'raytAln.phy'
# The associated rayt sequence file.

# Read sequences
raytseqs=readDNAStringSet(raytseqFile,format="fasta")

# If no rayt sequences found, return empty.
if(length(raytseqs) == 0) {
  logging::logwarn(paste0("No RAYT sequences found in ", raytseqFile))
}


# Run muscle alignment.
logging::logdebug("Running muscle...")
aln=muscle(raytseqs, verbose=T, log='/tmp/muscle.log')
logging::logdebug("...done.")

# Write alignment to file.
write.phylip(aln,raytAlnFile)

