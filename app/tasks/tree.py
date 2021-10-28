import subprocess
import os
import sys
import shlex
import shutil

def tree_task(run_dir, treefile='tmptree.nwk'):
    """ Generate a phylogenetic tree from all DNA sequence files in given directory.

    :param run_dir: The directory containing the sequence files.
    :type  run_dir: str (path)

    :return: The return code and log from the tree generation.
    :rtype: tuple

    :raises RuntimeError: Directory does not contain any sequence files in fasta format.

    """

    inputs = [os.path.join(run_dir, f) for f in os.listdir(run_dir) if f.split(".")[-1] in ["fas", "fna", "fn", "fasta", "fastn"]]

    # Check if treefile exists.
    tf = os.path.join(run_dir, treefile)
    outdir = os.path.join(run_dir, 'out')
    if os.path.isfile(tf):
        shutil.copy(tf, os.path.join(outdir, 'tmptree.nwk'))

        return 0, "Copied {} to {}".format(tf, os.path.join(outdir, 'tmptree.nwk')), ""

    # else (no treefile exists.)
    command = "andi -j {} | clustDist > {}".format(" ".join(inputs), treefile)

    proc = subprocess.Popen(shlex.split(command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=False)

    log, _ = proc.communicate()

    # Append stdout and stderr to logfile.
    with open(os.path.join(outdir, 'rarefan.log'), 'ab') as fh:
        fh.write(log)

    return proc.returncode, log

