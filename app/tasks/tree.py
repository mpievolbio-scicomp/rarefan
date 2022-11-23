import subprocess
import os
import sys
import shlex
import shutil
from app.models import Job as DBJob
from rq import get_current_job
from app import app
app.app_context().push()

def tree_task(run_dir, treefile=None):
    """ Generate a phylogenetic tree from all DNA sequence files in given directory.

    :param run_dir: The directory containing the sequence files.
    :type  run_dir: str (path)

    :param treefile: Name of the treefile (default: 'tmptree.nwk')
    :type  treefile: str

    :return: The return code and log from the tree generation.
    :rtype: tuple

    :raises RuntimeError: Directory does not contain any sequence files in fasta format.

    """
    logger = app.logger

    inputs = [os.path.join(run_dir, f) for f in os.listdir(run_dir) if f.split(".")[-1] in ["fas", "fna", "fn", "fasta", "fastn"]]

    for f in inputs:
        logger.debug("Found sequence file %s.", f)

    # Check if treefile exists.
    if treefile is None or treefile == "None":
        treefile = 'tmptree.nwk'
    in_treefile = os.path.join(run_dir, treefile)
    outdir = os.path.join(run_dir, 'out')
    out_treefile = os.path.join(outdir, treefile)
    if os.path.isfile(in_treefile):
        shutil.copy(in_treefile, out_treefile)

        return 0, "Copied {} to {}".format(in_treefile, out_treefile)

    redis_job =  get_current_job()
    dbjob = DBJob.objects.get(run_id=redis_job.meta['run_id'])
    dbjob.set_status('tree')

    # else (no treefile exists.)
    andi_command = "andi -j {}".format(" ".join(inputs))

    logger.debug("tree generation command: %s", andi_command)

    proc = subprocess.Popen(shlex.split(andi_command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=False)

    # andi writes distance matrix to stdout and log to stderr.
    dist, log = proc.communicate()

    with open(os.path.join(outdir, 'tmptree.dist'), 'wb') as ofh:
        ofh.write(dist)

    logger.info(dist)
    logger.info(log)

    # clustDist command
    cd_command = "clustDist {}".format(os.path.join(outdir, 'tmptree.dist'))

    proc = subprocess.Popen(shlex.split(cd_command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=False)

    # clustDist writes distance matrix to stdout and log to stderr.
    nwk, cd_log = proc.communicate()

    logger.info(nwk)
    logger.info(cd_log)

    with open(os.path.join(outdir, 'tmptree.nwk'), 'wb') as ofh:
        ofh.write(nwk)

    # Append stdout and stderr to logfile.
    with open(os.path.join(outdir, 'tree.log'), 'ab') as fh:
        fh.write(log + cd_log)

    return {'returncode': proc.returncode,
            'log': log + cd_log
            }


def empty_task():
    return {'returncode': 2,
            'log': """A phylogeny is only computed for more than 3 submitted genomes.
             Plotting is disabled."""
            }
