import subprocess
import os
import sys
import shlex
import shutil
from app.models import Job as DBJob
from rq import get_current_job

import logging
logger = logging.getLogger('rarefan')

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

    # else (no treefile exists.)
    command = "andi -j {} | clustDist > {}".format(" ".join(inputs), out_treefile)

    redis_job =  get_current_job()
    dbjob = DBJob.objects.get(run_id=redis_job.meta['run_id'])
    dbjob.set_status('tree')

    logger.debug("tree generation command: %s", command)

    proc = subprocess.Popen(command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=True)

    log, _ = proc.communicate()

    # Append stdout and stderr to logfile.
    with open(os.path.join(outdir, 'rarefan.log'), 'ab') as fh:
        fh.write(log)

    return {'returncode': proc.returncode,
            'log': log
            }


def empty_task():
    return {'returncode': 0,
            'log': "No treefile as per job parameters."
            }
