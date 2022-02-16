""" :module rayt_phylo: Hosting the redis task to compute rayt phylogenies."""

import os, sys, shutil, shlex
import subprocess

from app.models import Job as DBJob
from rq import get_current_job
from app import app

app.app_context().push()
logger = app.logger

from app.utilities.rarefan_cli import rarefan_command

def phylogeny_task(**kwargs):
    """Run the rarefan java code with arguments."""
    pass

def alignment_task(**kwargs):
    """ Run the RAYT alignment """

    for k,v in kwargs.items():
        logger.debug("%s = %s", k, str(v))

    redis_job = get_current_job()
    dbjob = DBJob.objects.get(run_id=redis_job.meta['run_id'])
    dbjob.set_status('rarefan')

    log, ret = run_alignment(run_dir)

    # Append stdout and stderr to logfile.
    with open(os.path.join(kwargs['tmpdir'], 'out', 'rarefan.log'), 'ab') as fh:
        fh.write(log)

    return {'returncode': ret,
            'log': log
            }

def run_alignment(run_dir):
    """ Workhorse function to run the muscle aligner. """
    input_fname = os.path.join(run_dir, 'out', 'repin_rayt_association.txt.fas')
    output_fname = os.path.join(run_dir, 'out', 'raytAln.phy')
    command = 'muscle -in {} -out {} -phyiout'.format(input_fname, output_fname)

    proc = subprocess.Popen(shlex.split(command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, shell=False)

    log, _ = proc.communicate()

    return proc.returncode, log

def run_phyml(run_dir, seed=None):
    """ Workhorse function to run the phyml tool.

    :param run_dir: The rarefan run directory containing the user submitted input data and the out/ directory created by rarefan.
    :param seed: Random seed for phyml. Mainly used for testing.

    """

    input_fname = os.path.join(run_dir, 'out', 'raytAln.phy')

    command = 'phyml --quiet -i {} -m GTR'.format(input_fname)

    if seed is not None:
        command = command + " --r_seed {}".format(int(seed))

    proc = subprocess.Popen(shlex.split(command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, shell=False,
                            cwd=os.path.join(run_dir, 'out'))

    log, _ = proc.communicate()

    return proc.returncode, log



