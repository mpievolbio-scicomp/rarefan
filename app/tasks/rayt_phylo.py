""" :module rayt_phylo: Hosting the redis task to compute rayt phylogenies."""

import os, sys, shutil, shlex
import subprocess

from app.models import Job as DBJob
from rq import get_current_job
import logging
from app import app

app.app_context().push()
logger = app.logger
logger.setLevel(logging.DEBUG)

from app.utilities.rarefan_cli import rarefan_command
utilities_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..', 'utilities'))


def phylogeny_task(**kwargs):
    """ Run the RAYT phylogeny calculator. """

    for k,v in kwargs.items():
        logger.debug("%s = %s", k, str(v))

    redis_job = get_current_job()
    dbjob = DBJob.objects.get(run_id=redis_job.meta['run_id'])
    dbjob.set_status('rayt_phylogeny')
    run_dir = dbjob.setup['tmpdir']

    log, ret = run_phyml(run_dir)

    # Append stdout and stderr to logfile.
    with open(os.path.join(run_dir, 'out', 'rarefan.log'), 'ab') as fh:
        fh.write(log)

    return {'returncode': ret,
            'log': log
            }


def alignment_task(**kwargs):
    """ Run the RAYT alignment """

    for k,v in kwargs.items():
        logger.debug("%s = %s", k, str(v))

    redis_job = get_current_job()
    dbjob = DBJob.objects.get(run_id=redis_job.meta['run_id'])
    dbjob.set_status('rayt_alignment')
    run_dir = dbjob.setup['tmpdir']

    log, ret = run_alignment(run_dir)

    # Append stdout and stderr to logfile.
    with open(os.path.join(run_dir, 'out', 'rarefan.log'), 'ab') as fh:
        fh.write(log)

    return {'returncode': ret,
            'log': log
            }

def run_alignment(run_dir):
    """ Workhorse function to run the muscle aligner. """
    input_fname = os.path.join(run_dir, 'out', 'repin_rayt_association.txt.fas')
    output_fname = os.path.join(run_dir, 'out', 'raytAln.phy')
    R = shutil.which('Rscript')
    rscript = os.path.join(utilities_path, 'run_muscle.R')
    logger.debug(rscript)
    command = '{} {} -i {}'.format(R, rscript,  input_fname)

    logger.debug(command)

    proc = subprocess.Popen(shlex.split(command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, shell=False,
                            cwd=os.path.join(run_dir, 'out'))

    log, _ = proc.communicate()
    logger.debug(log)

    return log, proc.returncode

def run_phyml(run_dir, seed=None):
    """ Workhorse function to run the phyml tool.

    :param run_dir: The rarefan run directory containing the user submitted input data and the out/ directory created by rarefan.
    :param seed: Random seed for phyml. Mainly used for testing.

    """

    input_fname = os.path.join(run_dir, 'out', 'raytAln.phy')

    command = 'phyml --quiet -i {} -m GTR'.format(input_fname)


    if seed is not None:
        command = command + " --r_seed {}".format(int(seed))

    logger.debug(command)

    proc = subprocess.Popen(shlex.split(command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, shell=False,
                            cwd=os.path.join(run_dir, 'out'))

    log, _ = proc.communicate()
    logger.debug(log)

    return log, proc.returncode



