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
    pass

