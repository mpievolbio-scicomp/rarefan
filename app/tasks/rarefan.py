"""module rarefan: Implementation of the function that runs the rarefan java code as a system process. Can be called as a script with arguments or from a python session including
   as a redis task."""

import os, sys, shutil, shlex
import subprocess
import logging

from app.models import Job as DBJob
from rq import get_current_job

from app.utilities.rarefan_cli import rarefan_command

def rarefan_task(**kwargs):
    """Run the rarefan java code with arguments."""

    for k,v in kwargs.items():
        logging.debug("%s = %s", k, str(v))

    java_command = rarefan_command(**kwargs) 

    logging.info("Java command: %s", java_command)

    redis_job = get_current_job()
    logging.debug("In rarefan.py, redis job id = %s", redis_job.id)
    logging.debug("In rarefan.py, redis job meta = %s", str(redis_job.meta))
    dbjob = DBJob.objects.get(run_id=redis_job.meta['run_id'])
    dbjob.set_status('rarefan')

    proc = subprocess.Popen(shlex.split(java_command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, shell=False)

    log, _ = proc.communicate()

    # Append stdout and stderr to logfile.
    with open(os.path.join(kwargs['tmpdir'], 'out', 'rarefan.log'), 'ab') as fh:
        fh.write(log)

    return {'returncode': proc.returncode,
            'log': log
            }

