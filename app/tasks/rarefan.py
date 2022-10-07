"""
module rarefan: Implementation of the function that runs the rarefan java code as a system process.
"""

import os, sys, shutil, shlex
import subprocess

from app.models import Job as DBJob
from rq import get_current_job
from app import app

app.app_context().push()
logger = app.logger

from app.utilities.rarefan_cli import rarefan_command

def rarefan_task(**kwargs):
    """Run the rarefan java code with arguments."""

    for k,v in kwargs.items():
        logger.debug("%s = %s", k, str(v))

    java_command = rarefan_command(**kwargs) 

    logger.info("Java command: %s", java_command)

    redis_job = get_current_job()
    run_id=redis_job.meta['run_id']
    logger.debug("In rarefan.py, redis job id = %s", redis_job.id)
    logger.debug("In rarefan.py, redis job meta = %s", str(redis_job.meta))
    dbjob = DBJob.objects.get(run_id=run_id)
    dbjob.set_status('rarefan')

    proc = subprocess.Popen(shlex.split(java_command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, shell=False, )

    log, _ = proc.communicate()
    logger.debug(log)

    log = log.replace(b'Wrong letter in DNA sequence: |', b'')

    # Append stdout and stderr to logfile.
    with open(os.path.join(kwargs['tmpdir'], 'out', 'rarefan.log'), 'ab') as fh:
        run_id_str = "# RAREFAN run {}\n".format(run_id)
        fh.write(run_id_str.encode('ascii'))
        fh.write(log)

    return {'returncode': proc.returncode,
            'log': log
            }

