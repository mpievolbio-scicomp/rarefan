import subprocess
import os
import sys
import shlex
import shutil

from app.models import Job as DBJob
from rq import get_current_job
from app import app
app.app_context().push()
logger = app.logger


def zip_task(run_dir):
    """ Compress all output data from a given run

    :param run_dir: The directory containing the output directory 'out/'
    :type  run_dir: str (path)

    :return: The return code and log from the zip process.
    :rtype: tuple

    """
    oldwd = os.getcwd()
    logger.debug("Saving old WD %s", oldwd)
    os.chdir(run_dir)
    logger.debug("Chdir to %s", os.getcwd())
    run_id = os.path.split(run_dir)[-1]
    zip_command = " ".join(["zip",
                                "-r",
                                "--quiet",
                                 run_id + "_out.zip",
                                'out'
                                ]
                               )

    logger.debug("zip command: %s", zip_command)
    redis_job =  get_current_job()
    dbjob = DBJob.objects.get(run_id=redis_job.meta['run_id'])
    dbjob.set_status('zip')

    proc = subprocess.Popen(shlex.split(zip_command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=False)

    log, err = proc.communicate()

    logger.debug("proc.stdout: %s", log)
    logger.debug("proc.stderr: %s", err)

    # Append only stderr to logfile.
    with open(os.path.join('out', 'zip.log'), 'ab') as fh:
        fh.write(err)

    os.chdir(oldwd)
    logger.debug("Chdir to %s", os.getcwd())

    return {'returncode': proc.returncode, "log": log}

