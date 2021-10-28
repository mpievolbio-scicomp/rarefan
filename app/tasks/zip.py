import subprocess
import os
import sys
import shlex
import shutil


def zip_task(run_dir):
    """ Compress all output data from a given run

    :param run_dir: The directory containing the output directory 'out/'
    :type  run_dir: str (path)

    :return: The return code and log from the zip process.
    :rtype: tuple

    """
    oldwd = os.getcwd()
    os.chdir(run_dir)
    run_id = os.path.split(run_dir)[-1]
    zip_command = " ".join(["zip",
                                "-r",
                                 run_id + "_out.zip",
                                'out'
                                ]
                               )


    proc = subprocess.Popen(shlex.split(zip_command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=False)

    log, _ = proc.communicate()

    # Append stdout and stderr to logfile.
    with open(os.path.join('out', 'rarefan.log'), 'ab') as fh:
        fh.write(log)

    os.chdir(oldwd)

    return proc.returncode, log

