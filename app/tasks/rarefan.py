"""module rarefan: Implementation of the function that runs the rarefan java code as a system process. Can be called as a script with arguments or from a python session including
   as a redis task."""

import os, sys, shutil, shlex
import subprocess
import logging

JAR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar'))
mcl_threads = max(os.cpu_count()//4, 1)


def rarefan_task(**kwargs):
    """Run the rarefan java code with arguments."""

    for k,v in kwargs.items():
        logging.debug("%s = %s", k, str(v))

    java_command = " ".join(['java',
                            '-Dcom.sun.management.jmxremote',
                            '-Dcom.sun.management.jmxremote.port=9010',
                            '-Dcom.sun.management.jmxremote.local.only=true',
                            '-Dcom.sun.management.jmxremote.authenticate=false',
                            '-Dcom.sun.management.jmxremote.ssl=false',
                            '-Xmx10g',
                            '-jar',
                            JAR,
                            kwargs['tmpdir'],
                            kwargs['outdir'],
                            kwargs['reference_strain'],
                            '{0:s}'.format(kwargs['min_nmer_occurence']),
                            '{0:s}'.format(kwargs['nmer_length']),
                            kwargs['query_rayt_fname'],
                            kwargs['treefile'],
                            '{0:s}'.format(kwargs['e_value_cutoff']),
                            {"y": "true", None: "false"}[kwargs['analyse_repins']],
                            '{0:d}'.format(mcl_threads),
                            ]
                           )

    logging.info("Java command: %s", java_command)

    proc = subprocess.Popen(shlex.split(java_command),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, shell=False)

    log, _ = proc.communicate()

    # Append stdout and stderr to logfile.
    with open(os.path.join(kwargs['tmpdir'], 'rarefan.log'), 'ab') as fh:
        fh.write(log)

    return {'returncode': proc.returncode,
            'log': log
            }


if __name__ == "main":
    rarefan_task(*sys.argv)
