import os
import rq
import copy

from app import db
from app.models import Job as DBJob
from app.utilities.checkers import parse_results
from app.tasks.email import email_task

import logging
logging.basicConfig(level=logging.DEBUG)

def rarefan_on_success(job, connection, result, *args, **kwargs):

    # Get db object.
    run_id = job.meta['run_id']
    dbjob_id = job.meta['dbjob_id']
    redis_job_id = job.id
    dbjob = DBJob.objects.get(id=dbjob_id)

    logging.debug(run_id)
    logging.debug(dbjob_id)
    logging.debug(redis_job_id)
    logging.debug(list(result.keys()))

    # Update returncode.
    dbjob.update(set__stages__rarefan__results__returncode=result['returncode'])

    # Parse output files and report counts.
    parsed = parse_results(dbjob['setup']['outdir'], dbjob['setup']['reference_strain'])

    dbjob.update(set__stages__rarefan__results__data_sanity__rayts=parsed['status']['rayts'])
    dbjob.update(set__stages__rarefan__results__data_sanity__nmers=parsed['status']['nmers'])
    dbjob.update(set__stages__rarefan__results__data_sanity__repins=parsed['status']['repins'])
    dbjob.update(set__stages__rarefan__results__counts__rayts=parsed['counts']['rayts'])
    dbjob.update(set__stages__rarefan__results__counts__nmers=parsed['counts']['nmers'])
    dbjob.update(set__stages__rarefan__results__counts__repins=sum(parsed['counts']['repins'].values()))

    dbjob.save()



def on_failure(job, connection, type, value, traceback):

    logging.debug(job.id)
    logging.debug(job.meta)
    logging.debug(type)
    logging.debug(value)
    logging.debug(traceback)

    dbjob_id = job.meta['dbjob_id']
    dbjob = DBJob.objects.get(id=dbjob_id)

    email_task(dbjob)
    raise RuntimeError
