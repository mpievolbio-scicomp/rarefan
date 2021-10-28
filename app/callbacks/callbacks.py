import os
import rq
import copy

from app import db
from app.models import Job as DBJob
from app.utilities.checkers import parse_results

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

    rc = copy.deepcopy(result['returncode'])
    dbjob.update(set__stages__rarefan__results__returncode=rc)
    # dbjob.update(set__stages__rarefan__results__log=result['log'])

    results = parse_results(dbjob['setup']['tmpdir'], dbjob['setup']['reference_strain'])

    dbjob.update(set__stages__rarefan__results__data_sanity__rayts=results['status']['rayts'])
    dbjob.update(set__stages__rarefan__results__data_sanity__seeds=results['status']['seeds'])
    dbjob.update(set__stages__rarefan__results__data_sanity__repins=results['status']['repins'])
    dbjob.update(set__stages__rarefan__results__counts__rayts=results['counts']['rayts'])
    dbjob.update(set__stages__rarefan__results__counts__seeds=results['counts']['seeds'])
    dbjob.update(set__stages__rarefan__results__counts__repins=results['counts']['repins'])



def on_failure(job, connection, type, value, traceback):

    logging.debug(job.id)
    logging.debug(job.meta)
    logging.debug(type)
    logging.debug(value)
    logging.debug(traceback)
