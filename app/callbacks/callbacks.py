""":module callbacks: Hosting all callback routines for rarefan tasks."""

from app.models import Job as DBJob
from app.utilities.checkers import parse_results
from app.tasks.email import email_task
from app.tasks.zip import zip_task
from app import app

logger = app.logger

def rarefan_on_success(job, result):
    """ Callback for the 'rarefan' task if completed successfully. """

    # Get db object.
    dbjob_id = job.meta['dbjob_id']
    dbjob = DBJob.objects.get(id=dbjob_id)

    # Update returncode.
    dbjob.update(set__stages__rarefan__results__returncode=result['returncode'])

    # Parse output files and report counts.
    parsed = parse_results(dbjob['setup']['outdir'], dbjob['setup']['reference_strain'])

    dbjob.update(set__stages__rarefan__results__data_sanity__rayts=parsed['status']['rayts'])
    dbjob.update(set__stages__rarefan__results__data_sanity__nmers=parsed['status']['nmers'])
    dbjob.update(set__stages__rarefan__results__data_sanity__repins=parsed['status']['repins'])
    dbjob.update(set__stages__rarefan__results__counts__rayts=parsed['counts']['rayts'])
    dbjob.update(set__stages__rarefan__results__counts__nmers=parsed['counts']['nmers'])
    dbjob.update(set__stages__rarefan__results__counts__repins__0=parsed['counts']['repins'].get(0, 0))
    dbjob.update(set__stages__rarefan__results__counts__repins__1=parsed['counts']['repins'].get(1, 0))
    dbjob.update(set__stages__rarefan__results__counts__repins__2=parsed['counts']['repins'].get(2, 0))
    dbjob.update(set__stages__rarefan__results__counts__repins__3=parsed['counts']['repins'].get(3, 0))
    dbjob.update(set__stages__rarefan__results__counts__repins__4=parsed['counts']['repins'].get(4, 0))
    dbjob.update(set__stages__rarefan__results__counts__repins__5=parsed['counts']['repins'].get(5, 0))
    dbjob.update(set__stages__rarefan__results__counts__repins__6=parsed['counts']['repins'].get(6, 0))
    dbjob.update(set__stages__rarefan__results__counts__repins__7=parsed['counts']['repins'].get(7, 0))
    dbjob.update(set__stages__rarefan__results__counts__repins__8=parsed['counts']['repins'].get(8, 0))
    dbjob.update(set__stages__rarefan__results__counts__repins__total=sum(parsed['counts']['repins'].values()))

    dbjob.save()


def on_success(job, connection, result, *args, **kwargs):
    """ Generic callback on success. """

    # Get db object.
    dbjob_id = job.meta['dbjob_id']
    stage = job.meta['stage']
    dbjob = DBJob.objects.get(id=dbjob_id)

    if stage == "rarefan":
        rarefan_on_success(job, result)
    elif stage == "tree":
        dbjob.update(set__stages__tree__results__returncode=result['returncode'])
    elif stage == "zip":
        dbjob.update(set__stages__zip__results__returncode=result['returncode'])

    dbjob.set_overall()

    dbjob.save()


def on_failure(job, connection, type, value, traceback):
    """ Generic failure callback for all tasks. """

    logger.debug(job.id)
    logger.debug(job.meta)
    logger.debug(type)
    logger.debug(value)
    logger.debug(traceback)

    dbjob_id = job.meta['dbjob_id']
    stage = job.meta['stage']
    dbjob = DBJob.objects.get(id=dbjob_id)

    if stage == "rarefan":
        dbjob.update(set__stages__rarefan__results__returncode=1)
    elif stage == "tree":
        dbjob.update(set__stages__tree__results__returncode=1)

        # Even though getting the tree was not successfull, we still want to generate the zip file.
        zip_results = zip_task(dbjob.setup['tmpdir'])
        dbjob.update(set__stages__zip__results__returncode=zip_results['returncode'])
    elif stage == "zip":
        dbjob.update(set__stages__zip__results__returncode=1)

    dbjob.set_overall()
    dbjob.save()

    email_task(dbjob.run_id)
