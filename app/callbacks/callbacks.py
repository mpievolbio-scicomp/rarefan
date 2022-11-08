""":module callbacks: Hosting all callback routines for rarefan tasks."""

from app.models import Job as DBJob
from app.utilities.checkers import parse_results
from app.tasks.email import email_task
from app.tasks.zip import zip_task
from app import app

logger = app.logger
def convert_dict_to_update(data, roots=None, return_dict=None):
    """
      Returns a new dict that can be used in Document.update(**dict),
      this is used for updating MongoEngine documents with dictionary
      serialized from json.
      >>> data = {'email' : 'email@example.com'}
      >>> convert_dict_to_update(data)
      {'set__email': 'email@example.com'}
      :param data: dictionary with update parameters
      :param roots: roots of nested documents - used for recursion
      :param return_dict: used for recursion
      :return: new dict

      :credits: https://gist.github.com/Visgean/e536e466207bf439983a
    """
    if return_dict is None:
        return_dict = {}
    if roots is None:
        roots = []

    for key, value in data.items():
        # Convert to str.
        key = str(key)
        if isinstance(value, dict):
            roots.append(key)
            convert_dict_to_update(value, roots=roots, return_dict=return_dict)
            roots.remove(key)  # go one level down in the recursion
        else:
            if roots:
                set_key_name = 'set__{roots}__{key}'.format(
                    roots='__'.join(roots), key=key)
            else:
                set_key_name = 'set__{key}'.format(key=key)
            return_dict[set_key_name] = value

    return return_dict

def rarefan_on_success(job, result):
    """ Callback for the 'rarefan' task if completed successfully. """

    # Get db object.
    dbjob_id = job.meta['dbjob_id']
    dbjob = DBJob.objects.get(id=dbjob_id)

    # Update returncode.
    dbjob.update(set__stages__rarefan__results__returncode=result['returncode'])

    # Parse output files and report counts.
    parsed = parse_results(dbjob['setup']['outdir'], dbjob['setup']['reference_strain'])

    dbjob.update(**convert_dict_to_update(parsed['status'], roots=['stages', 'rarefan', 'results', 'data_sanity']))
    dbjob.update(**convert_dict_to_update(parsed['counts'], roots=['stages', 'rarefan', 'results', 'counts']))
    # dbjob.update(set__stages__rarefan__results__data_sanity__rayts=parsed['status']['rayts'])
    # dbjob.update(set__stages__rarefan__results__data_sanity__nmers=parsed['status']['nmers'])
    # dbjob.update(set__stages__rarefan__results__data_sanity__repins=parsed['status']['repins'])
    # dbjob.update(set__stages__rarefan__results__counts__rayts=parsed['counts']['rayts'])
    # dbjob.update(set__stages__rarefan__results__counts__nmers=parsed['counts']['nmers'])
    # dbjob.update(set__stages__rarefan__results__counts__repins__0=parsed['counts']['repins'].get(0, 0))
    # dbjob.update(set__stages__rarefan__results__counts__repins__1=parsed['counts']['repins'].get(1, 0))
    # dbjob.update(set__stages__rarefan__results__counts__repins__2=parsed['counts']['repins'].get(2, 0))
    # dbjob.update(set__stages__rarefan__results__counts__repins__3=parsed['counts']['repins'].get(3, 0))
    # dbjob.update(set__stages__rarefan__results__counts__repins__4=parsed['counts']['repins'].get(4, 0))
    # dbjob.update(set__stages__rarefan__results__counts__repins__5=parsed['counts']['repins'].get(5, 0))
    # dbjob.update(set__stages__rarefan__results__counts__repins__6=parsed['counts']['repins'].get(6, 0))
    # dbjob.update(set__stages__rarefan__results__counts__repins__7=parsed['counts']['repins'].get(7, 0))
    # dbjob.update(set__stages__rarefan__results__counts__repins__8=parsed['counts']['repins'].get(8, 0))

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
        dbjob.update(set__stages__tree__results__log=result['log'])
    elif stage == "rayt_alignment":
        dbjob.update(set__stages__rayt_alignment__results__returncode=result['returncode'])
        dbjob.update(set__stages__rayt_alignment__results__log=result['log'])
    elif stage == "rayt_phylogeny":
        dbjob.update(set__stages__rayt_phylogeny__results__returncode=result['returncode'])
        dbjob.update(set__stages__rayt_phylogeny__results__log=result['log'])
    elif stage == "zip":
        dbjob.update(set__stages__zip__results__returncode=result['returncode'])
        dbjob.update(set__stages__zip__results__log=result['log'])

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
    elif stage == "rayt_alignment":
        dbjob.update(set__stages__rayt_alignment__results__returncode=1)
    elif stage == "rayt_phylogeny":
        dbjob.update(set__stages__rayt_phylogeny__results__returncode=1)
    elif stage == "zip":
        dbjob.update(set__stages__zip__results__returncode=1)

    dbjob.set_overall()
    dbjob.save()

    email_task(dbjob.run_id)
