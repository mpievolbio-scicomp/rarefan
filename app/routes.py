from flask import render_template
from flask import request
from flask import session
from flask import redirect
from flask import url_for
from flask import abort
from flask import send_from_directory
from flask import flash

from werkzeug.utils import secure_filename

from redis import Redis
import rq
from rq.job import Job as RQJob
from rq.exceptions import NoSuchJobError

from app.views import SubmitForm, AnalysisForm, UploadForm, ReturnToResultsForm, RunForm
from app import app, db
from app.models import Job
from app.tasks.rarefan import rarefan_task
from app.tasks.tree import tree_task
from app.tasks.zip import zip_task
from app.tasks.email import email_task, email_test

from app.callbacks.callbacks import rarefan_on_success, on_failure

import copy
import os
import shlex
import shutil
import stat
import subprocess
import tempfile
import logging
from Bio import SeqIO

import datetime
import logging
logging.basicConfig(level=logging.DEBUG)

def get_logger():
    logger = logging.getLogger(__name__)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(module)s: %(message)s')
    timestamp = datetime.datetime.now().strftime(format="%Y%m%d-%H%M%S")
    handler = logging.FileHandler("/tmp/rarefan_{}.log".format(timestamp))
    handler.setFormatter(formatter)
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)

    return logger

logger = get_logger()


def get_status_code(run_id_path):
    # Check if the run has finished.
    is_started = ".start.stamp" in os.listdir(run_id_path)
    is_java_finished = ".java.stamp" in os.listdir(run_id_path)
    is_zip_finished = ".zip.stamp" in os.listdir(run_id_path)

    status = is_started * 1 + is_java_finished * 10 + is_zip_finished * 100

    return status


def validate_fasta(filename):
    """ Validates input from passed file as fasta formatted sequence data.

    :param filename: The filename of the file to validate.

    """
    logger.info("Validating fasta file %s.", filename)
    with open(filename, 'r') as fp:
        fasta = SeqIO.parse(fp, "fasta")
        is_fasta = any(fasta)

        if not is_fasta:
            logger.warning("%s is not a valid fasta file.", filename)
        return is_fasta


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/upload', methods=['GET', 'POST'])
def upload():
    logging.debug("upload/%s", request.method)
    if request.method == 'POST': #upload_form.validate_on_submit():
        session['tmpdir'] = tempfile.mkdtemp(
            suffix=None,
            prefix="",
            dir=app.config["UPLOAD_DIR"]
        )

        seqs = [v for k,v in request.files.items() if k.startswith('file')]
        logger.info("Uploading %s.", str(seqs))

        dna_extensions = ['fn', 'fna', 'fastn', 'fas', 'fasta']
        aa_extensions = ['fa', 'faa']
        tree_extensions = ['nwk']
        fnames = [os.path.join(session['tmpdir'], secure_filename(seq.filename)) for seq in seqs]
        [seq.save(fname) for seq, fname in zip(seqs, fnames)]

        # Check if valid fasta files.
        for f in fnames:
            if f.split('.')[-1] in tree_extensions:
                continue
            is_fasta = validate_fasta(f)
            if not is_fasta:
                flash("{} is neither a valid fasta file nor a tree file and will be ignored.".format(os.path.basename(f)))
                fnames.remove(f)
                os.remove(f)

        basenames = [os.path.basename(f) for f in fnames]
        strain_names = [bn for bn in basenames if bn.split(".")[-1] in dna_extensions]
        rayt_names = [bn for bn in basenames if bn.split(".")[-1] in aa_extensions]
        tree_names = [bn for bn in basenames if bn.split(".")[-1] in tree_extensions]

        session['strain_names'] = strain_names
        session['rayt_names'] = rayt_names
        session['tree_names'] = tree_names

        for k,v in session.items():
            logging.debug("session[%s] = %s", k, str(v))

        return redirect(url_for('submit', _method='GET'))

        logging.error("How on earth did you get here?????")

    else:
        form = RunForm()
        session['tmpdir'] = None
        session['strain_names'] = None
        session['rayt_names'] = None
        session['tree_names'] = None
        session['outdir'] = None
        session['reference_strain'] = None
        session['query_rayt'] = None
        session['treefile'] = None
        session['min_nmer_occurence'] = None
        session['nmer_length'] = None
        session['e_value_cutoff'] = None
        session['analyse_repins'] = None
        session['email'] = None

        return render_template(
            'upload.html',
            title="Upload sequences",
            confirmation_form=form
        )


@app.route('/submit', methods=['GET', 'POST'])
def submit():

    logging.debug("submit/%s", request.method)
    submit_form = SubmitForm()

    strain_names = session.get('strain_names')
    submit_form.reference_strain.choices.extend(strain_names)
    submit_form.query_rayt.choices.extend(session.get('rayt_names'))
    submit_form.treefile.choices.extend(["None"] + session.get('tree_names'))

    if submit_form.validate_on_submit():
        tmpdir = session['tmpdir']
        session['outdir'] = os.path.join(tmpdir, 'out')

        if os.path.isdir(session['outdir']):
            logging.warning("Rerun")
        else:
            os.mkdir(session['outdir'])

        session['reference_strain'] = request.form.get('reference_strain')
        session['query_rayt'] = request.form.get('query_rayt')
        session['min_nmer_occurence'] = request.form.get('min_nmer_occurence')
        treefile = request.form.get('treefile', None)
        run_id = os.path.basename(session['tmpdir'])

        session['treefile'] = treefile
        session['nmer_length'] = request.form.get('nmer_length')
        session['e_value_cutoff'] = request.form.get('e_value_cutoff')
        session['analyse_repins'] = request.form.get('analyse_repins')
        session['email'] = request.form.get('email', None)

        logger.info("Session parameters:")
        logger.info("tmpdir: %s", session['tmpdir'])
        logger.info("outdir: %s", session['outdir'])
        logger.info("reference_strain: %s", session['reference_strain'])
        logger.info("treefile: %s", session['treefile'])
        logger.info("email: %s", session['email'])

        # Store session in db.
        run_id = os.path.basename(session['tmpdir'])
        dbjob = Job(run_id=run_id)
        dbjob.setup = copy.deepcopy(session)

        dbjob.save()

        # If one of the server provided rayt files was selected,  copy it to the working dir. In the dropdown menu,
        # the server provided rayts are listed without filename extension, so have to append that here.
        query_rayt_fname = os.path.join(session['tmpdir'], session['query_rayt'])
        if session['query_rayt'] in ['yafM_Ecoli', 'yafM_SBW25']:
            query_rayt_fname = query_rayt_fname + ".faa"
            src = os.path.join(app.static_folder, "rayts", session['query_rayt'] + ".faa")
            logging.debug("Copying rayt from %s to %s.", src, query_rayt_fname)
            shutil.copyfile(src, query_rayt_fname)

            if not os.path.isfile(query_rayt_fname):
                raise IOError("Copying %s to %s failed." % (src, query_rayt_fname))

        rarefan_job = RQJob.create(
            rarefan_task,
            connection=app.redis,
            on_success=rarefan_on_success,
            on_failure=on_failure,
            meta={'run_id': 'run_id', 'dbjob_id': dbjob.id},
            kwargs={
                "tmpdir": session['tmpdir'],
                "outdir": session['outdir'],
                "reference_strain": session['reference_strain'],
                "min_nmer_occurence": session['min_nmer_occurence'],
                "nmer_length": session['nmer_length'],
                "query_rayt_fname": query_rayt_fname,
                "treefile": session['treefile'],
                "e_value_cutoff": session['e_value_cutoff'],
                "analyse_repins": session['analyse_repins'],
                }
        )


        run_tree_task = len(strain_names) >= 4
        if run_tree_task:
            tree_job = RQJob.create(tree_task,
                                    depends_on=[rarefan_job],
                                    on_failure=on_failure,
                                    connection=app.redis,
                                    meta={'run_id': 'run_id', 'dbjob_id': dbjob.id},
                                    kwargs={
                                         "run_dir": session['tmpdir'],
                                         "treefile": session['treefile'],
                                        }
                                    )
        else:
            tree_job = RQJob.create(lambda x: None, depends_on=rarefan_job)

        zip_job = RQJob.create(zip_task,
                               depends_on=[rarefan_job, tree_job],
                               on_failure=on_failure,
                               meta={'run_id': 'run_id', 'dbjob_id': dbjob.id},
                               connection=app.redis,
                               kwargs={'run_dir':session['tmpdir']},
                                    )

        dbjob.stages = {'rarefan': {'redis_job_id': rarefan_job.get_id(),
                                    'results': {'returncode': None,
                                                'counts': {'rayts': None,
                                                           'seeds': None,
                                                           'repins': None,
                                                          },
                                                'data_sanity': {'rayts': None,
                                                           'seeds': None,
                                                           'repins': None,
                                                          },
                                    },
                                    'status': rarefan_job.get_status()
                                    },
                      'tree': {'redis_job_id': tree_job.get_id(),
                               'results': {'returncode': None, 'log': ""},
                               'status': tree_job.get_status()
                               },
                      'zip': {'redis_job_id': zip_job.get_id(),
                               'results': {'returncode': None, 'log': ""},
                               'status': zip_job.get_status()
                              }
                      }
        dbjob.save()

        # Enqueue the jobs
        app.queue.enqueue_job(rarefan_job)
        app.queue.enqueue_job(tree_job)
        app.queue.enqueue_job(zip_job)

        # Enqueue email job.
        # Emails should be sent out on overall success or failure. It must somehow be a monitoring job, or a callback.
        # Use a callback. on_success on zip task, on_failure on general. On failure must check dbjob if mail has already been sent.
        app.queue.enqueue(email_task,
                          dbjob,
                          depends_on=['rarefan_job',
                                      'tree_job',
                                      'zip_job',
                                      ]
                          )


        return redirect(url_for('results', run_id=run_id))

    return render_template(
                    'submit.html',
                    title='Submit',
                    submit_form=submit_form,
                    )


@app.route('/results', methods=['GET', 'POST'])
def results():

    args = request.args
    results_form = AnalysisForm()

    if 'run_id' in args.keys():

        run_id = args['run_id']

        dbjob = Job.objects.get_or_404(run_id=run_id)

        stages = dbjob.stages

        stati = {}

        for stage in stages.keys():
            redis_job_id = stages[stage]['redis_job_id']
            try:
                redis_job = rq.job.Job.fetch(redis_job_id, connection=app.redis)
                redis_job_status = redis_job.get_status()
            except NoSuchJobError:
                redis_job_status = "complete"
            except:
                redis_job_status = "failed"

            stati[stage] = redis_job_status

        if any([stage == "failed" for stage in stati.values()]):
            stati['overall'] = 'failed'
        elif all([stage in ['complete', 'finished'] for stage in stati.values()]):
            stati['overall'] = "complete"
        elif stati['rarefan'] == "queued":
            stati['overall'] = "queued"
        else:
            stati['overall'] = "running"

        # Update database.
        dbjob.update(set__stages__rarefan__status=stati['rarefan'])
        dbjob.update(set__stages__tree__status=stati['tree'])
        dbjob.update(set__stages__zip__status=stati['zip'])

        # Only show plots if more than 3 strains.
        render_plots = len(dbjob.setup.get('strain_names', [])) > 3
        return render_template('results.html',
                               title="Results for RAREFAN run {}".format(run_id),
                               results_form=results_form,
                               run_id=run_id,
                               stati=stati,
                               results_data=dbjob.stages['rarefan']['results'],
                               render_plots=render_plots,
                               )

    return render_template("results_query.html",
                           results_form=results_form,
                           title="Results",
                           )

@app.route('/files/<path:req_path>')
def files(req_path):
    """
    The 'files' route generates a navigable directory listing of a given run directory.
    @param req_path: The requested path to display.
    """
    uploads_dir = os.path.join(app.static_folder, 'uploads')
    nested_file_path = os.path.join(uploads_dir, req_path)
    splits = nested_file_path.split('/')
    uploads_idx = splits.index('uploads')
    run_id = splits[uploads_idx+1]


    if os.path.isdir(nested_file_path):
        item_list = os.listdir(nested_file_path)

        # Move directories to a separate list.
        dirs = [item_list.pop(i) for (i,d) in enumerate(item_list) if os.path.isdir(os.path.join(nested_file_path, d))]

        # Sort files and dirs.
        item_list.sort()
        dirs.sort()

        # Concat dirs and files.
        item_list = [i for i in item_list if not "stamp" in i]

        # Leading '/'
        if not req_path.startswith("/"):
            req_path = "/" + req_path

        # Remove trailing '/'
        if req_path.endswith('/'):
            req_path = req_path[:-1]
        logger.warning("Request dir is %s in (%s).", req_path, os.path.dirname(req_path))

        # Save the target for the 'back to results' link.
        tmp_dir = session.get('tmpdir', None)
        if tmp_dir is not None:
            run_id = os.path.basename(tmp_dir)

        try:
            back_link = url_for('results', run_id=run_id)
        except:
            back_link = url_for('results')

        link_to_parent = True
        # Only insert link to parent dir if not at top level.
        if os.path.dirname(req_path) == "/":
            link_to_parent = False

        return render_template('files.html',
                               req_path=req_path,
                               files=item_list,
                               dirs=dirs,
                               link_to_parent=link_to_parent,
                               back_link=back_link
                               )

    else:
        # Serve the file.
        return send_from_directory(*os.path.split(nested_file_path))


@app.route('/check_tasks', methods=['GET'])
def queue():
    args = request.args
    job_id = request.args['job_id']

    redis_job = rq.job.Job.fetch(job_id, connection=app.redis)

    redis_job.refresh()

    return("Job with id {} is finished: {}".format(redis_job.id, redis_job.is_finished))

@app.route('/manual', methods=['GET'])
def manual():
    return render_template('manual.html')

@app.route('/test_mail')
def test_mail():

    app.queue.enqueue(email_test)
    return redirect(url_for('index'))
