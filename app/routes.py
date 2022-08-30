from flask import render_template
from flask import request
from flask import session
from flask import redirect
from flask import url_for
from flask import abort
from flask import send_from_directory
from flask import flash

from werkzeug.utils import secure_filename

import sys
import rq
from rq.job import Job as RQJob

from app.views import SubmitForm, AnalysisForm, UploadForm, ReturnToResultsForm, RunForm
from app import app, db
from app.models import Job as DBJob
from app.tasks.rarefan import rarefan_task
from app.tasks.tree import tree_task, empty_task
from app.tasks.rayt_phylo import alignment_task, phylogeny_task
from app.tasks.zip import zip_task
from app.tasks.email import email_task, email_test
from app.callbacks.callbacks import on_success, on_failure
from app.tasks import redis_tests

from Bio import SeqIO
import copy
import os
import shutil
import tempfile
import time


logger = app.logger


def get_status_code(run_id_path):
    """Legacy function (###REMOVEME)"""
    # Check if the run has finished.
    is_started = ".start.stamp" in os.listdir(run_id_path)
    is_java_finished = ".java.stamp" in os.listdir(run_id_path)
    is_zip_finished = ".zip.stamp" in os.listdir(run_id_path)

    status = is_started * 1 + is_java_finished * 10 + is_zip_finished * 100

    return status


def validate_fasta(filename):
    """
    Validate input from passed file as fasta formatted sequence data.

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
    """Route to homepage."""
    return render_template("index.html")


@app.route('/upload', methods=['GET', 'POST'])
def upload():
    """Upload files to server."""
    logger.debug("upload/%s", request.method)
    if request.method == 'POST':
        session['tmpdir'] = tempfile.mkdtemp(
            suffix=None,
            prefix="",
            dir=app.config["UPLOAD_DIR"]
        )

        seqs = [v for k, v in request.files.items() if k.startswith('file')]
        logger.info("Uploading %s.", str(seqs))

        dna_extensions = ['fn', 'fna', 'fastn', 'fas', 'fasta']
        aa_extensions = ['fa', 'faa']
        tree_extensions = ['nwk']

        # Generate filenames: Intermediate "." by "_" and secure.
        fnames = [seq.filename for seq in seqs]
        fnames = ["_".join(fname.split(".")[:-1]) + "." + fname.split(".")[-1] for fname in fnames]
        fnames = [os.path.join(session['tmpdir'], secure_filename(fname)) for fname in fnames]
        [seq.save(fname) for seq, fname in zip(seqs, fnames)]

        # Check if valid fasta files.
        for f in fnames:

            # Convert to unix filetype (see gh issue #9).
            with open(f, 'rb') as infile:
                content = infile.read()
            with open(f, 'wb') as output:
                for line in content.splitlines():
                    output.write(line + b'\n')

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

        for k, v in session.items():
            logger.debug("session[%s] = %s", k, str(v))

        return redirect(url_for('submit', _method='GET'))

    form = RunForm()
    session['tmpdir'] = None
    session['strain_names'] = None
    session['rayt_names'] = None
    session['tree_names'] = None
    session['outdir'] = None
    session['reference_strain'] = None
    session['query_rayt'] = None
    session['treefile'] = None
    session['min_nmer_occurrence'] = None
    session['distance_group_seeds'] = None
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
    """Submit job to server."""
    logger.debug("submit/%s", request.method)
    submit_form = SubmitForm()

    strain_names = session.get('strain_names')
    submit_form.reference_strain.choices.extend(strain_names)
    submit_form.query_rayt.choices.extend(session.get('rayt_names'))
    submit_form.treefile.choices.extend(["None"] + session.get('tree_names'))

    if submit_form.validate_on_submit():
        tmpdir = session['tmpdir']
        session['outdir'] = os.path.join(tmpdir, 'out')

        logger.debug("tmpdir: %s", tmpdir)
        logger.debug("tmpdir: %s", session['tmpdir'])

        # If there is already an outdir, this must be a rerun!
        old_run_id = None
        if os.path.isdir(session['outdir']):
            logger.warning("Rerun, creating new directory")
            rerun_tmpdir = tempfile.mkdtemp(
                suffix=None,
                prefix="",
                dir=app.config["UPLOAD_DIR"]
            )
            # Copy old run_dir
            old_run_id = os.path.basename(session['tmpdir'])
            shutil.copytree(session['tmpdir'], rerun_tmpdir, dirs_exist_ok=True)
            # Rm old results
            shutil.rmtree(os.path.join(rerun_tmpdir, 'out'))

            # Reset values.
            session['tmpdir'] = rerun_tmpdir
            session['outdir'] = os.path.join(rerun_tmpdir, 'out')

        os.mkdir(session['outdir'])

        session['reference_strain'] = request.form.get('reference_strain')
        session['query_rayt'] = request.form.get('query_rayt')
        session['min_nmer_occurrence'] = request.form.get('min_nmer_occurrence')
        treefile = request.form.get('treefile', None)
        run_id = os.path.basename(session['tmpdir'])

        session['treefile'] = treefile
        session['nmer_length'] = request.form.get('nmer_length')
        session['distance_group_seeds'] = request.form.get('distance_group_seeds', 15)
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

        # Create new Job instance.
        dbjob = DBJob(run_id=run_id,
                      stages={"rarefan": {"redis_job_id": None,
                                          "status": 'setup',
                                          "results": {"returncode": None,
                                                      "counts": {
                                                          "rayts": None,
                                                          "nmers": None,
                                                          "repins": None},
                                                      "data_sanity": {
                                                          "rayts": None,
                                                          "nmers": None,
                                                          "repins": None},
                                                      }
                                          },
                              "tree": {"redis_job_id": None,
                                       "status": 'setup',
                                       "results": {"returncode": None, "log": ""}},
                              "rayt_alignment": {
                                  "redis_job_id": None,
                                  "status": 'setup',
                                  "results": {"returncode": None,
                                              "log": None},
                              },
                              "rayt_phylogeny": {
                                  "redis_job_id": None,
                                  "status": 'setup',
                                  "results": {"returncode": None,
                                              "log": None},
                              }, 
                              "zip": {"redis_job_id": None,
                                      "status": 'setup',
                                      "results": {"returncode": None, "log": ""}}
                              },
                    setup=copy.deepcopy(session),
                    overall_status="setup",
                    notification_is_sent=False,
                    parent_run=old_run_id
                    )
        logger.debug("Constructed dbjob with job ID %s.", dbjob.run_id)
        logger.debug("Attempting to save dbjob in DB.")
        success = dbjob.save()
        logger.debug("Return code is %s", str(success))

        # If one of the server provided rayt files was selected, copy it to the working dir. In the dropdown menu,
        # the server provided rayts are listed without filename extension, so have to append that here.
        query_rayt_fname = os.path.join(session['tmpdir'], session['query_rayt'])
        if session['query_rayt'] in ['yafM_Ecoli', 'yafM_SBW25']:
            query_rayt_fname = query_rayt_fname + ".faa"
            src = os.path.join(app.static_folder, "rayts", session['query_rayt'] + ".faa")
            logger.debug("Copying rayt from %s to %s.", src, query_rayt_fname)
            shutil.copyfile(src, query_rayt_fname)

            if not os.path.isfile(query_rayt_fname):
                raise IOError("Copying %s to %s failed." % (src, query_rayt_fname))

        rarefan_job = RQJob.create(
            rarefan_task,
            connection=app.redis,
            on_success=on_success,
            on_failure=on_failure,
            timeout='24h',
            meta={'run_id': run_id, 'dbjob_id': dbjob.id, 'stage': 'rarefan'},
            kwargs={
                "tmpdir": session['tmpdir'],
                "outdir": session['outdir'],
                "reference_strain": session['reference_strain'],
                "min_nmer_occurrence": session['min_nmer_occurrence'],
                "nmer_length": session['nmer_length'],
                "distance_group_seeds": session.get('distance_group_seeds', 15),
                "query_rayt_fname": query_rayt_fname,
                "treefile": session['treefile'],
                "e_value_cutoff": session['e_value_cutoff'],
                "analyse_repins": session['analyse_repins'],
            }
        )
        logger.debug("Constructed rarefan job %s.", str(rarefan_job))

        run_tree_task = len(dbjob.setup['strain_names']) >= 4
        if run_tree_task:
            tree_job = RQJob.create(tree_task,
                                    depends_on=[rarefan_job],
                                    on_success=on_success,
                                    on_failure=on_failure,
                                    connection=app.redis,
                                    meta={'run_id': run_id, 'dbjob_id': dbjob.id, "stage": 'tree'},
                                    kwargs={
                                        "run_dir": session['tmpdir'],
                                        "treefile": session['treefile'],
                                    }
                                    )

        else:
            tree_job = RQJob.create(empty_task,
                                    on_success=on_success,
                                    on_failure=on_failure,
                                    depends_on=rarefan_job,
                                    meta={'run_id': run_id, 'dbjob_id': dbjob.id, "stage": 'tree'},
                                    connection=app.redis,
            )

        rayt_alignment_job = RQJob.create(
            alignment_task,
            depends_on=rarefan_job,
            meta={'run_id': run_id, "stage":'rayt_alignment'},
            connection=app.redis,
            kwargs={'run_id': run_id},
        )
        rayt_phylogeny_job = RQJob.create(
                    phylogeny_task,
                    depends_on=rayt_alignment_job,
                    meta={'run_id': run_id, "stage":'rayt_phylogeny'},
                    connection=app.redis,
                    kwargs={'run_id': run_id},
                )
    
        logger.debug("Constructed tree job %s.", str(tree_job))
        zip_job = RQJob.create(zip_task,
                               depends_on=[rarefan_job, tree_job],
                               on_success=on_success,
                               on_failure=on_failure,
                               meta={'run_id': run_id, 'dbjob_id': dbjob.id, "stage": 'zip'},
                               connection=app.redis,
                               kwargs={'run_dir': session['tmpdir']},
                               )
        logger.debug("Constructed zip job %s.", str(zip_job))

        email_job = RQJob.create(email_task,
                                 depends_on=['rarefan_job',
                                             'tree_job',
                                             'zip_job',
                                             ],
                                 meta={'run_id': run_id, 'dbjob_id': dbjob.id, "stage": 'email'},
                                 connection=app.redis,
                                 kwargs={'run_id': run_id},
                                 )

        logger.debug("Constructed email job %s.", str(email_job))

        # Enqueue the jobs
        app.queue.enqueue_job(rarefan_job)
        app.queue.enqueue_job(tree_job)
        app.queue.enqueue_job(rayt_alignment_job)
        app.queue.enqueue_job(rayt_phylogeny_job)
        app.queue.enqueue_job(zip_job)
        app.queue.enqueue_job(email_job)

        dbjob.update(set__stages__rarefan__redis_job_id=rarefan_job.id)
        dbjob.update(set__stages__tree__redis_job_id=tree_job.id)
        dbjob.update(set__stages__tree__redis_job_id=tree_job.id)
        dbjob.update(set__stages__tree__redis_job_id=rayt_alignment_job.id)
        dbjob.update(set__stages__tree__redis_job_id=rayt_phylogeny_job.id)
        dbjob.update(set__stages__zip__redis_job_id=zip_job.id)

        time.sleep(2)

        return redirect(url_for('results',
                                run_id=run_id,
                                _method='GET',)
                        )

    logger.debug("Form not validated, rendering submit template.")

    return render_template(
        'submit.html',
        title='Submit',
        submit_form=submit_form,
    )


@app.route('/results', methods=['GET', 'POST'])
def results():
    """Results summary page."""

    results_form = AnalysisForm()

    if request.method == 'GET':
        run_id = request.args.get('run_id', None)

    else:
        if results_form.validate_on_submit():
            run_id = request.form.get('run_id')

    if run_id is not None:
        logger.debug(run_id)

        try:
            dbjob = DBJob.objects.get(run_id=run_id)
        except:
            flash("Run {} was not found in our records. Please provide a valid run ID.".format(
                run_id))
            return render_template("results_query.html",
                                   results_form=results_form,
                                   title="Results",
                                   )

        # Update stage status by querying rq.
        dbjob.set_overall()

        # Only show plots if more than 3 strains.
        render_plots = len(dbjob.setup.get('strain_names', [])) > 3

        return render_template('results.html',
                               title="Results for RAREFAN run {}".format(run_id),
                               run_id=run_id,
                               job=dbjob,
                               render_plots=render_plots,
                               )

    return render_template("results_query.html",
                           results_form=results_form,
                           title="Results",
                           )


@app.route('/files/<path:req_path>')
def files(req_path):
    """
    Generate a navigable directory listing of a given run directory.

    :param req_path: The requested path to display.
    """
    uploads_dir = os.path.join(app.static_folder, 'uploads')
    nested_file_path = os.path.join(uploads_dir, req_path)
    splits = nested_file_path.split('/')
    uploads_idx = splits.index('uploads')
    run_id = splits[uploads_idx + 1]

    if os.path.isdir(nested_file_path):
        item_list = os.listdir(nested_file_path)

        # Move directories to a separate list.
        dirs = [item_list.pop(i) for (i, d) in enumerate(item_list) if os.path.isdir(os.path.join(nested_file_path, d))]

        # Sort files and dirs.
        item_list.sort()
        dirs.sort()

        # Concat dirs and files.
        item_list = [i for i in item_list if "stamp" not in i]

        # Leading '/'
        if not req_path.startswith("/"):
            req_path = "/" + req_path

        # Remove trailing '/'
        if req_path.endswith('/'):
            req_path = req_path[:-1]
        logger.info("Request dir is %s in (%s).", req_path, os.path.dirname(req_path))

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

    # Serve the file.
    return send_from_directory(*os.path.split(nested_file_path))


@app.route('/check_tasks', methods=['GET'])
def queue():
    """Queue a job on the redis job queue."""
    args = request.args
    job_id = request.args['job_id']

    redis_job = rq.job.Job.fetch(job_id, connection=app.redis)

    redis_job.refresh()

    return("Job with id {} is finished: {}".format(redis_job.id, redis_job.is_finished))


@app.route('/manual', methods=['GET'])
def manual():
    """Return the rarefan manual."""
    return render_template('manual.html')


@app.route('/rerun')
def rerun():
    """Rerun a job."""
    run_id = request.args['run_id']
    do_repins = request.args.get('do_repins', None)
    dbjob = DBJob.objects.get_or_404(run_id=run_id)
    logger.debug("Found job %s", str(dbjob.id))
    logger.debug("Job run_id = %s", str(dbjob.run_id))

    submit_form = SubmitForm()

    session['tmpdir'] = dbjob.setup.get('tmpdir')
    session['strain_names'] = dbjob.setup.get('strain_names')
    submit_form.reference_strain.choices.extend(session['strain_names'])
    session['reference_strain'] = dbjob.setup.get('reference_strain')
    submit_form.reference_strain.data = session['reference_strain']

    session['rayt_names'] = dbjob.setup.get('rayt_names')
    submit_form.query_rayt.choices.extend(session['rayt_names'])
    session['query_rayt'] = dbjob.setup.get('query_rayt')
    submit_form.query_rayt.data = session['query_rayt']

    session['tree_names'] = dbjob.setup.get('tree_names')
    submit_form.treefile.choices.extend(["None"] + session['tree_names'])

    session['treefile'] = dbjob.setup.get('treefile')
    submit_form.treefile.data = session['treefile']

    session['min_nmer_occurrence'] = dbjob.setup.get('min_nmer_occurrence')
    submit_form.min_nmer_occurrence.data = dbjob.setup.get('min_nmer_occurrence')
    session['nmer_length'] = dbjob.setup.get('nmer_length')
    submit_form.nmer_length.data = dbjob.setup.get('nmer_length')
    session['distance_group_seeds'] = dbjob.setup.get('distance_group_seeds', 15)
    submit_form.distance_group_seeds.data = dbjob.setup.get('distance_group_seeds', 15)
    session['analyse_repins'] = dbjob.setup.get('analyse_repins')
    submit_form.analyse_repins.data = dbjob.setup.get('analyse_repins')
    session['e_value_cutoff'] = dbjob.setup.get('e_value_cutoff')
    submit_form.e_value_cutoff.data = dbjob.setup.get('e_value_cutoff')
    session['email'] = dbjob.setup.get('email')
    submit_form.email.data = dbjob.setup.get('email')

    logger.debug("DO_REPINS? %s", do_repins)
    logger.debug("analyse_repins= %s", submit_form.analyse_repins.data)
    # Update do_repins if requested.
    if do_repins is not None:
        if do_repins in ['y', '1', 1, True]:
            submit_form.analyse_repins.data = session['analyse_repins'] = 'y'
        else:
            submit_form.analyse_repins.data = session['analyse_repins'] = None

    logger.debug('session = %s', str(session))

    return render_template(
        'submit.html',
        title='Submit',
        submit_form=submit_form,
    )


@app.route('/plot')
def plot():
    """ Redirect to the shiny app for the run id given via the request. """

    run_id = request.args['run_id']

    return redirect('http://rarefan.evolbio.mpg.de/shiny/analysis?run_id={}'.format(run_id))

@app.route('/test_task')
def test_task():
    job  = app.queue.enqueue(redis_tests.example, 10)
    logger.info(job.result)
    return redirect(url_for('index'))

@app.route('/test_mail')
def test_mail():
    success, message = email_test()
    # logger.debug("Attempting to send mail throug redis queue.")
    # job = app.queue.enqueue(email_test)
    # logger.debug(job)
    # time.sleep(3
               # )

    return message
