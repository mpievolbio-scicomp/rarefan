from flask import render_template
from flask import request
from flask import session
from flask import redirect
from flask import url_for
from flask import abort
from flask import send_from_directory
from flask import flash

from werkzeug.utils import secure_filename
from app.views import SubmitForm, AnalysisForm, UploadForm, ReturnToResultsForm, RunForm
from app import app

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
logging.basicConfig(level=logging.INFO)

def get_logger():
    logger = logging.getLogger(__name__)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(module)s: %(message)s')
    timestamp = datetime.datetime.now().strftime(format="%Y%m%d-%H%M%S")
    handler = logging.FileHandler("/tmp/rarefan_{}.log".format(timestamp))
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)

    return logger

logger = get_logger()

logger.warning("RAREFAN")


def get_no_cpus():
    return os.cpu_count()


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
        session['min_nmer_occurence'] = None
        session['treefile'] = None
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
        os.mkdir(session['outdir'])
        session['reference_strain'] = request.form.get('reference_strain')
        session['query_rayt'] = request.form.get('query_rayt')
        session['min_nmer_occurence'] = request.form.get('min_nmer_occurence')
        treefile = request.form.get('treefile', None)
        run_id = os.path.basename(session['tmpdir'])

        if len(strain_names) >= 4:
            if treefile == "None":
                run_andi_clustdist = True
                treefile = "tmptree.nwk"
            else:
                run_andi_clustdist = False
        else:
            treefile = "tmptree.nwk"
            run_andi_clustdist = False

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

        # If one of the server provided rayt files was selected,  copy it to the working dir. In the dropdown menu,
        # the server provided rayts are listed without filename extension, so have to append that here.
        query_rayt_fname = os.path.join(session['tmpdir'], session['query_rayt'])
        if session['query_rayt'] in ['yafM_Ecoli', 'yafM_SBW25']:
            query_rayt_fname = query_rayt_fname+".faa"
            src = os.path.join(app.static_folder, "rayts", session['query_rayt']+".faa")
            logging.debug("Copying rayt from %s to %s.", src, query_rayt_fname)
            shutil.copyfile(src, query_rayt_fname)
            
            if not os.path.isfile(query_rayt_fname):
                raise IOError("Copying %s to %s failed." % (src, query_rayt_fname))

        # Copy R scripts
        shutil.copyfile(os.path.join(os.path.dirname(__file__),
                                     "..",
                                     'shinyapps',
                                     'analysis',
                                      "analysis.R"
                                     ),
                        os.path.join(session['outdir'],
                                     'analysis.R'
                                     )
                        )
        shutil.copyfile(os.path.join(os.path.dirname(__file__),
                                     "..",
                                     'shinyapps',
                                     'analysis',
                                     "run_analysis.R"
                                     ),
                        os.path.join(session['outdir'],
                                     'run_analysis.R'
                                     )
                        )

        oldwd = os.getcwd()
        os.chdir(tmpdir)

        start_stamp = os.path.join(session['tmpdir'], '.start.stamp')

        mcl_threads = max(get_no_cpus()//4, 1)

        logging.info("Detected %d CPUs, will utilize %d.", 2*mcl_threads, mcl_threads)
        java_command = " ".join(['java',
                                     '-Dcom.sun.management.jmxremote',
                                      '-Dcom.sun.management.jmxremote.port=9010',
                                      '-Dcom.sun.management.jmxremote.local.only=true',
                                      '-Dcom.sun.management.jmxremote.authenticate=false',
                                      '-Dcom.sun.management.jmxremote.ssl=false',
                                     '-jar',
                                     '-Xmx10g',
                                     os.path.abspath(
                                     os.path.join(os.path.dirname(app.root_path),
                                     'REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar',
                                     )
                                     ),
                                     session['tmpdir'],
                                     session['outdir'],
                                     session['reference_strain'],
                                     '{0:s}'.format(session['min_nmer_occurence']),
                                     '{0:s}'.format(session['nmer_length']),
                                     query_rayt_fname,
                                     treefile,
                                     '{0:s}'.format(session['e_value_cutoff']),
                                     {"y": "true", None: "false"}[session['analyse_repins']],
                                     '{0:d}'.format(mcl_threads),
                                        ]
                                )

        logging.info("Java command: %s", java_command)
        java_stamp = os.path.join(session['tmpdir'], '.java.stamp')

        R_command = " ".join(["Rscript",
                                  'displayREPINsAndRAYTs.R',
                                  session['outdir'],
                                  treefile
                                     ])
        logging.info("R command: %s", R_command)
        R_stamp = os.path.join(session['tmpdir'], '.R.stamp')

        if run_andi_clustdist:
            andi_inputs = [os.path.join(session['tmpdir'], f) for f in os.listdir() if f.split(".")[-1] in ["fas", "fna", "fn", "fasta", "fastn"]]
            distfile = "".join(session['treefile'].split('.')[:-1])+'.dist'
            distfile = os.path.join(session['outdir'], os.path.basename(distfile))
            andi_command = "andi -j {} > {}".format(" ".join(andi_inputs), distfile)
            logging.info("andi command: %s", andi_command)
            andi_stamp = os.path.join(session['tmpdir'], '.andi.stamp')

            clustdist_command = "clustDist {} > {}".format(distfile, os.path.join(session['outdir'],treefile))

        else:
            andi_command = "echo 'Not running andi.'"
            logging.info("andi command: %s", andi_command)
            andi_stamp = os.path.join(session['tmpdir'], '.andi.stamp')

            if len(session.get('strain_names')) > 3:
                clustdist_command = "ln -s {} {}".format(os.path.join(session['tmpdir'], treefile),
                                                         os.path.join(session['outdir'], 'tmptree.nwk'))

            else:
                clustdist_command = "echo 'Not running clustDist.'"

        logging.info("clustdist command: %s", clustdist_command)
        clustdist_stamp = os.path.join(session['tmpdir'], '.clustdist.stamp')

        # Zip results.
        zip_command = " ".join(["zip",
                                "-rv",
                                os.path.split(session['tmpdir'])[-1] + "_out.zip",
                                'out'
                                ]
                               )
        zip_stamp = os.path.join(session['tmpdir'], '.zip.stamp')
        logging.info("zip command: %s", zip_command)

        email_command = get_email_command(session)
        email_stamp = os.path.join(session['tmpdir'], '.email.stamp')

        command_lines = [
            "touch {} &&".format(start_stamp),
            "{} && touch {}".format(java_command, java_stamp),
            "{} && touch {}".format(andi_command, andi_stamp),
            "{} && touch {}".format(clustdist_command, clustdist_stamp),
            "{} && touch {}".format(zip_command, zip_stamp),
            "{} && touch {}".format(email_command, email_stamp)
        ]

        with open(os.path.join(tmpdir,'job.sh'), 'w') as fp:
            fp.write(r"#! /bin/bash")
            fp.write('\n')
            fp.write("export LD_LIBRARY_PATH={}".format(os.environ["LD_LIBRARY_PATH"]))
            fp.write('\n')
            for line in command_lines:
                fp.write(line)
                fp.write('\n')
            fp.write('\n')

        os.chmod('job.sh', stat.S_IRWXU)

        # Write batch script to submit the job.
        with open(os.path.join(tmpdir, 'batch.sh'), 'w') as fp:
            fp.write(r"#! /bin/bash")
            fp.write('\n')
            fp.write('echo "./job.sh > out/rarefan.log 2>&1" | batch')
            fp.write('\n')

        os.chmod('batch.sh', stat.S_IRWXU)

        shell_command = os.path.join(tmpdir, 'batch.sh')
        proc = subprocess.Popen(shlex.split(shell_command), shell=False)

        os.chdir(oldwd)

        return redirect(url_for('results', run_id=run_id))

    return render_template(
                    'submit.html',
                    title='Submit',
                    submit_form=submit_form,
                    )

def get_email_command(session):

    # Aggregate the run path.
    run_id_path = session["tmpdir"]
    run_id = os.path.basename(run_id_path)

    # List of recipients. Always send to rarefan.
    recipients = ["rarefan@evolbio.mpg.de"]
    if session['email'] is None or session['email'] == "":
        logging.debug("No email set.")

    recipients.append(session["email"])

    email_subject = "Your RAREFAN run {0:s} has completed.".format(run_id)
    email_body = """Hallo,
your job on rarefan.evolbio.mpg.de with ID {0:s} has completed.
You can browse and download the results at this link:
http://rarefan.evolbio.mpg.de/results?run_id={0:s}.

In case of problems, please reply to this email and leave the email subject as is.

Thank you for using RAREFAN.

http://rarefan.evolbio.mpg.de
""".format(run_id)

    email_command = 'printf "Subject: {0:s}\n\n{1:s}" | msmtp {2:s} >> {3:s}'.format(
        email_subject,
        email_body,
        " ".join(recipients),
        os.path.join(
            run_id_path,
            'out',
            'rarefan.log'
        )
    )
    return email_command



@app.route('/results', methods=['GET', 'POST'])
def results():

    args = request.args
    results_form = AnalysisForm()

    if 'run_id' in args.keys():
        
        run_id = args['run_id']
        # Check if this is a valid run id.

        run_id_path = os.path.join(app.static_folder, "uploads", run_id)
        is_valid_run_id = os.path.isdir(run_id_path)

        status = get_status_code(run_id_path)

        if status < 1:
            flash("Your job {} is queued, please wait for page to refresh.".format(run_id))
        elif status == 1:
            flash("Your job {} is running, please wait for page to refresh.".format(run_id))
        elif status == 11:
            flash("Your job {} is finished. Preparing run files for download".format(run_id))
        elif status == 101:
            flash("Your job {} has failed. Preparing run files for download.".format(run_id))
        elif status == 111:
            flash("Your job {} has finished. Results and download links below.".format(run_id))
        elif status in [10, 100]:
            flash("Your job {} has failed. Please inspect the run files and resubmit your data.".format(run_id))
        else:
            flash("Your job {} has failed with an unexpected failure.".format(run_id))

        # Only show plots if more than 3 strains uploaded.
        strain_names = session.get('strain_names', None)
        if strain_names is None:
            render_plots = True
        else:
            render_plots = len(strain_names) > 3
        return render_template('results.html',
                               title="Run {} results".format(run_id),
                               results_form=results_form,
                               run_id=run_id,
                               status=status,
                               render_plots=render_plots,
                               )

    flash("Not a valid run ID.")

    return render_template("results_query.html",
                       results_form=results_form,
                           title="Results")

@app.route('/files/<path:req_path>')
def files(req_path):
    """"""
    uploads_dir = os.path.join(app.static_folder, 'uploads')
    nested_file_path = os.path.join(uploads_dir, req_path)
    #
    # if os.path.realpath(nestedFilePath) != nestedFilePath:
    #     return "no directory traversal please."

    if os.path.isdir(nested_file_path):
        item_list = os.listdir(nested_file_path)

        # Move directories to a separate list.
        dirs = [item_list.pop(i) for (i,d) in enumerate(item_list) if os.path.isdir(os.path.join(nested_file_path, d))]

        # Sort files and dirs.
        item_list.sort()
        dirs.sort()

        # Concat dirs and files.
        item_list = [i for i in item_list if not "stamp" in i]

        if not req_path.startswith("/"):
            req_path = "/" + req_path
        if req_path.endswith('/'):
            req_path = req_path[:-1]
        logger.warning("Request dir is %s in (%s).", req_path, os.path.dirname(req_path))

        tmp_dir = session.get('tmpdir', None)
        if tmp_dir is not None:
            back_link = url_for('results', run_id=os.path.basename(tmp_dir))
        else:
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

@app.route('/manual', methods=['GET'])
def manual():
    return render_template('manual.html')