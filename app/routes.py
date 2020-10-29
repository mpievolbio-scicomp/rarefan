from flask import render_template
from flask import request
from flask import session
from flask import redirect
from flask import url_for
from flask import abort
from flask import send_from_directory
from flask import flash

from werkzeug.utils import secure_filename
import re
from wtforms.validators import ValidationError
from app.views import SubmitForm, AnalysisForm, UploadForm
from app import app

import os
import shlex
import shutil
import stat
import subprocess
import tempfile
from threading import Thread

from Bio import SeqIO

def validate_fasta(filename):
    """ Validates input from passed file as fasta formatted sequence data.

    :param filename: The filename of the file to validate.

    """
    with open(filename, 'r') as fp:
        fasta = SeqIO.parse(fp, "fasta")
        is_fasta = any(fasta)

        return is_fasta


@app.route('/')
def index():
    return render_template("index.html")

@app.route('/upload', methods=['GET', 'POST'])
def upload():
    upload_form = UploadForm()

    if upload_form.validate_on_submit():
        seqs = request.files.getlist(upload_form.sequences.name)
        print(seqs)


        session['tmpdir'] = tempfile.mkdtemp(
            suffix=None,
            prefix="",
            dir=app.config["UPLOAD_DIR"]
        )

        dna_extensions = ['fn', 'fastn', 'fas', 'fasta']
        aa_extensions = ['fa', 'faa']
        tree_extensions = ['nwk']
        fnames = [os.path.join(session['tmpdir'], secure_filename(seq.filename)) for seq in seqs]
        [seq.save(fname) for seq, fname in zip(seqs, fnames)]

        # Check if valid fasta files.
        for f in fnames:
            if f.split('.')[-1] in tree_extensions:
                continue
            is_fasta = validate_fasta(f)
            if not validate_fasta(f):
                flash("{} is neither a valid fasta file nor a tree file and will be ignored.".format(os.path.basename(f)))
                fnames.remove(f)
                os.remove(f)

        basenames = [os.path.basename(f) for f in fnames]
        strain_names = [".".join(bn.split(".")[:-1]) for bn in basenames if bn.split(".")[-1] in dna_extensions]
        rayt_names = [".".join(bn.split(".")[:-1]) for bn in basenames if bn.split(".")[-1] in aa_extensions]
        tree_names = [bn for bn in basenames if bn.split(".")[-1] in tree_extensions]

        session['strain_names'] = strain_names
        session['rayt_names'] = rayt_names
        session['tree_names'] = tree_names

        return redirect(url_for('submit'))

    return render_template(
        'upload.html',
        title="Upload sequences",
        upload_form=upload_form,
    )

@app.route('/submit', methods=['GET', 'POST'])
def submit():

    submit_form = SubmitForm()
    submit_form.reference_strain.choices.extend(session.get('strain_names'))
    submit_form.query_rayt.choices.extend(session.get('rayt_names'))
    submit_form.treefile.choices.extend(["None"] + session.get('tree_names'))

    if submit_form.validate_on_submit():
        tmpdir = session['tmpdir']
        session['outdir'] = os.path.join(tmpdir, 'out')
        os.mkdir(session['outdir'])
        session['reference_strain'] = request.form.get('reference_strain')
        session['query_rayt'] = request.form.get('query_rayt')
        session['min_nmer_occurence'] = request.form.get('min_nmer_occurence')
        treefile = request.form.get('treefile')
        if treefile == "None":
            treefile = "tmptree.nwk"
        session['treefile'] = treefile
        session['nmer_length'] = request.form.get('nmer_length')
        session['e_value_cutoff'] = request.form.get('e_value_cutoff')
        session['analyse_repins'] = request.form.get('analyse_repins')
        session['email'] = request.form.get('email', None)

        # copy query rayt to working dir
        query_rayt_fname = os.path.join(session['tmpdir'], session['query_rayt']+".faa")
        if session['query_rayt'] in ['yafM_Ecoli', 'yafM_SBW25']:
                                                                     src=os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                                     "..",
                                                                     'data',
                                                                     session['query_rayt']+".faa"
                                                                     )
                                                                     )
                                                                     shutil.copyfile(src, query_rayt_fname)

        # Copy R script
        shutil.copyfile(os.path.join(os.path.dirname(__file__),
                                     "..", "displayREPINsAndRAYTs.R"
                                     ),
                        os.path.join(session['tmpdir'],
                                     'displayREPINsAndRAYTs.R'
                                     )
                        )

        oldwd = os.getcwd()
        os.chdir(tmpdir)

        start_stamp = os.path.join(session['tmpdir'], '.start.stamp')

        java_command = " ".join(['java',
                                     '-jar',
                                     '-Xmx14g',
                                     os.path.abspath(
                                     os.path.join(os.path.dirname(app.root_path),
                                     'REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar',
                                     )
                                     ),
                                     session['tmpdir'],
                                     session['outdir'],
                                     session['reference_strain']+".fas",
                                     '{0:s}'.format(session['min_nmer_occurence']),
                                     '{0:s}'.format(session['nmer_length']),
                                     query_rayt_fname,
                                     treefile,
                                     '{0:s}'.format(session['e_value_cutoff']),
                                     {"y": "true", None: "false"}[session['analyse_repins']],
                                        ])

        java_stamp = os.path.join(session['tmpdir'], '.java.stamp')

        R_command = " ".join(["Rscript",
                                  'displayREPINsAndRAYTs.R',
                                  session['outdir'],
                                  treefile
                                     ])
        R_stamp = os.path.join(session['tmpdir'], '.R.stamp')

        andi_inputs = [os.path.join(session['tmpdir'], f) for f in os.listdir() if f.split(".")[-1] in ["fas", "fna"]]
        distfile = "".join(session['treefile'].split('.')[:-1])+'.dist'
        distfile = os.path.join(session['outdir'], os.path.basename(distfile))
        andi_command = "andi -j {} > {}".format(" ".join(andi_inputs), distfile)
        andi_stamp = os.path.join(session['tmpdir'], '.andi.stamp')

        clustdist_command = "clustDist {} > {}".format(distfile, os.path.join(session['outdir'],treefile))
        clustdist_stamp = os.path.join(session['tmpdir'], '.clustdist.stamp')

        # Zip results.
        zip_command = " ".join(["zip",
                                "-rv",
                                os.path.split(session['tmpdir'])[-1] + "_out.zip",
                                'out'
                                ]
                               )
        zip_stamp = os.path.join(session['tmpdir'], '.zip.stamp')

        command_lines = [
            "touch {} &&".format(start_stamp),
            "{} && touch {}".format(java_command, java_stamp),
            "{} && touch {}".format(andi_command, andi_stamp),
            "{} && touch {}".format(clustdist_command, clustdist_stamp),
            "{} && touch {}".format(R_command, R_stamp),
            "{} && touch {}".format(zip_command, zip_stamp)
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

        os.chmod('job.sh', stat.S_IRWXU )

        # Write batch script to submit the job.
        with open(os.path.join(tmpdir,'batch.sh'), 'w') as fp:
            fp.write(r"#! /bin/bash")
            fp.write('\n')
            fp.write('echo "./job.sh > out/rarefan.log 2>&1" | batch')
            fp.write('\n')

        os.chmod('batch.sh', stat.S_IRWXU )

        shell_command = os.path.join(tmpdir, 'batch.sh')
        proc = subprocess.Popen(shlex.split(shell_command), shell=False)

        os.chdir(oldwd)

        return redirect(url_for('results', run_id=os.path.basename(session['tmpdir'])))

    return render_template(
                    'submit.html',
                    title='Submit',
                    submit_form=submit_form,
                    )

def send_email(run_id, status_code, recipient):

    # Aggregate the run path.
    run_id_path = os.path.join(app.static_folder, "uploads", run_id)

    # Check if email notification was requested.
    if recipient is None or recipient == "":
        return

    # Check if an email has already been sent.
    if '.email.stamp' in os.listdir(run_id_path):
        return

    recipients = [recipient]

    # Job failed.
    if status_code in [101, 1011, 1001]:
        email_subject = "Your RAREFAN run {0:s} has failed.".format(os.path.basename(run_id_path))
        email_body = """Hallo,
your job on rarefan.evolbio.mpg.de with ID {0:s} has failed.
You can browse and download the run files at this link:
http://rarefan.evolbio.mpg.de/results?run_id={0:s}.

Pay attention to the log file under out/rarefan.log as it may provide further information about the failure.
Please feel free to seek our support at mailto:computing.evolbio.mpg.de.

Thank you for using RAREFAN. We hope to see you soon again.

Kind regards,

RAREFAN.

http://rarefan.evolbio.mpg.de
""".format(os.path.basename(run_id_path))

        # Include admin as recipient.
        recipients.append('computing@evolbio.mpg.de')

    # Job success.
    elif status_code == 1111:
        email_subject = "Your RAREFAN run {0:s} has finished.".format(os.path.basename(run_id_path))
        email_body = """Hallo,
your job on rarefan.evolbio.mpg.de with ID {0:s} has finished.
You can browse and download the results at this link:
http://rarefan.evolbio.mpg.de/results?run_id={0:s}.

Thank you for using RAREFAN. We hope to see you soon again.

Kind regards,

RAREFAN.

http://rarefan.evolbio.mpg.de
""".format(os.path.basename(run_id_path))

    # All other cases (job still running or queued).
    else:
        return

    # Send mail to all recipients.
    for recipient in recipients:
        email_command = 'printf "Subject: {0:s}\n\n{1:s}" | msmtp -v {2:s} > {3:s}'.format(email_subject,
                                                                                email_body, recipient, os.path.join(run_id_path, 'out', 'rarefan.log'))
        proc = subprocess.Popen(email_command, shell=True)
    #
    # Generate email stamp.
    proc = subprocess.Popen(shlex.split("touch {}/.email.stamp".format(run_id_path)))


@app.route('/results', methods=['GET', 'POST'])
def results():

    args = request.args
    
    results_form = AnalysisForm()

    if 'run_id' in args.keys():
        
        run_id = args['run_id']
        # Check if this is a valid run id.

        run_id_path = os.path.join(app.static_folder, "uploads", run_id)
        is_valid_run_id = os.path.isdir(run_id_path)
         
        if is_valid_run_id:
            # Check if the run has finished.
            is_started = ".start.stamp" in os.listdir(run_id_path)
            is_java_finished = ".java.stamp" in os.listdir(run_id_path)
            is_R_finished = ".R.stamp" in os.listdir(run_id_path)
            is_zip_finished = ".zip.stamp" in os.listdir(run_id_path)

            status = is_started*1 + is_java_finished*10 + is_R_finished*100 + is_zip_finished*1000
            # flash("DEBUG: Status={}".format(status))

            if status < 1:
                flash("Your job {} is queued, please wait for page to refresh.".format(run_id))
            elif status == 1:
                flash("Your job {} is running, please wait for page to refresh.".format(run_id))
            elif status == 11:
                flash("Your job {} has finished, postprocessing.".format(run_id))
            elif status == 111:
                flash("Your job {} and postprocessing have finished. Preparing run files for download".format(run_id))
            elif status == 101:
                flash("Your job {} has failed. Preparing run files for download.".format(run_id))
            elif status == 1111:
                flash("Your job {} has finished. Results and download links below.".format(run_id))
            elif status == 1011:
                flash("Your job {} has finished but postprocessing failed. Download files below.".format(run_id))
            elif status == 1001:
                flash("Your job {} has failed. Please inspect the run files and resubmit your data.".format(run_id))
            else:
                flash("Your job {} has failed with an unexpected failure.".format(run_id))

            send_email(run_id, status, session['email'])
            return render_template('results.html',
                                   title="Run {} results".format(run_id),
                                   results_form=results_form,
                                   run_id=run_id,
                                   status=status
                                   )

        else:
            flash("Not a valid run ID.")

    return render_template("results_query.html",
                       results_form=results_form,
                           title="Results")

@app.route('/files/<path:req_path>')
def files(req_path):
    """"""
    """ Only a stub, file listing will be taken care of by AutoIndex."""
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

        # Prepend the parent dir.
        dirs.insert(0, '..')
        if not req_path.startswith("/"):
            req_path = "/" + req_path
        if req_path.endswith('/'):
            req_path = req_path[:-1]
        return render_template('files.html', req_path=req_path, files=item_list, dirs=dirs)
    else:
        # Serve the file.
        return send_from_directory(*os.path.split(nested_file_path))
