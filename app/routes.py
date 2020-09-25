from flask import render_template,\
                  request,\
                  session,\
                  redirect,\
                  url_for,\
                  abort,\
                  send_from_directory

from app.views import SubmitForm, AnalysisForm
from app import app

import os
import shutil
import stat
import subprocess
import tempfile


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/submit', methods=['GET', 'POST'])
def submit():
    submit_form = SubmitForm()

    if submit_form.upload.data:
        seqs = request.files.getlist(submit_form.sequences.name)
        
        session['tmpdir'] = tempfile.mkdtemp(
                              suffix=None,
                              prefix="",
                              dir=app.config["UPLOAD_DIR"]
                              )
        if seqs:
            basenames = [seq.filename for seq in seqs]
            strain_names = [".".join(bn.split(".")[:-1]) for bn in basenames if bn.split(".")[-1] == "fas"]
            rayt_names = [".".join(bn.split(".")[:-1]) for bn in basenames if bn.split(".")[-1] == "faa"]
            tree_names = [bn for bn in basenames if bn.split(".")[-1] == "nwk"]
            fnames = [os.path.join(session['tmpdir'], bn) for bn in basenames]
            [seq.save(fname) for seq, fname in zip(seqs, fnames)]
            submit_form.reference_strain.choices = strain_names
            submit_form.query_rayt.choices += rayt_names
            submit_form.treefile.choices = ["None"] + tree_names
        
    if submit_form.go.data:
        tmpdir = session['tmpdir']
        session['outdir'] = os.path.join(tmpdir, 'out')
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
        
        oldwd = os.getcwd()
        os.chdir(tmpdir)

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
                   os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'displayREPINsAndRAYTs.R')),
                   session['outdir'],
                   treefile
                   ])
        R_stamp = os.path.join(session['tmpdir'], '.R.stamp')

        # Zip results.
        zip_command = " ".join(["zip",
                   "-rv",
                   os.path.split(session['tmpdir'])[-1]+"_out.zip",
                   'out'
                   ])

        zip_stamp = os.path.join(session['tmpdir'], '.zip.stamp')
        
        command_lines = [java_command+" && ",
                        "touch {} && ".format(java_stamp),
                        R_command+" && ",
                        "touch {} && ".format(R_stamp),
                        zip_command+" && ",
                        "touch {}".format(zip_stamp)
                        ]
        
        with open(os.path.join(tmpdir,'job.sh'), 'w') as fp:
            fp.write(r"#! /bin/sh") 
            fp.write('\n')
            for line in command_lines:
                fp.write(line)
            fp.write('\n')

        os.chmod('job.sh', stat.S_IRWXU )
        shell_command = [os.path.join(tmpdir, 'job.sh'), '> {}'.format(os.path.join(tmpdir, 'rarefan.log')), '2>&1']
        print(" ".join(shell_command))
        proc = subprocess.Popen(shell_command)

        os.chdir(oldwd)
        
        return redirect(url_for('results', run_id=os.path.basename(session['tmpdir'])))

    return render_template(
                    'submit.html',
                    title='Submit',
                    submit_form=submit_form, 
                    )
          
@app.route('/results', methods=['GET', 'POST'])
def results():

    args = request.args
    
    results_form = AnalysisForm()

    is_valid_run_id = False
    is_java_finished = False
    is_R_finished = False
    is_zip_finished = False

    if 'run_id' in args.keys():
        
        run_id = args['run_id']
        # Check if this is a valid run id.

        run_id_path = os.path.join(app.static_folder, "uploads", run_id)
        is_valid_run_id = os.path.isdir(run_id_path)
         
        if is_valid_run_id:
            # Check if the run has finished.
            is_java_finished = ".java.stamp" in os.listdir(run_id_path)
            is_R_finished = ".R.stamp" in os.listdir(run_id_path)
            is_zip_finished = ".zip.stamp" in os.listdir(run_id_path)
        
        print(run_id, is_valid_run_id, is_java_finished, is_R_finished, is_zip_finished)                        
        return render_template('results.html',
                title="Results",
                results_form=results_form,
                run_id=run_id,
                is_valid_run_id=is_valid_run_id,
                is_java_finished=is_java_finished,
                is_R_finished=is_R_finished,
                is_zip_finished=is_zip_finished,
                )
    return render_template("results_query.html",
                           results_form=results_form)

@app.route('/files/<path:req_path>')
def files(req_path):

    base_dir = os.path.join(app.static_folder, 'uploads')
    abs_path = os.path.join(base_dir, req_path)
    print(base_dir, abs_path)

    if not os.path.exists(abs_path):
        print("{} not found".format(abs_path))
        return abort(404)

    if os.path.isfile(abs_path):
        print("serving file")
        return send_from_directory(base_dir, req_path)

    fnames = os.listdir(abs_path)
    print(fnames)

    return render_template('files.html', files=fnames)

 
