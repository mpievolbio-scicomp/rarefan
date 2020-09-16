from flask import render_template,\
                  request,\
                  session,\
                  redirect,\
                  url_for,\
                  abort,\
                  send_from_directory

from app.views import SubmitForm, RunForm, AnalysisForm
from app import app

import os, shutil, sys
import subprocess
import tempfile
#from werkzeug import secure_filename

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
            [seq.save(fname) for seq,fname in zip(seqs,fnames)]
            submit_form.reference_strain.choices = strain_names
            submit_form.query_rayt.choices += rayt_names
            submit_form.treefile.choices = ["None"] + tree_names
        
    if submit_form.go.data:
        tmpdir = session['tmpdir']
        session['outdir'] = outdir = os.path.join(tmpdir, 'out')
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
        query_rayt_fname = os.path.join(session['tmpdir'],session['query_rayt']+".faa")
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

        command = ['java',
                '-jar',
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
                ]
        print(" ".join(command))

        # Open log stream.
        with open(os.path.join(tmpdir, 'repinpop.log'), 'wb') as log:
            # Start subprocess as context manager.
            proc = subprocess.Popen(command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT)
            for line in iter(proc.stdout.readline, b''):
                sys.stdout.write(line.decode('ascii'))
                log.write(line)

            command = ["Rscript",
                       os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'displayREPINsAndRAYTs.R')),
                       session['outdir'],
                       treefile
                       ]

            proc = subprocess.Popen(command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT)

            proc.wait()

            # Zip results.
            command = ["zip",
                       "-rv",
                       os.path.split(session['tmpdir'])[-1]+"_out.zip",
                       'out'
                       ]

            proc = subprocess.Popen(command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT)

            proc.wait()

        os.chdir(oldwd)
        return redirect(url_for('results', run_id=os.path.basename(tmpdir)))

    return render_template(
                    'submit.html',
                    title='Submit',
                    submit_form=submit_form, 
                    )
          
@app.route('/results', methods=['GET', 'POST'])
def results():

    args = request.args
    print(args)

    results_form = AnalysisForm()

    if 'run_id' in args.keys():

        return render_template('results.html',
                title="Results",
                results_form=results_form,
                args=args,
                )
    else:
        return render_template('results_query.html',
                title="Results",
                results_form=results_form
                )


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

 
