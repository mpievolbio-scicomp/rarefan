from flask import render_template, request, session
from app.views import SubmitForm, RunForm
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
    upload_form = SubmitForm()

    

    

    print(upload_form.upload.data)
    print(upload_form.go.data)
    if upload_form.upload.data:
        seqs = request.files.getlist(upload_form.sequences.name)
        
        session['tmpdir'] = tempfile.mkdtemp(suffix=None,
                              prefix=None,
                              dir=app.config["UPLOAD_DIR"]
                              )
        if seqs:
            basenames = [seq.filename for seq in seqs]
            strain_names = [".".join(bn.split(".")[:-1]) for bn in basenames if bn.split(".")[-1] == "fas"]
            rayt_names = [".".join(bn.split(".")[:-1]) for bn in basenames if bn.split(".")[-1] == "faa"]
            tree_names = [bn for bn in basenames if bn.split(".")[-1] == "nwk"]
            fnames = [os.path.join(session['tmpdir'], bn) for bn in basenames]
            [seq.save(fname) for seq,fname in zip(seqs,fnames)]
            upload_form.reference_strain.choices = strain_names
            upload_form.query_rayt.choices += rayt_names
            upload_form.treefile.choices = ["None"] + tree_names
        
    if upload_form.go.data:
        tmpdir = session['tmpdir']
        session['reference_strain'] = request.form.get('reference_strain')
        session['query_rayt'] = request.form.get('query_rayt')
        session['min_nmer_occurence'] = request.form.get('min_nmer_occurence')
        treefile = request.form.get('treefile')
        if treefile == "None":
            treefile = "tmptree.nwk"
        session['treefile'] = treefile
        session['nmer_length'] = request.form.get('nmer_length')
        session['e_value_cutoff'] = request.form.get('e_value_cutoff')

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
        os.chdir(app.config["UPLOAD_DIR"])

        command = ['java',
                '-jar',
                os.path.abspath(
                    os.path.join(os.path.dirname(__file__),
                        '..',
                        'REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar',
                        )
                    ),
                session['tmpdir'],
                session['reference_strain']+".fas",
                '{0:s}'.format(session['min_nmer_occurence']),
                '{0:s}'.format(session['nmer_length']),
                query_rayt_fname,
                os.path.join(tmpdir, treefile),
                '{0:s}'.format(session['e_value_cutoff'])
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

        os.chdir(oldwd)

    return render_template('submit.html',
                    title='Submit',
                    upload_form=upload_form, 
                    )
          
@app.route('/run_repinpop', methods=['GET', 'POST'])
def run_repinpop():
    print(session["tmpdir"])
    print(session["reference_strain"])
    print(session["nmer_length"])
    
    return "`'/-"


        #
    # return "Running REPINPop. Stay tuned...."
