from flask import render_template, request, session
from app.views import SubmitForm, RunForm
from app import app

import os
import subprocess
import tempfile
#from werkzeug import secure_filename

@app.route('/')
def index():
    return render_template("index.html")

@app.route('/submit', methods=['GET', 'POST'])
def submit():
    upload_form = SubmitForm()
    tmpdir = tempfile.mkdtemp(suffix=None,
                prefix=None,
                dir=app.config["UPLOAD_DIR"]
                )

    print(upload_form.upload.data)
    print(upload_form.go.data)
    if upload_form.upload.data:
        # print("Upload validated")

        seqs = request.files.getlist(upload_form.sequences.name)

        if seqs:
            basenames = [seq.filename for seq in seqs]
            strain_names = [".".join(bn.split(".")[:-1]) for bn in basenames]
            fnames = [os.path.join(tmpdir, bn) for bn in basenames]
            [seq.save(fname) for seq,fname in zip(seqs,fnames)]
            upload_form.reference_strain.choices = strain_names
        
    if upload_form.go.data:
        session['tmpdir'] = tmpdir
        session['reference_strain'] = request.form.get('reference_strain')
        session['query_rayt'] = request.form.get('query_rayt')
        session['min_nmer_occurence'] = request.form.get('min_nmer_occurence')
        session['nmer_length'] = request.form.get('nmer_length')
        session['e_value_cutoff'] = request.form.get('e_value_cutoff')

        print(session["tmpdir"])
        print(session["reference_strain"])
        print(session["nmer_length"])

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
                tmpdir,
                session['reference_strain'],
                '{0:s}'.format(session['min_nmer_occurence']),
                '{0:s}'.format(session['nmer_length']),
                os.path.join(tmpdir, session['query_rayt']),
                os.path.join(tmpdir, 'tree.nwk'),
                '{0:s}'.format(session['e_value_cutoff'])
                ]
        print(command)

        with subprocess.Popen(command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT) as proc:
            print(proc.stdout.read())

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
