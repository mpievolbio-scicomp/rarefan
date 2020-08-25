from flask import render_template, request
from app.views import SubmitForm
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
    # if upload_form.validate_on_submit():
    if True:
        tmpdir = tempfile.mkdtemp(suffix=None,
                    prefix=None,
                    dir=app.config["UPLOAD_DIR"]
                    )

        if upload_form.upload.data:
            print(upload_form.upload.data)
            # print("Upload validated")

            seqs = request.files.getlist(upload_form.sequences.name)

            if seqs:
                basenames = [seq.filename for seq in seqs]
                strain_names = [".".join(bn.split(".")[:-1]) for bn in basenames]
                fnames = [os.path.join(tmpdir, bn) for bn in basenames]
                [seq.save(fname) for seq,fname in zip(seqs,fnames)]
                upload_form.reference_strain.choices = strain_names

            else:
                print("ERROR: no sequence files found.")

            # return render_template('submit.html',
            #         title='Upload',
            #         upload_form=upload_form, 
            #         )

        if upload_form.go.data:
            print("Run")

            # oldwd = os.getcwd()
            # os.chdir(app.config["UPLOAD_DIR"])
            # command = ['java',
            #         '-jar',
            #         os.path.abspath(
            #             os.path.join(os.path.dirname(__file__),
            #                 '..',
            #                 'REPIN_ecology/REPIN_ecology/build/libs/REPIN_ecology.jar',
            #                 )
            #             ),
            #         tmpdir,
            #         request.form.get('reference_strain'),
            #         # '{0:d}'.format(request.form.get['min_nmer_occurence']),
            #         # '{0:d}'.format(upload_form.nmer_length),
            #         # os.path.join(tmpdir, upload_form.query_rayt),
            #         # os.path.join(tmpdir, 'tree.nwk'),
            #         # '{0:e}'.format(upload_form.e_value_cutoff)
            #         ]
            #
            # with subprocess.Popen(command,
            #                       stdout=subprocess.PIPE,
            #                       stderr=subprocess.STDOUT) as proc:
            #
            #
            #     print(proc.stdout.read())
            #
            # os.chdir(oldwd)

            # return render_template('submit.html',
            #         title='Submit',
            #         upload_form=upload_form, 
            #         reference_strain_form=reference_strain_form
            #         )
        

            # return render_template('submit.html',
            #         title='Submit',
            #         upload_form=upload_form, 
            #         )
    return render_template('submit.html',
            title='Submit',
            upload_form=upload_form, 
            )

@app.route('/run_repinpop', methods=['POST'])
def run_repinpop():
    print("Ajax Run")

    return "Nothing"
