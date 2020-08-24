from flask import render_template, request
from app.views import SubmitForm, SequenceUploadForm, ReferenceStrainForm
from app import app

import os
#from werkzeug import secure_filename

@app.route('/')
def index():
    return render_template("index.html")


@app.route('/submit', methods=['GET', 'POST'])
def submit():
    sequence_upload_form = SequenceUploadForm()
    reference_strain_form = ReferenceStrainForm()

    if sequence_upload_form.validate_on_submit():
        seqs = request.files.getlist(sequence_upload_form.sequences.name)
        reference_strain_form.choices = seqs

        if seqs:
            for seq in seqs:
                print(seq)
                fname = os.path.join(app.config["UPLOAD_DIR"], seq.filename) 
                print(fname)
                seq.save(fname)
        else:
            print("ERROR: no sequence files found.")

        return render_template('submit.html', title='Submit', form=form)
    return render_template('submit.html', title='Submit', form=form)
