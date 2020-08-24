# views.py

from flask import Flask, render_template, flash, request

from flask_wtf import FlaskForm
from wtforms import StringField,\
                    PasswordField,\
                    BooleanField,\
                    SubmitField,\
                    MultipleFileField,\
                    SelectField

from wtforms.validators import DataRequired
from app import app

from wtforms import Form, TextField, TextAreaField, validators, StringField, SubmitField

class SequenceUploadForm(FlaskForm):
    sequences = MultipleFileField('Sequences', validators=[DataRequired(),])
    upload = SubmitField('Upload')

class ReferenceStrainForm(FlaskForm):
    reference_strain = SelectField('Reference sequence',
            choices=[],
            validate_choice=True)

    rayt = StringField('Rayt')
class SubmitForm(FlaskForm):
    submit = SubmitField('Submit')


