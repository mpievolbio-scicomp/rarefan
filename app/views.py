# views.py

from flask import Flask, render_template, flash, request, session

from flask_wtf import FlaskForm
from wtforms import BooleanField
from wtforms import SelectField
from wtforms import MultipleFileField
from wtforms import StringField
from wtforms import SubmitField
from wtforms import IntegerField
from wtforms import FloatField
from wtforms import validators


from wtforms.validators import DataRequired
from flask_wtf.file import FileRequired

import tempfile
from app import app


class SubmitForm(FlaskForm):
    sequences = MultipleFileField('Sequences',
            )
    upload = SubmitField("Upload!")

    
    reference_strain = SelectField(
            'Reference sequence',
            choices=[],
            )

    query_rayt = SelectField(
            "Query rayt",
            choices = ['yafM_Ecoli',
                       'yafM_SBW25',
                       ], 
            validators=[DataRequired(),]
            )

    treefile = SelectField(
            "Tree file",
            choices=[],
            )

    min_nmer_occurence = IntegerField("Min. nmer occurence", 
                          default=55,
                          validators=[DataRequired(),]
                          )

    nmer_length = IntegerField("Nmer length", 
                             default=21, 
                             validators=[DataRequired(),]
                             )

    e_value_cutoff = FloatField("e value cutoff",
                                default=1.0e-30,
                                validators=[DataRequired(),]
                                )

    analyse_repins = BooleanField("Analyse REPINs",
                                  default=True,
                                  description="Leave unchecked to analyse REPs only."
                                  )

    email = StringField("Optional: provide your email address to receive a notification once your job is done.")

    go = SubmitField("Go!")

class RunForm(FlaskForm):
    go = SubmitField("GoGo!")


class AnalysisForm(FlaskForm):
    run_id = StringField("Enter run ID")
    go = SubmitField("Download zip")

