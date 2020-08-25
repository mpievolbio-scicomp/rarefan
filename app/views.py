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
from flask_wtf.file import FileRequired
from app import app

from wtforms import Form,\
                    TextField,\
                    TextAreaField,\
                    StringField,\
                    SubmitField,\
                    IntegerField,\
                    FloatField,\
                    validators

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

    go = SubmitField("Go!")

