# views.py

from flask import Flask, render_template, flash, request

from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField
from wtforms.validators import DataRequired

from app import app

from wtforms import Form, TextField, TextAreaField, validators, StringField, SubmitField
class LoginForm(FlaskForm):
    sequences = StringField('Sequences', validators=[DataRequired()])
    reference_sequence = StringField('Reference sequence', validators=[DataRequired()])
    rayt = StringField('Rayt')
    submit = SubmitField('Submit')
