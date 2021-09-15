# views.py

from flask_wtf import FlaskForm
from flask_wtf.file import FileRequired, FileAllowed
from wtforms import BooleanField
from wtforms import SelectField
from wtforms import MultipleFileField
from wtforms import StringField
from wtforms import SubmitField
from wtforms import IntegerField
from wtforms import FloatField
from wtforms import validators
from wtforms.validators import ValidationError

# Needed to validate fasta files.
# https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta


class OptionalEmail(validators.Email):
    def __call__(self, form, field):
        if field.data is None or field.data == "":
            return None

        super().__call__(form, field)


class SequenceFilesValidator():
    allowed_extensions = ['fasta', 'fastn', 'fa', 'fas', 'fna', 'fn', 'faa', 'nwk']
    def __call__(self, form, field):
        extensions = [f.filename.split(".")[-1] for f in field.data]
        if not all([ext in self.allowed_extensions for ext in extensions]):
            raise ValidationError("Only files with these filename extensions are allowed: {}".format(["."+ext for ext in self.allowed_extensions]))

        if len([True for ext in extensions if ext == 'nwk'] ) > 1:
            raise ValidationError("At max 1 tree file can be uploaded.")


class UploadForm(FlaskForm):
    sequences = MultipleFileField('File upload',
                                   validators=[validators.DataRequired(),
                                               validators.Length(min=1, message="Please select at least one sequence file (fasta format) and (optionally) one tree file."),
                                               SequenceFilesValidator(),
                                              ]
                                  )

    upload = SubmitField("Upload!")


class SubmitForm(FlaskForm):

    reference_strain = SelectField(
            'Reference sequence',
            choices=[],
            )

    query_rayt = SelectField(
            "Query rayt",
            choices = ['yafM_Ecoli',
                       'yafM_SBW25',
                       ], 
            validators=[validators.DataRequired()]
            )

    treefile = SelectField(
            "Tree file",
            choices=[],
            )

    min_nmer_occurence = IntegerField("Min. nmer occurence", 
                          default=55,
                          validators=[validators.DataRequired(message="Please enter the minimal nmer occurence as an integer!")]
                          )

    nmer_length = IntegerField("Nmer length", 
                             default=21, 
                             validators=[validators.DataRequired(message="Please enter the nmer length as an integer!")]
                             )

    e_value_cutoff = FloatField("e value cutoff",
                                default=1.0e-30,
                                validators=[validators.DataRequired(message="Please enter the e value cutoff in scientific notation (e.g. 1e-30)")]
                                )

    analyse_repins = BooleanField("Analyse REPINs",
                                  default=True,
                                  description="Leave unchecked to analyse REPs only."
                                  )

    email = StringField("Optional: Your email address.",
                        validators=[OptionalEmail()]
                        )

    submit = SubmitField("Go!")


class RunForm(FlaskForm):
    confirm = SubmitField("Confirm file upload")
    proceed = SubmitField("Proceed to run configuration")


class ReturnToResultsForm(FlaskForm):
    back_button = SubmitField("Back to Results")


class AnalysisForm(FlaskForm):
    run_id = StringField("Enter run ID")
    go = SubmitField("Download zip")

