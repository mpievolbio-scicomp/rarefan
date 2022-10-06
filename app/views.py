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
        description="Select the reference strain from your uploaded genome sequences.",
            )

    query_rayt = SelectField(
            "Query rayt",
            choices = ['yafM_Ecoli',
                       'yafM_SBW25',
                       ], 
            description="Select the RAYT protein sequence to use for identfying RAYT sequences in the provided genomes. 'yafM_Ecoli' and 'yafM_SBW25' are provided by default. You can supply your own RAYTs in the 'Upload' step.",
            validators=[validators.DataRequired()]
            )

    treefile = SelectField(
            "Tree file",
            choices=[],
        description="Optional: If your uploaded data contains a '.nwk'  phylogenetic tree file, you may select it here and it will be used in the postprocessing step. Inferred RAYT and REP/REPIN population frequencies will be plotted against the phylogeny. Select 'None' to calculate the tree from the submitted genomes.",
            )

    min_nmer_occurrence = IntegerField("Min. seed sequence occurrence", 
                          default=55,
        description="Set the seed sequence occurrence cutoff. Only seeds of 'Seed length' basepairs (see below) that occur more frequently than this number will be considered in the REP/REPIN analysis.",
                          validators=[validators.DataRequired(message="Please enter the minimal seed occurrence as an integer!")]
                          )

    nmer_length = IntegerField("Seed length", 
                             default=21, 
        description="Set the Seed length (in basepairs). Only sequences of this length that occur more frequently than 'min. Seed occurrence (see above) will be considered in the REP/REPIN analysis.",
                             validators=[validators.DataRequired(message="Please enter the nmer length as an integer!")]
                             )
    distance_group_seeds= IntegerField("Distance group seeds", 
                             default=15, 
        description="Set the group seeds distance. Determines whether REPINs are grouped accurately. ",
                             validators=[validators.DataRequired(message="Please enter an integer number > 0.")]
                             )

    e_value_cutoff = FloatField("e value cutoff",
                                default=1.0e-30,
        description="Set the e-value cutoff for 'tblastn'  alignment of query RAYT protein sequence (selected above in 'Query RAYT') to input genomes. The e-value should be given in scientific 'e' notation.",
                                validators=[validators.DataRequired(message="Please enter the e value cutoff in scientific notation (e.g. 1e-30)")]
                                )

    analyse_repins = BooleanField("Analyse REPINs",
                                  default=True,
        description="Toggle REPIN analysis. If unchecked, only REPs will be analysed.",
                                  )

    email = StringField("Optional: Your email address.",
                        validators=[OptionalEmail()],
        description="Provide your email address to receive a notification when the job is finished.",
                        )

    submit = SubmitField("Go!",
        description="Click here to submit the job.",
    )


class RunForm(FlaskForm):
    confirm = SubmitField("Confirm file upload")
    proceed = SubmitField("Proceed to run configuration")


class ReturnToResultsForm(FlaskForm):
    back_button = SubmitField("Back to Results")


class AnalysisForm(FlaskForm):
    run_id = StringField("Enter run ID")
    go = SubmitField("Ok")

