# models.py
""" :module models: This module hosts all models of the rarefan web app. """

import os
import sys

from app import db

# class Result(db.EmbeddedDocument):
#     returncode = db.IntField()
#     log = db.StringField()

# class Stage(db.EmbeddedDocument):
#     redis_job_id = db.StringField()
#     result = db.EmbeddedDocumentField(Result)
#     status = db.StringField()

# class Stages(db.EmbeddedDocument):
#     rarefan = db.EmbeddedDocumentField(Stage)
#     tree = db.EmbeddedDocumentField(Stage)
#     zip = db.EmbeddedDocumentField(Stage)

# class Job(db.EmbeddedDocument):
#     """ :class Job: Represents all parameters, statuses, and results for a rarefan job. """
#     run_id = db.StringField()
#     setup = db.DynamicField() # Session cookie
#     stages = db.EmbeddedDocumentField(Stages)

class Job(db.Document):
    """ :class Job: Represents all parameters, statuses, and results for a rarefan job. """
    run_id = db.StringField()
    setup = db.DynamicField() # Session cookie
    stages = db.DictField()




