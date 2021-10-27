# models.py
""" :module models: This module hosts all models of the rarefan web app. """

import os
import sys

from app import db

class Job(db.Document):
    """ :class Job: Represents all parameters, statuses, and results for a rarefan job. """
    run_id = db.StringField()
    setup = db.DynamicField() # Session cookie
    stages = db.DynamicField()
    results = db.DynamicField()


