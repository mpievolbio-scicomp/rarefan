# models.py
""" :module models: This module hosts all models of the rarefan web app. """

import os
import sys
import rq
from rq.exceptions import NoSuchJobError
from app import app, db
logger = app.logger


class Job(db.Document):
    """ :class Job: Represents all parameters, statuses, and results for a rarefan job. """
    session = db.StringField()
    parent_run = db.StringField()
    run_id = db.StringField()
    setup = db.DictField() # Session cookie
    stages = db.DictField()
    overall_status = db.StringField(default='setup')
    notification_is_sent = db.BooleanField(default=False)

    def get_status_from_redis(self, stage):
        """ Infer the status of a job stage. """
        redis_job_id = self.stages.get(stage)['redis_job_id']
        try:
            redis_job = rq.job.Job.fetch(redis_job_id, connection=app.redis)
            redis_job_status = redis_job.get_status()
        except NoSuchJobError:
            redis_job_status = "complete"
        except:
            redis_job_status = "failed"

        return redis_job_status

    def set_status(self, stage, value=None):

        if stage not in self.stages.keys():
            raise AttributeError("{} is not a rarefan job stage. Valid stages are {}".format(stage, str(list(self.stages.keys()))))

        if value is None:
            value = self.get_status_from_redis(stage)
        if value != self.get_status_from_redis(stage):
                raise RuntimeError("Supplied status '{}' is not consistent with status reported from redis queue ({}).".format(value, self.get_status_from_redis(stage)))

        if stage == 'rarefan': self.update(set__stages__rarefan__status=value)
        if stage == 'rayt_alignment': self.update(set__stages__rayt_alignment__status=value)
        if stage == 'rayt_phylogeny': self.update(set__stages__rayt_phylogeny__status=value)
        if stage == 'tree': self.update(set__stages__tree__status=value)
        if stage == 'zip': self.update(set__stages__zip__status=value)


    def set_overall(self):
        # Update stati
        stage_names =  self.stages.keys()
        for stage in stage_names:
            self.set_status(stage)
            logger.debug("%s status= %s", stage, self.stages[stage]['status'] )

        overall = 'setup'
        if all([ stage['status'] == "none" for stage in self.stages.values()]):
            overall = "setup"
        elif self.stages['rarefan']['status'] == "queued":
            overall = "queued"
        elif self.stages['rarefan']['status'] in ["started", "running"]:
            overall = 'running'
        elif self.stages['rarefan']['status'] == "failed":
            overall = 'failed'

        elif self.stages['rarefan']['status'] in ["complete", "finished"]:
            if any([self.stages[stage]['status'] in ["started", "running"] for stage in ['rayt_alignment', 'rayt_phylogeny', 'tree', 'zip']]):
                overall = "postprocessing"

            else:
                overall = "complete"

        self.overall_status = overall
        self.update(set__overall_status=overall)














