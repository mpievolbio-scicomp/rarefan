import subprocess
import copy
import os
import sys
import shlex
import shutil
import re
import subprocess

from app import app
from app.utilities import checkers
from app import mail, views

from flask_mail import Message
from flask import render_template
from app.views import AnalysisForm
from app.models import Job as DBJob

import jinja2
import logging
logging.basicConfig(level=logging.DEBUG)

app.app_context().push()

def email_test():

    recipients = ['grotec@evolbio.mpg.de']

    message = Message("Test from flask",
                      sender="rarefan@evolbio.mpg.de",
                      recipients=recipients,
                      body="It works!"
                      )

    success = mail.send(message)

    logging.debug("Mail was sent: %s", str(success))
    logging.debug("Mail message was: %s", str(message))

    return success, message


def email_task(run_id):

    dbjob = DBJob.objects.get(run_id=run_id)
    # Get session.
    session = dbjob.setup

    # List of recipients
    recipients = [session['email']]

    # Overall status.
    status = dbjob.overall_status

    if status == "failed":
        recipients.append(app.config["MAIL_USERNAME"])

    # Only non-empty strings in recipients list.
    recipients = [rec for rec in recipients if rec != ""]
    logging.debug("Recipients: %s", recipients)

    # Simply return if no recipients configured.
    if len(recipients) == 0:
        return

    # Get the email subject.
    email_subject = "Your RAREFAN run {0:s} is complete.".format(dbjob.run_id)

    # Assemble email body.
    template_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', app.template_folder))
    templateLoader = jinja2.FileSystemLoader(searchpath=template_dir)
    templateEnv = jinja2.Environment(loader=templateLoader, autoescape=True)
    template_file = "email.md"
    template = templateEnv.get_template(template_file)

    email_body = template.render(job=dbjob)  # this is where to put args to the template renderer


    sender = app.config['DEFAULT_MAIL_SENDER']

    message = Message(subject=email_subject,
                      body=email_body,
                      sender=sender,
                      recipients=recipients,
                      )
    logging.debug("Message: %s", message)

    is_sent = copy.deepcopy(dbjob.notification_is_sent)
    while not is_sent:
        mail_status = mail.send(message)
        logging.debug("Mail status %s", str(mail_status))
        logging.debug("Mail status type: %s", str(type(mail_status)))

        if mail_status is None:
            logging.debug("Updating Mail status.")
            is_sent = True
    dbjob.update(set__notification_is_sent=copy.deepcopy(is_sent))
    dbjob.save()

    logging.debug("Mail sent? %s", str(is_sent))
