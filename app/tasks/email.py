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
from app import mail

from flask_mail import Message

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


def email_task(dbjob):

    session = dbjob.setup
    results = checkers.parse_results(dbjob.setup['outdir'],
                                     dbjob.setup['reference_strain'])
    subject, body  = get_email(session, results)
    recipients = [session['email']]

    status = any([st not in ['finished', 'completed'] for st in [dbjob.stages[stage]['status'] for stage in ['rarefan', 'tree', 'zip']]])

    if status is False:
        recipients.append(app.config["MAIL_USERNAME"])

    recipients = [rec for rec in recipients if rec != ""]
    logging.debug("Recipients: %s", recipients)
    # Simply return if no recipients configured.
    if len(recipients) == 0:
        return

    sender = app.config['DEFAULT_MAIL_SENDER']

    message = Message(subject=subject,
                      body=body,
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


def get_email(session, results):
    # Aggregate the run path.
    run_id_path = session['tmpdir']
    run_id = os.path.basename(run_id_path)

    counts = results['counts']
    status = results['status']

    status_msg = {0: "passed", 1: "failed"}

    email_subject = "Your RAREFAN run {0:s} is complete.".format(run_id)
    email_body = f"""Hallo,
your job on rarefan.evolbio.mpg.de with ID {run_id} is complete.

Job Summary
===========

    RAYTs
    -----
    Data sanity check: {status_msg[status['rayts']]}.

    We discovered {counts['rayts']} RAYTs using tblastn with
    {session["query_rayt"]} at an e-value threshold of {session["e_value_cutoff"]}.

    Nmers
    -----
    Data sanity check: {status_msg[status['nmers']]}.

    There are {counts['nmers']} {session['nmer_length']}bp long sequences that
    occur more frequently than {session['min_nmer_occurence']} times.

    REPINs
    ------
    Data sanity check: {status_msg[status['repins']]}.

    We detected {sum(counts['repins'].values())} REPINs.

You can browse and download the results at http://rarefan.evolbio.mpg.de/results?run_id={run_id}.

In case of problems, please reply to this email and leave the email subject as is.

Thank you for using RAREFAN.

http://rarefan.evolbio.mpg.de
"""

    return email_subject, email_body


if __name__ == "__main__":
    print(email_test())
