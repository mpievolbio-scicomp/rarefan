import subprocess
import os
import sys
import shlex
import shutil
import subprocess

from app import app
from app.utilities import checkers
from app import mail

from flask_mail import Message
from flask import g

import logging
logging.basicConfig(level=logging.DEBUG)

app.app_context().push()

def email_test():

    # app.config["TESTING"] = False
    print(app.config)
    logging.debug("TESTING is %s", app.config['TESTING'])

    recipients = ['grotec@evolbio.mpg.de']

    message = Message("Test from flask",
                      sender="rarefan@evolbio.mpg.de",
                      recipients=recipients,
                      body="It works!"
                      )

    success = mail.send(message)

    logging.debug("Mail was sent: %s", str(success))
    logging.debug("Mail message was: %s", str(message))

    # assert g.outbox[0].subject == 'Test from flask'
    # app.config["TESTING"] = True
    return success, message


def email_task(dbjob):
    pass

#     session = dbjob.setup
#     results = checkers.parse_results(dbjob.setup['outdir'],
#                                      dbjob.setup['reference_strain'])
#     subject, body  = get_email(session, results)
#     recipients = []

#     user_email = session['email']
#     if user_email is not None or user_email != "None":
#         recipients.append(user_email)

#     status = any([st not in ['finished', 'completed'] for st in [dbjob.stages[stage]['status'] for stage in ['rarefan', 'tree', 'zip']]])

#     if status is False:
#         recipients.append(app.config["MAIL_USERNAME"])

#     message = Message(subject,
#                       sender=app.config['MAIL_USERNAME'],
#                       recipients=recipients,
#                       body=body)

#     mail.send(msg)


# def get_email(session, results):
#     # Aggregate the run path.
#     run_id_path = session['tmpdir']
#     run_id = os.path.basename(run_id_path)

#     # List of recipients. Always send to rarefan.
#     if session['email'] is None or session['email'] == "":
#         logging.debug("No email set.")

#     counts = results['counts']
#     status = results['status']

#     status_msg = {0: "OK", 1: "ERROR"}

#     # Append support to recipients if there was an error.
#     if sum(status.values()) > 0:
#         recipients.append('rarefan@evolbio.mpg.de')

#     if len(recipients) == 0:
#         return None

#     email_subject = "Your RAREFAN run {0:s} is complete.".format(run_id)
#     email_body = f"""Hallo,
# your job on rarefan.evolbio.mpg.de with ID {run_id} is complete.

# Job Summary
# ===========

#     RAYTs
#     -----
#     Exit status: {status_msg[status['rayts']]}.

#     We discovered {counts['rayts']} RAYTs using tblastn with
#     {session["query_rayt"]} at an e-value threshold of {session["e_value_cutoff"]}.

#     NMERs
#     -----
#     Exit status: {status_msg[status['overreps']]}.

#     There are {counts['overreps']} {session['nmer_length']} bp long sequences that
#     occur more frequently than {session['min_nmer_occurence']} times.

#     REPINs
#     ------
#     Exit status: {status_msg[status['repins']]}.

#     We detected {sum(counts['repins'].values())} REPINs.

# You can browse and download the results at this link:
# http://rarefan.evolbio.mpg.de/results?run_id={run_id}.

# In case of problems, please reply to this email and leave the email subject as is.

# Thank you for using RAREFAN.

# http://rarefan.evolbio.mpg.de
# """

#     return email_body, email_subject


if __name__ == "__main__":
    print(email_test())
