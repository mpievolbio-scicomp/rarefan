import unittest

# Import modules and functions to be tested.
import sys
import os

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
print(project_dir)

sys.path.insert(0, project_dir)
import app.routes as routes
import logging

routes.logger.setLevel(logging.DEBUG)

class RoutesTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_get_status_code(self):

        status_code = routes.get_status_code("../data/neisseria_completed")
        self.assertEqual(111, status_code )  # add assertion here

    def test_get_email_command_incomplete(self):

        session = {
                'tmpdir'    :   '../data/neisseria',
                'email'     :   'no.name@no.host.xyz',
        }
        cmd = routes.get_email_command(session)

        expected_cmd = """printf "Subject: Your RAREFAN run neisseria has finished.

Hallo,
your job on rarefan.evolbio.mpg.de with ID neisseria has finished.
You can browse and download the results at this link:
http://rarefan.evolbio.mpg.de/results?run_id=neisseria.

In case of problems, please reply to this email and leave the email subject as is.

Thank you for using RAREFAN.

http://rarefan.evolbio.mpg.de
" | msmtp no.name@no.host.xyz >> ../data/neisseria/out/rarefan.log"""

        self.assertEqual(expected_cmd, cmd)


    def test_get_email_command_complete(self):

        session = {
                'tmpdir'    :   '../data/neisseria_completed',
                'email'     :   'no.name@no.host.xyz',
        }
        cmd = routes.get_email_command(session)

        expected_cmd = """printf "Subject: Your RAREFAN run neisseria_completed has finished.

Hallo,
your job on rarefan.evolbio.mpg.de with ID neisseria_completed has finished.
You can browse and download the results at this link:
http://rarefan.evolbio.mpg.de/results?run_id=neisseria_completed.

In case of problems, please reply to this email and leave the email subject as is.

Thank you for using RAREFAN.

http://rarefan.evolbio.mpg.de
" | msmtp no.name@no.host.xyz >> ../data/neisseria_completed/out/rarefan.log"""

        self.assertEqual(expected_cmd, cmd)


if __name__ == '__main__':
    unittest.main()
