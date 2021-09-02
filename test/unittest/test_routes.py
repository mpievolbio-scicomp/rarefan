import unittest

# Import modules and functions to be tested.
import sys
import os

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
print(project_dir)

sys.path.insert(0, project_dir)
import app.routes as routes

class RoutesTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_get_status_code(self):

        status_code = routes.get_status_code("../data/neisseria")
        self.assertEqual(111, status_code )  # add assertion here

    def test_get_email_command(self):

        session = {
                'tmpdir'    :   '../data/neisseria',
                'email'     :   'no.name@no.host.xyz',
        }
        cmd = routes.get_email_command(session)

        expected_cmd = """printf "Subject: Your RAREFAN run neisseria_completed has finished.

Hallo,
your job on rarefan.evolbio.mpg.de with ID neisseria_completed has finished.
You can browse and download the results at this link:
http://rarefan.evolbio.mpg.de/results?run_id=neisseria_completed.

Thank you for using RAREFAN. We hope to see you soon again.

Kind regards,

RAREFAN.

http://rarefan.evolbio.mpg.de
" | msmtp no.name@no.host.xyz >> ../data/neisseria_completed/out/rarefan.log"""

        self.assertEqual(cmd, expected_cmd)


if __name__ == '__main__':
    unittest.main()
