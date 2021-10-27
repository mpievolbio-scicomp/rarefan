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
        self.test_data_dir = os.path.abspath(os.path.join(project_dir, 'test/data/neisseria_completed'))
        self.session = {
                'tmpdir'    :   '',
                'outdir'    :   '',
                'email'     :   'no.name@no.host.xyz',
                'reference_strain' : "Nmen_2594",
                'e_value_cutoff': 1e-30,
                'nmer_length': 21,
                'min_nmer_occurence': 51,
                'query_rayt': '',
        }


    def tearDown(self):
        pass

    def test_get_status_code(self):

        status_code = routes.get_status_code(self.test_data_dir)
        self.assertEqual(111, status_code )  # add assertion here

    def test_get_email_command_incomplete(self):

        session = self.session
        session['tmpdir'] = os.path.abspath(os.path.join(project_dir, 'test/data/neisseria'))
        session['outdir'] = os.path.join(session['tmpdir'], 'out')
        session['reference_strain'] = "Nmen_2594"
        session['query_rayt'] = "yafM_Ecoli"

        cmd = routes.get_email_command(session)

        print(cmd)

        expected_cmd = """printf "Subject: Your RAREFAN run neisseria is complete.

Hallo,
your job on rarefan.evolbio.mpg.de with ID neisseria is complete.
    Job Summary
    ===========
        RAYTs
        -----
        Exit status: ERROR.

        We discovered 0 RAYTs using tblastn with
        yafM_Ecoli at an e-value threshold of 1e-30.

        NMERs
        -----
        Exit status: ERROR.

        There are 0 21 bp long sequences that
        occur more frequently than 51 times.

        REPINs
        ------
        Exit status: ERROR.

        We detected 0 REPINs.

    You can browse and download the results at this link:
    http://rarefan.evolbio.mpg.de/results?run_id=neisseria.

    In case of problems, please reply to this email and leave the email subject as is.

    Thank you for using RAREFAN.

    http://rarefan.evolbio.mpg.de
    " | msmtp rarefan@evolbio.mpg.de no.name@no.host.xyz >> /home/grotec/Repositories/RepinPop/test/data/neisseria/out/rarefan.log"""

        self.assertEqual(expected_cmd, cmd)


    def test_get_email_command_complete(self):

        session = self.session
        session['tmpdir'] = os.path.abspath(os.path.join(project_dir, 'test/data/neisseria_completed'))
        session['outdir'] = os.path.join(session['tmpdir'], 'out')
        session['reference_strain'] = "Nmen_2594"
        session['query_rayt'] = "yafM_Ecoli"

        cmd = routes.get_email_command(session)

        print(cmd)

        expected_cmd = """printf "Subject: Your RAREFAN run neisseria_completed is complete.

Hallo,
your job on rarefan.evolbio.mpg.de with ID neisseria_completed is complete.
    Job Summary
    ===========
        RAYTs
        -----
        Exit status: OK.

        We discovered 60 RAYTs using tblastn with
        yafM_Ecoli at an e-value threshold of 1e-30.

        NMERs
        -----
        Exit status: OK.

        There are 184 21 bp long sequences that
        occur more frequently than 51 times.

        REPINs
        ------
        Exit status: OK.

        We detected 48 REPINs.

    You can browse and download the results at this link:
    http://rarefan.evolbio.mpg.de/results?run_id=neisseria_completed.

    In case of problems, please reply to this email and leave the email subject as is.

    Thank you for using RAREFAN.

    http://rarefan.evolbio.mpg.de
    " | msmtp rarefan@evolbio.mpg.de no.name@no.host.xyz >> /home/grotec/Repositories/RepinPop/test/data/neisseria_completed/out/rarefan.log"""

        self.assertEqual(expected_cmd, cmd)


if __name__ == '__main__':
    unittest.main()
