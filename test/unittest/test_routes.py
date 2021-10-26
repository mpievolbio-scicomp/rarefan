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

        expected_cmd = """ """

        self.assertEqual(expected_cmd, cmd)


    def test_get_email_command_complete(self):

        session = self.session
        session['tmpdir'] = os.path.abspath(os.path.join(project_dir, 'test/data/neisseria_completed'))
        session['outdir'] = os.path.join(session['tmpdir'], 'out')
        session['reference_strain'] = "Nmen_2594"
        session['query_rayt'] = "yafM_Ecoli"

        cmd = routes.get_email_command(session)

        print(cmd)

        expected_cmd = """ """

        self.assertEqual(expected_cmd, cmd)


if __name__ == '__main__':
    unittest.main()
