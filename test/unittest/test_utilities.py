import unittest

# Import modules and functions to be tested.
import sys
import os

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
sys.path.insert(0, project_dir)

import app.utilities.check_errors as test_object
import logging
logging.getLogger().setLevel(logging.DEBUG)


class UtilitiesTest(unittest.TestCase):
    """ Testing the app utilities."""

    def setUp(self):
        """ Setup a test."""
        self.neisseria = os.path.normpath("../data/neisseria_completed/out")

    def tearDown(self):
        """ Teardown a test."""
        pass

    def test_rayt_repin_counts(self):
        """ Test that number of rayts, reps, and repins are correctly counted from output."""

        counts = test_object.rayt_repin_counts(self.neisseria, "Nmen_2594")

        # Check return type
        self.assertIsInstance(counts, dict)
        self.assertIn("number_of_rayts", counts.keys())
        self.assertIn("number_of_overreps", counts.keys())
        self.assertIn("number_of_repins", counts.keys())

        self.assertEqual(counts['number_of_rayts'], 60)
        self.assertEqual(counts['number_of_overreps'], 184)
        self.assertEqual(counts['number_of_repins'][0], 0)
        self.assertEqual(counts['number_of_repins'][1], 0)
        self.assertEqual(counts['number_of_repins'][2], 0)
        self.assertEqual(counts['number_of_repins'][3], 0)
        self.assertEqual(counts['number_of_repins'][4], 0)

        logging.info(counts)


if __name__ == '__main__':
    unittest.main()
