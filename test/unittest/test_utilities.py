import unittest

# Import modules and functions to be tested.
import sys
import os

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
print(project_dir)

sys.path.insert(0, project_dir)
import app.utilities.check_errors as test_object
import logging
logging.getLogger().setLevel(logging.DEBUG)

class UtilitiesTest(unittest.TestCase):
    def setUp(self):
        self.test_data_dir = os.path.normpath("../data/neisseria_completed/out")

    def tearDown(self):
        pass

    def test_rayt_repin_counts(self):
        """ Test that number of rayts, reps, and repins are correctly counted from output."""

        counts = test_object.rayt_and_repin_counts(self.test_data_dir, "Nmen_2594")

        # Check return type
        self.assertIsInstance(counts, dict)
        self.assertIn("number_of_rayts", counts.keys())
        self.assertIn("number_of_overreps", counts.keys())
        self.assertIn("number_of_repins", counts.keys())

        logging.info(counts)

if __name__ == '__main__':
    unittest.main()
