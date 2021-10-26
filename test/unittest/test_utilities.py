import unittest

# Import modules and functions to be tested.
import sys
import os

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
sys.path.insert(0, project_dir)

from app.utilities import checkers
import logging
logging.getLogger().setLevel(logging.DEBUG)


class UtilitiesTest(unittest.TestCase):
    """ Testing the app utilities."""

    def setUp(self):
        """ Setup a test."""
        self.neisseria = os.path.join(project_dir, "test/data/neisseria_completed/out")

    def tearDown(self):
        """ Teardown a test."""
        pass

    def test_rayt_repin_counts_neisseria(self):
        """ Test that counting util for the reduced neisseria dataset."""

        results = checkers.parse_results(self.neisseria, "Nmen_2594")

        # Check return type
        self.assertIsInstance(results, dict)
        self.assertIn("counts", results.keys())
        self.assertIn("status", results.keys())

        counts = results['counts']
        status = results['status']

        self.assertIn("rayts", counts.keys())
        self.assertIn("overreps", counts.keys())
        self.assertIn("repins", counts.keys())

        self.assertIn("rayts", status.keys())
        self.assertIn("overreps", status.keys())
        self.assertIn("repins", status.keys())

        self.assertEqual(counts['rayts'], 60)
        self.assertEqual(counts['overreps'], 184)
        self.assertEqual(counts['repins'][0], 48)
        self.assertEqual(counts['repins'][1], 0)
        self.assertEqual(counts['repins'][2], 0)
        self.assertEqual(counts['repins'][3], 0)
        self.assertEqual(counts['repins'][4], 0)

        self.assertEqual(status['rayts'], 0)
        self.assertEqual(status['overreps'], 0)
        self.assertEqual(status['repins'], 0)

if __name__ == '__main__':
    unittest.main()
