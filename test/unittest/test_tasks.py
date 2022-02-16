
import unittest

# Import modules and functions to be tested.
import sys
import os
import shutil
import difflib

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
sys.path.insert(0, project_dir)

import logging
logging.getLogger().setLevel(logging.DEBUG)

from app.tasks.tree import tree_task
from app.tasks.rayt_phylo import alignment_task, phylogeny_task, run_alignment, run_phyml

class TasksTest(unittest.TestCase):
    """ Testing the app utilities."""

    def setUp(self):
        """ Setup a test."""
        self.neisseria = os.path.join(project_dir, "test/data/neisseria_small")
        self.run_dir = os.path.join('/tmp', 'neisseria_small')
        self.reference_out_dir = os.path.join(self.neisseria, 'out')
        self.out_dir = os.path.join(self.run_dir, 'out')

        # Copy input data.
        shutil.copytree(self.neisseria, self.run_dir)

        # Remove some files that will be generated during testing.
        for f in ['tmptree.nwk',
                  'raytAln.phy',
                  'raytAln.phy_phyml_stats.txt',
                  'raytAln.phy_phyml_tree.txt',
                  ] :
            if os.path.exists(f):
                os.remove(f)


        self.__files_to_remove = []
        self.__dirs_to_remove = [self.run_dir]

    def tearDown(self):
        """ Teardown a test."""

        for f in self.__files_to_remove:
            if os.path.isfile(f):
                os.remove(f)

        for d in self.__dirs_to_remove:
            if os.path.isdir(d):
                shutil.rmtree(d)


    def test_tree_task_default(self):
        """ Test the tree generation task with default treefile name.
        """

        ret, log = tree_task(self.run_dir, None)

        self.assertIn('tmptree.nwk', os.listdir(self.out_dir))

        with open(os.path.join(self.out_dir, 'tmptree.nwk')) as fh:
            tree = "".join(fh.readlines())

        self.assertEqual(tree, "(Nmen_14-563:0.008237,(Nmen_331401:0.001783,Nmen_2594:0.001417):0.007763,(Nmen_38277:0.008887,Ngon_NJ1711654:0.027612):0.000413);\n")

    def test_tree_task(self):
        """ Test the tree generation task
        """

        ret, log = tree_task(self.run_dir, 'tmptree.nwk')

        self.assertIn('tmptree.nwk', os.listdir(self.out_dir))

        with open(os.path.join(self.out_dir, 'tmptree.nwk')) as fh:
            tree = "".join(fh.readlines())
        self.assertEqual(tree, "(Nmen_14-563:0.008237,(Nmen_331401:0.001783,Nmen_2594:0.001417):0.007763,(Nmen_38277:0.008887,Ngon_NJ1711654:0.027612):0.000413);\n")

    def test_rayt_phylo_run_alignment(self):
        """ Test the task for computing the rayt alignment ."""

        ret, log = run_alignment(self.run_dir)

        # Check expected output file is present.
        expected_out_fname = 'raytAln.phy'
        self.assertIn(expected_out_fname, os.listdir(self.out_dir))

        # Check output data is equal to reference data.
        test_out_fname = os.path.join(self.out_dir, 'raytAln.phy')
        reference_out_fname = os.path.join(self.reference_out_dir, 'raytAln.phy') 
        
        shutil.copy(test_out_fname, '/tmp/raytAln.phy')

        with open(test_out_fname, 'r') as ifh:
            test_alignment = ifh.read()
        with open(reference_out_fname, 'r') as ifh:
            reference_alignment = ifh.read()

        # Get the diff between reference data and test run data.
        diff = [line for line in difflib.unified_diff(test_alignment, reference_alignment)]

        self.assertEqual(len(diff), 0)


    def test_rayt_phylo_run_phyml(self):
        """ Test the task for computing the rayt phylogeny ."""

        ret, log = run_phyml(self.run_dir, seed=1645040163)

        # Check tree and stats files are present.
        expected_tree_fname = 'raytAln.phy_phyml_tree.txt'
        expected_stats_fname = 'raytAln.phy_phyml_stats.txt'
        self.assertIn(expected_tree_fname, os.listdir(self.out_dir))
        self.assertIn(expected_stats_fname, os.listdir(self.out_dir))

        # Compare generated data and reference data.
        test_tree_fname = os.path.join(self.out_dir, expected_tree_fname)
        reference_tree_fname = os.path.join(self.reference_out_dir, expected_tree_fname)

        with open(test_tree_fname, 'r') as ifh:
            test_tree_data = ifh.read()
        with open(reference_tree_fname, 'r') as ifh:
            reference_tree_data = ifh.read()

        # Get diffs.
        tree_diff = [line for line in difflib.unified_diff(test_tree_data, reference_tree_data)]

        self.assertSequenceEqual(tree_diff, [])

if __name__ == '__main__':
    unittest.main()
