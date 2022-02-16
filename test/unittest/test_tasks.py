
import unittest

# Import modules and functions to be tested.
import sys
import os
import shutil

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
sys.path.insert(0, project_dir)

import logging
logging.getLogger().setLevel(logging.DEBUG)

from app.tasks.tree import tree_task
from app.tasks.rayt_phylo import alignment_task, phylogeny_task

class TasksTest(unittest.TestCase):
    """ Testing the app utilities."""

    def setUp(self):
        """ Setup a test."""
        self.neisseria = os.path.join(project_dir, "test/data/neisseria_small")
        self.run_dir = os.path.join('/tmp', 'neisseria_small')
        self.out_dir = os.path.join(self.run_dir, 'out')
        shutil.copytree(self.neisseria, self.run_dir)
        os.mkdir(os.path.join(self.out_dir))

    def tearDown(self):
        """ Teardown a test."""
        shutil.rmtree('/tmp/neisseria_small')

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

    def test_rayt_align_task(self):
        """ Test the task for computing the rayt alignment ."""

        ret, log = alignment_task(self.run_dir)

        self.assertIn('raytAln.phy', os.listdir(self.out_dir))

    def test_rayt_phylogeny_task(self):
        """ Test the task for computing the rayt phylogeny ."""

        ret, log = phylogeny_task(self.run_dir)

        self.assertIn('raytAln.phy_phyml_tree.txt', os.listdir(self.out_dir))
        self.assertIn('raytAln.phy_phyml_stats.txt', os.listdir(self.out_dir))

if __name__ == '__main__':
    unittest.main()
