import glob
import os
import unittest

import blazar_clean
from blazar_properties import *


class TestClean(unittest.TestCase):
    def test_clean_missing_dir(self):
        self.assertFalse(blazar_clean.clean(data_folder="i-do-not-exist", parameter_folder="i-also-do-not-exist"))

    def test_clean_no_action(self):
        data_folder = "blazar_mcmc/tests/test_data"
        data_full = BASE_PATH + data_folder
        param_folder = "blazar_mcmc/tests/test_params"
        param_full = BASE_PATH + param_folder

        # create data file to test clearing
        if not os.path.exists(data_full):
            os.mkdir(data_full)
        with open(data_full + "/test.dat", 'w') as f:
            f.write("testing")
        # create file to test clearing
        if not os.path.exists(param_full):
            os.mkdir(param_full)
        with open(param_full + "/test.txt", 'w') as f:
            f.write("testing")

        num_files_data = len(glob.glob(data_full))
        num_files_params = len(glob.glob(data_full))

        self.assertTrue(blazar_clean.clean(data=False, data_folder=data_folder, parameter_files=False,
                                           parameter_folder=param_folder))

        self.assertTrue(os.path.exists(data_full + "/test.dat"))
        self.assertTrue(os.path.exists(param_full + "/test.txt"))
        self.assertEqual(num_files_data, len(glob.glob(data_full)))
        self.assertEqual(num_files_params, len(glob.glob(param_full)))

    def test_clean(self):
        data_folder = "blazar_mcmc/tests/test_data"
        data_full = BASE_PATH + data_folder
        param_folder = "blazar_mcmc/tests/test_params"
        param_full = BASE_PATH + param_folder

        # create data file to test clearing
        if not os.path.exists(data_full):
            os.mkdir(data_full)
        with open(data_full + "/test.dat", 'w') as f:
            f.write("testing")
        # create file to test clearing
        if not os.path.exists(param_full):
            os.mkdir(param_full)
        with open(param_full + "/test.txt", 'w') as f:
            f.write("testing")

        # it's not supposed to delete this
        if not os.path.exists(param_full + "/params.txt"):
            with open(param_full + "/params.txt", 'w') as f:
                f.write("testing")
        self.assertTrue(blazar_clean.clean(data=True, data_folder=data_folder, parameter_files=True,
                                           parameter_folder=param_folder))

        # folders not deleted
        self.assertTrue(os.path.exists(data_full))
        self.assertTrue(os.path.exists(param_full))
        # folders empty except for params.txt
        self.assertTrue(os.path.exists(param_full + "/params.txt"))
        self.assertEqual(len(glob.glob(data_full + "/*")), 0)
        self.assertEqual(len(glob.glob(param_full + "/*")), 1)


if __name__ == '__main__':
    unittest.main()
