#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# MIT License
# 
# Copyright (c) 2022 Reuter Group
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


""" Tests related to Preprocessing, especially running cath-superpose

__author__ = ["Phillippe Samer"]
__organization__ = "Computational Biology Unit, Universitetet i Bergen"
__copyright__ = "Copyright (c) 2022 Reuter Group"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Phillippe Samer"
__email__ = "samer@uib.no"
__status__ = "Prototype"
"""

import unittest
import pytest

from src.settings import Settings
from src.dataset_manager import DatasetManager
from src.ibs_tagging import IBSTagging

class TestPreprocessing(unittest.TestCase):

    def test_fetch_cath_superpose_binary(self):
        # TO DO: all
        assert True

    def test_run_cath_superpose(self):
        # TO DO: all
        assert True

    def test_orient_on_z_axis(self):
        # TO DO: all
        assert True

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)