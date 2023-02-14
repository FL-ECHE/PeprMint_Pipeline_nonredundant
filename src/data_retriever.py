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


""" Data download (previously on Notebook #1)

__author__ = ["Thibault Tubiana", "Phillippe Samer"]
__organization__ = "Computational Biology Unit, Universitetet i Bergen"
__copyright__ = "Copyright (c) 2022 Reuter Group"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Phillippe Samer"
__email__ = "samer@uib.no"
__status__ = "Prototype"
"""

import pandas as pd
from collections import defaultdict
import os

import shutil
import urllib
import tarfile
import gzip

from settings import Settings

'''
TO DO: not including IPython and related packages, nor settings e.g. activating
progress bars. 
Consider making an adequate interface to use everything in notebooks
'''

class DataRetriever:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings
        self.UPDATE = False

        # try to use pandarallel (only on GNU/Linux and OS X?)
        try:
            from pandarallel import pandarallel
            pandarallel.initialize(nb_workers=8, progress_bar=True)
            self.PARALLEL = True
        except:
            self.PARALLEL = False

    def retrieve_cath_domains(self):
        pass

    def retrieve_uniprot_to_pdb_correspondence(self):
        pass

    def retrieve_prosite(self):
        pass

    def retrieve_cath_pdb_files(self):
        pass

    def fetch_data_from_all(self):
        pass
