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


""" Methods to download/use AlphaFold structures (previously on Notebook #3)

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
import numpy as np
import math
import seaborn as sns
import urllib
import glob
import os
from urllib.error import HTTPError

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import MDAnalysis as mda
import nglview as nv
from scipy.spatial import ConvexHull

from tqdm import tnrange, tqdm

from src.settings import Settings
from src.notebook_handle import NotebookHandle

class AlphaFoldUtils:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings
        self.dataset = None

        self._libs_setup()

    def _libs_setup(self):
        # Seaborn
        sns.set_style("darkgrid")

        # Pandas
        pd.options.mode.chained_assignment = (
            None  # remove warning when adding a new column; default='warn'
        )
        pd.set_option("display.max_columns", None)
        tqdm.pandas()   # activate tqdm progressbar for pandas

        # Numpy
        np.seterr(divide='ignore', invalid='ignore')

        # IPython
        if self.settings.USING_NOTEBOOK:
            self.settings.NOTEBOOK_HANDLE.alphafold_utils_options()

    def run(self, df: pd.DataFrame):
        # runs each step to fetch and prepare alphafold data (as of Notebook #3)
        self.dataset = df
        self.process_all_domains()
        self.align_fasta_files()
        self.save_PH_domains_from_alphafold()

    """
    ### All methods below just encapsulate the steps in Notebook #3
    """

    def process_all_domains(self):
        pass

    def _fetch_pdb_alfafold(self):
        pass

    def _get_prosite_boundaries_dict(self):
        pass

    def _get_json(self):
        pass

    def _get_domain_fragment_query(self):
        pass

    def align_fasta_files(self):
        pass

    def save_PH_domains_from_alphafold(self):
        pass
