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


""" Tagging the IBS in the dataset (previously on notebooks under "tools")

Methods to tag the interfacial binding sites (IBS) of proteins in the dataset, 
with or without using the alphafold models.

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

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.spatial import ConvexHull
import MDAnalysis as mda

from tqdm import tnrange, tqdm

import urllib
import glob
from urllib.error import HTTPError

from src.settings import Settings
from src.notebook_handle import NotebookHandle

class IBSTagging:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings
        self.dataset = None

        self._libs_setup()

        # TO DO:
        # read settings: type of data, clustering level, uniref level, z axis level, COMPARISONMODE

    def _libs_setup(self):
        # IPython
        if self.settings.USING_NOTEBOOK:
            self.settings.NOTEBOOK_HANDLE.ibs_options()
        
        # Pandas
        pd.options.mode.chained_assignment = (
            None  # remove warning when adding a new column; default='warn'
        )
        pd.set_option("display.max_columns", None)
        tqdm.pandas()   # activate tqdm progressbar for pandas

        # Numpy
        np.seterr(divide='ignore', invalid='ignore')

        # Seaborn
        sns.set_style("darkgrid")

    def run(self, df: pd.DataFrame):
        self.dataset = df
        self._PH_domain_preprocessing()
        self._tag_PH()
        self._tag_C2()
        self._tag_START()
        self._tag_C1()
        self._tag_C2DIS()
        self._tag_PX()
        self._tag_PLD()
        self._tag_ANNEXIN()
        self._tag_PLA()
        self._merge_datasets()

    """
    ### All methods below just encapsulate the steps in the original notebook
    """

    # Set of Uniref structures to take for comparison test (between AF and Cath)
    def _get_uniprot_in_common(self, domain):
        pass

    # Prepare Exclusion of PTB and other domains
    def _PH_domain_preprocessing(self):
        pass

    def _tag_PH(self):
        pass

    def _tag_C2(self):
        pass

    def _tag_START(self):
        pass

    def _tag_C1(self):
        pass

    def _tag_C2DIS(self):
        pass

    def _tag_PX(self):
        pass

    def _tag_PLD(self):
        pass

    def _tag_ANNEXIN(self):
        pass

    def _tag_PLA(self):
        pass

    def _merge_datasets(self):
        pass
