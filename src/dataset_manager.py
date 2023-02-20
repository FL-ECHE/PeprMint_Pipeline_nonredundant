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


""" PePrMInt dataset creation (previosuly on Notebook #2)

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
import os
import importlib

from tqdm import trange, tqdm
from typing import Optional

from src.settings import Settings
from src.notebook_handle import NotebookHandle

class DatasetManager:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings

        self.RECALCULATION = self.settings.config_file.getboolean(
            'DATASET_MANAGER', 'recalculate')

        self._libs_setup()
        self._pepr2ds_setup()

    def _libs_setup(self):
        # Pandas
        pd.options.mode.chained_assignment = (
            None  # remove warning when adding a new column; default='warn'
        )
        pd.set_option("display.max_columns", None)
        tqdm.pandas()   # activate tqdm progressbar for pandas

        # IPython
        if self.settings.USING_NOTEBOOK:
            self.settings.NOTEBOOK_HANDLE.dataset_manager_options()

    def _pepr2ds_setup(self):
        import pepr2ds.builder.Builder as builderEngine
        importlib.reload(builderEngine)
        self.builder = builderEngine.Builder(self.settings.SETUP, 
                                             recalculate = self.RECALCULATION,
                                             update = False,
                                             notebook = self.settings.USING_NOTEBOOK,
                                             core = 1)

    def run(self, recalculate: Optional[bool] = None):
        if recalculate is not None:
            self.RECALCULATION = recalculate
        self.clean()
        self.compute_protusion()
        self.add_cluster_structural_info()
        self.add_uniprot_basic_info()
        self.add_prosite_info()
        self.add_sequences_without_structure()
        self.add_uniprot_protein_sheet_info()
        self.add_cluster_uniref_info()
        self.add_conservation()
        self.save_dataset()

    def clean(self):
        self.builder.structure.clean_all_pdbs()
        self.DATASET = self.builder.structure.build_structural_dataset()
        self.DATASET.data_type.unique()

    def compute_protusion(self):
        pass

    def add_cluster_structural_info(self):
        pass

    def add_uniprot_basic_info(self):
        pass

    def add_prosite_info(self):
        pass

    def add_sequences_without_structure(self):
        pass

    def add_uniprot_protein_sheet_info(self):
        pass

    def add_cluster_uniref_info(self):
        pass

    def add_conservation(self):
        pass

    def save_dataset(self):
        pass
