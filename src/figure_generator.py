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


""" Methods to generate figures for data analyses (previously on Notebook #4)

In particular, this implementation was used to generate figures in the 2022 
paper "Dissecting peripheral protein-membrane interfaces" by Tubiana, Sillitoe, 
Orengo, Reuter: https://doi.org/10.1371/journal.pcbi.1010346

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
from src.dataset_manager import DatasetManager
from src.notebook_handle import NotebookHandle

class FigureGenerator:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings
        self.dataset = None

        self.RECALCULATE = self.settings.config_file.getboolean(
            'FIGURE_GENERATION', 'recalculate')
        self.INCLUDE_AF_FROM_START = self.settings.config_file.getboolean(
            'FIGURE_GENERATION', 'include_alphafold_from_the_beginning')

        self.FILESUFFIX = "-AF" if self.INCLUDE_AF_FROM_START else ""
        
        self._libs_setup()

        # TO DO: continue from here when done with tools/tag_IBS_domains
        # - Load/recalculate data
        # - Preparation/cleaning
        # - Palette and conf

    def _libs_setup(self):
        # IPython
        if self.settings.USING_NOTEBOOK:
            self.settings.NOTEBOOK_HANDLE.figure_generator_options()
        
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

    """
    ### All methods below just encapsulate the steps in Notebook #4
    """

    # Figure 2 in the 2022 paper: composition of the exposed IBS for all proteins
    def make_figure_composition_of_exposed_IBS(self):
        pass

    # Figure 3 in the 2022 paper: protrusions and hydrophobic protrusions in the dataset
    def make_figure_protrusions(self):
        pass

    # corresponding statistical tests
    def statistical_tests_on_protusion(self):
        pass

    # percentage of the dataset without HP
    def percentage_without_HP(self):
        pass

    # Figure 4 in the 2022 paper: composition of the IBS/nonIBS for protein with HP at IBS
    def make_figure_composition_for_proteins_with_HP_at_IBS(self):
        pass

    # Figure 5 in the 2022 paper: neighbourhood composition
    def make_figure_neighbourhood_composition(self):
        pass

    # Figure 6 in the 2022 paper: number of structures with/without HP at IBS and comparison of both datasets
    def make_figure_number_of_structures_w_and_wo_HP_at_IBS(self):
        pass

    # raw values and statistical test
    def statistical_tests_on_HP_at_IBS(self):
        pass

    # Figure 7 in the 2022 paper: composition of the IBS and the protrusion neighbourhood for proteins WITHOUT HP AT IBS
    def make_figure_composition_for_proteins_without_HP_at_IBS(self):
        pass

    # Preparation of AlphaFold data insertion
    def prepare_alphafold_data_insertion(self):
        pass

    # Figure 8 in the 2022 paper: superfamily decomposition of Exposed environment of hydrophobic protrusions at the IBS
    def make_figure_superfamily_decomposition(self):
        pass

    # Table 1: data collection and processing
    def make_table_data_collection_and_processing(self):
        pass

    # Table 2: classification of amino acids
    def make_table_classification_of_aminoacids(self):
        pass

    # Table 3: description of the properties observed (A) among sub-dataset (B) used in the odds ratio calculation (equation X).
    def make_table_description_of_properties(self):
        pass

    # Table 4: p-values for wilcoxon signed rank test
    def make_table_p_values_for_wilcoxon_test(self):
        pass

    # Supplement Figure 1: composition of the exposed IBS (DATASET augmented)
    def make_figure_composition_of_exposed_IBS_on_augmented_dataset(self):
        pass

    # Supplement Figure 2: amino acid composition per secondary structures
    def make_figure_aminoacid_composition_per_secondary_structures(self):
        pass

    # Supplement Figure 3: odds ratio for secondary structures
    def make_figure_odds_ratio_for_secondary_structures(self):
        pass

    # Supplement Figure 4: increasing our dataset size by integrating alphafold2 models
    def make_figure_dataset_size_increase_with_alphafold(self):
        pass

    # Supplement Figure 5: number of structures per domains with alphafold2
    def make_figure_structures_per_domains_with_alphafold(self):
        pass

    # Supplement Table 1: origin of the data
    def make_table_origin_of_data(self):
        pass

    # Supplement Figure 6: hydrophobic protrusions exposed environment at the IBS?
    def make_figure_hydro_protrusions_per_domain_at_IBS(self): 
        pass

    # Export dataset for journal
    def export_dataset(self):
        pass

    # Export format file for journal
    def export_format_file(self):
        pass
