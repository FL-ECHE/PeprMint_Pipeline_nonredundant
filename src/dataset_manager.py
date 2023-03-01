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


""" PePrMInt dataset creation (previously on Notebook #2)

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
        self._FULL_DATASET_FILENAME = "DATASET_peprmint_allatoms_d25"
        self._LIGHT_DATASET_FILENAME = "DATASET_peprmint_d25"

        import pepr2ds.builder.Builder as builderEngine
        importlib.reload(builderEngine)
        self.builder = builderEngine.Builder(self.settings.SETUP, 
                                             recalculate = self.RECALCULATION,
                                             update = False,
                                             notebook = self.settings.USING_NOTEBOOK,
                                             core = 1)

    def load_full_dataset(self) -> bool:
        # loads the dataset built in a previous run (this full version is NOT 
        # used in the original notebooks)
        full_path = self.settings.WORKDIR + self._FULL_DATASET_FILENAME + ".pkl"
        if not os.path.isfile(full_path):
            print("Could not find the full dataset file ({0})".format(full_path))
            print("Use build() from DatasetManager to create it")
            return False
        else:
            self.DATASET = pd.read_pickle(full_path)
            print("Dataset (full version) loaded successfully")
            return True

    def load_light_dataset(self) -> bool:
        # loads the dataset built in a previous run (this light version is the 
        # one used throughout the original notebooks)
        full_path = self.settings.WORKDIR + self._LIGHT_DATASET_FILENAME + ".pkl"
        if not os.path.isfile(full_path):
            print("Could not find the light dataset file ({0})".format(full_path))
            print("Use build() from DatasetManager to create it")
            return False
        else:
            self.DATASET = pd.read_pickle(full_path)
            print("Dataset (light version) loaded successfully")
            return True

    def build(self, recalculate: Optional[bool] = None):
        # runs each step in the creation of the dataset (as of Notebook #2)
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

        print("\nDataset built successfully")
        print("Dataset domains: ")
        print(list(self.DATASET.domain.unique()))

    """
    ### All methods below just encapsulate the steps in Notebook #2
    """

    def clean(self):
        self.builder.structure.clean_all_pdbs()
        self.DATASET = self.builder.structure.build_structural_dataset()
        self.DATASET.data_type.unique()

    def compute_protusion(self):
        self.DATASET = self.builder.structure.add_protrusions(self.DATASET)

    def add_cluster_structural_info(self):
        # TO DO: can we avoid recreating the builder object!?
        import pepr2ds.builder.Builder as builderEngine
        importlib.reload(builderEngine)
        self.builder = builderEngine.Builder(self.settings.SETUP, 
                                             recalculate = self.RECALCULATION, 
                                             update=False, 
                                             notebook = self.settings.USING_NOTEBOOK, 
                                             core=1)
        self.DATASET = self.builder.structure.add_structural_cluster_info(self.DATASET)

    def add_uniprot_basic_info(self):
        self.DATASET = self.builder.sequence.add_uniprotId_Origin(self.DATASET)

    def add_prosite_info(self):
        self.DATASET = self.builder.sequence.match_residue_number_with_alignment_position(self.DATASET)

    def add_sequences_without_structure(self):
        self.DATASET = self.builder.sequence.add_sequence_in_dataset(self.DATASET)

    def add_uniprot_protein_sheet_info(self):
        self.builder.sequence.download_uniprot_data(self.DATASET)
        self.DATASET = self.builder.sequence.add_info_from_uniprot(self.DATASET)

    def add_cluster_uniref_info(self):
        self.DATASET = self.builder.sequence.add_cluster_info(self.DATASET)

    def add_conservation(self):
        # TO DO: can we avoid recreating the builder object!?
        import pepr2ds.builder.Builder as builderEngine
        importlib.reload(builderEngine)
        self.builder = builderEngine.Builder(self.settings.SETUP, 
                                             recalculate = self.RECALCULATION, 
                                             notebook = self.settings.USING_NOTEBOOK)

        self.DATASET = self.builder.sequence.add_conservation(self.DATASET,
                                                              gapcutoff=0.8)

    def save_dataset(self):
        self.DATASET = self.builder.optimize_size(self.DATASET)
        self.DATASET = self.DATASET.drop_duplicates(subset=['atom_number',
                                                            'atom_name',
                                                            'residue_name',
                                                            'residue_number',
                                                            'cathpdb',
                                                            'chain_id'])

        # save two versions of the dataset: a complete one (with all PDB atoms)
        # and a light one (only with `CA` and `CB` atoms).  
        self.builder.save_checkpoint_dataset(self.DATASET,
                                             self._FULL_DATASET_FILENAME)
        self.builder.save_checkpoint_dataset(self.DATASET.query("atom_name in ['CA','CB']"),
                                             self._LIGHT_DATASET_FILENAME)

    # TO DO: does not seem meaningful; consider removing in the future
    def test_alignment_for_C2DIS_domain(self):
        # generate alignment file list for C2DIS domain (S95, because everything is just too slow, too much structure)
        c2dis = self._selectUniquePerCluster(self.DATASET.query("domain== 'PH'"), 
                                             'S95', 
                                             'uniref90', 
                                             withAlignment=False, 
                                             pdbreference='2da0A00')
        #pdblist = c2dis.cathpdb.dropna().unique()
        #print(c2dis)
        print(self.DATASET.query("atom_name == 'CA' and domain =='PH'").columns)
        #print(self.DATASET.query("atom_name == 'CA' and domain =='PH' and data_type == 'cathpdb'")[["ASA_res_freesasa_florian","RSA_freesasa_florian","ASA_total_freesasa","ASA_mainchain_freesasa","ASA_sidechain_freesasa","RSA_sidechain_freesasa","RSA_total_freesasa_tien","RSA_sidechain_freesasa_tien"]])

    # TO DO: does not seem necessary; move to an ad-hoc module
    def _selectUniquePerCluster(self,
                                df, 
                                cathCluster, 
                                Uniref, 
                                withAlignment=True, 
                                pdbreference=None,
                                removeStrand=False):
        # return a dataset with only 1 datum per choosed cluster
        if cathCluster not in ["S35", "S60", "S95", "S100"]:
            raise ValueError('CathCluster given not in ["S35","S60","S95","S100"]')

        if Uniref not in ["uniref50", "uniref90", "uniref100"]:
            raise ValueError('CathCluster given not in ["uniref50","uniref90","uniref100"]')

        if withAlignment:
            df = df[~df.alignment_position.isnull()]

        cathdf = df.query("data_type == 'cathpdb'")
        seqdf = df.query("data_type == 'prosite'")

        def selectUniqueCath(group):
            uniqueNames = group.cathpdb.unique()
            if pdbreference:
                if pdbreference in uniqueNames:
                    select = pdbreference
                else:
                    select = uniqueNames[0]
            else:
                select = uniqueNames[0]

            # return group.query("cathpdb == @select")
            return select

        def selectUniqueUniref(group, exclusion):
            uniqueNames = group.uniprot_acc.unique()
            select = uniqueNames[0]
            # return group.query("uniprot_acc == @select")
            if select not in exclusion:
                return select

        # structures prior to sequences
        dfReprCathNames = cathdf.groupby(["domain", cathCluster]).apply(selectUniqueCath).to_numpy()
        print(dfReprCathNames)
        excludeUniref = df.query("cathpdb in @dfReprCathNames").uniprot_acc.unique()

        dfReprUnirefNames = seqdf.groupby(["domain", Uniref]).apply(selectUniqueUniref,exclusion=excludeUniref).to_numpy()
        dfReprCath = cathdf.query("cathpdb in @dfReprCathNames")
        dfReprUniref = seqdf.query("uniprot_acc in @dfReprUnirefNames")

        return (pd.concat([dfReprCath, dfReprUniref]))
