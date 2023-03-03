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
import re

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import MDAnalysis as mda
from scipy.spatial import ConvexHull

from biopandas.pdb import PandasPdb
from Bio import AlignIO

import os
import glob
import json
import urllib
from urllib.error import HTTPError, URLError
import requests

from tqdm import tnrange, tqdm
from typing import Optional

from src.settings import Settings
from src.notebook_handle import NotebookHandle

class AlphaFoldUtils:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings
        self.dataset = None

        self._libs_setup()

        self.REBUILD = self.settings.config_file.getboolean(
            'ALPHAFOLD_UTILS', 'rebuild')
        self.use_uniprot_boundaries = self.settings.config_file.getboolean(
            'ALPHAFOLD_UTILS', 'use_uniprot_boundaries')
        self.use_all_AFmodels = self.settings.config_file.getboolean(
            'ALPHAFOLD_UTILS', 'use_all_AFmodels')

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

    def run(self,
            df: pd.DataFrame,
            EXCLUDE_LIST: Optional[list] = None,
            EXCLUDE_DOMAIN: Optional[list] = None):
        # runs each step to fetch and prepare alphafold data (as of Notebook #3)
        self.dataset = df

        self.REGEX = re.compile("^(\w+)\|(\w+)\/(\d+)-(\d+)")

        print(f"Preparing AlphaFold data")
        if EXCLUDE_LIST is not None:
            print("User option: excluding list ", end='')
            print(EXCLUDE_LIST)
        else:
            EXCLUDE_LIST = []

        if EXCLUDE_DOMAIN is not None:
            print("User option: excluding domain ", end='')
            print(EXCLUDE_DOMAIN)
        else:
            EXCLUDE_DOMAIN = []

        self.process_all_domains(EXCLUDE_LIST, EXCLUDE_DOMAIN)
        self.align_fasta_files()
        self.save_PH_domains_from_alphafold()

    """
    ### All methods below just encapsulate the steps in Notebook #3
    """

    def process_all_domains(self, EXCLUDE_LIST, EXCLUDE_DOMAIN):
        # TO DO: switch these assignments back
        #domains = self.dataset.domain.unique()
        domains = ['PLA']

        for domain in domains:
            if domain in EXCLUDE_DOMAIN:
                continue

            print(f"----- PROCESSING DOMAIN {domain} -----")

            group = self.dataset.query("domain == @domain")

            uniprot_acc_cathpdb = group.query("data_type == 'cathpdb'").uniprot_acc.unique()

            seqs_no_pdb = group[group["pdb"].isnull()].uniprot_acc.unique()

            boundaries_prosite = self._get_prosite_boundaries_dict(domain)

            if self.use_all_AFmodels:
                prosite_uniprot_acc = list(boundaries_prosite.keys()) 
                uniprot_acc_cathpdb = [acc for acc in uniprot_acc_cathpdb if acc in prosite_uniprot_acc]

                uniprot_acc_list = prosite_uniprot_acc + uniprot_acc_cathpdb

                seqs_with_model, seqs_without_model = self._fetch_pdb_alphafold(uniprot_acc_list, 
                                                                                domain)
            else:
                seqs_with_model, seqs_without_model = self._fetch_pdb_alphafold(seqs_no_pdb, 
                                                                                domain)

            for uniprot_id in tqdm(seqs_with_model, desc = "processing"):
                if uniprot_id in EXCLUDE_LIST:
                    continue

                pdbfile =  f"{self.settings.ALPHAFOLDFOLDER}/{domain}/raw/{uniprot_id}.pdb"
                # structure = PDBParser().get_structure('uniprot_id',)

                if os.path.isfile(pdbfile) and self.REBUILD == False:
                    continue   # skip the file if it already exists

                query = self._get_domain_fragment_query(uniprot_id, domain, boundaries_prosite)
                if query == None:
                    continue

                ppdb = PandasPdb().read_pdb(pdbfile)
                ppdb.df["ATOM"] = ppdb.df["ATOM"].query(f"{query}")
                ppdb.to_pdb(f"{self.settings.ALPHAFOLDFOLDER}/{domain}/extracted/{uniprot_id}.pdb")

    def _get_prosite_boundaries_dict(self, domain):
        boundaries = {}
        prosite_ids = self.settings.DOMAIN_PROSITE[domain]

        if type(prosite_ids) != type([]):
            prosite_ids = [prosite_ids]

        for msafile in prosite_ids:
            msafilepath = f"{self.settings.PROSITEFOLDER}/msa/{msafile}.msa"
            msa = AlignIO.read(msafilepath,'fasta')
            for record in msa:
                seqid = record.id
                match = self.REGEX.match(seqid)
                if match:
                    uniprot_id = match.group(2)
                    start = match.group(3)
                    end = match.group(4)
                    boundaries[uniprot_id] = (int(start),int(end))
        return boundaries

    def _fetch_pdb_alphafold(self, uniprotids, domain):
        nomodels = []
        withmodels = []

        outfolder = f"{self.settings.ALPHAFOLDFOLDER}/{domain}/raw"
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        
        extractedfolder = f"{self.settings.ALPHAFOLDFOLDER}/{domain}/extracted"
        if not os.path.exists(extractedfolder):
            os.makedirs(extractedfolder)
        else:
            if self.REBUILD == True:   #delete extracted files
                files = glob.glob(f"{extractedfolder}/*.pdb")
                for f in files:
                    os.remove(f)
        
        jsonfolder = f"{self.settings.ALPHAFOLDFOLDER}/{domain}/json"
        if not os.path.exists(jsonfolder):
            os.makedirs(jsonfolder)

        prefix = self.settings.config_file['ALPHAFOLD_UTILS']['AF_pdbs_url_prefix']
        suffix = self.settings.config_file['ALPHAFOLD_UTILS']['AF_pdbs_url_suffix']

        for uniprot_id in tqdm(uniprotids, desc="Downloading "):
            url = prefix + uniprot_id + suffix
            destination = f"{outfolder}/{uniprot_id}.pdb"
            if not os.path.isfile(destination): 
                try:
                    urllib.request.urlretrieve(url, destination)
                except urllib.error.HTTPError as err:
                    nomodels.append(uniprot_id)
                    continue
            withmodels.append(uniprot_id)

        rate = len(nomodels)/len(uniprotids) if len(uniprotids) > 0 else 0
        print(f"{len(nomodels)} out of {len(uniprotids)} without alphafold2 models ({rate*100:.2f}%)")
        return withmodels, nomodels

    def _get_domain_fragment_query(self, uniprot_acc, domain, boundaries_prosite):
        start_PS, end_PS = boundaries_prosite[uniprot_acc]
        starts_ends = [boundaries_prosite[uniprot_acc]]

        if self.settings.DOMAIN_INTERPRO_REFINE[domain] == True:
            if domain == "PLA":
                source = 'cathgene3d'
            else:
                source = 'ssf'

            interpro = self._get_json(uniprot_acc, domain, source)
            
            if interpro == None:
                return None

            QueryString = None
            
            for result in interpro["results"]:
                if result["metadata"]["accession"] == self.settings.DOMAIN_INTERPRO[domain]:
                    entry_protein_locations = result["proteins"][0]["entry_protein_locations"]
                    for entry in entry_protein_locations:
                        # get the number of truncation in the domain
                        nfrag = len(entry['fragments'])
                        
                        if domain == 'PLA':
                            # special case for PLA, we will ignore PROSITE annotation that are actually wrong
                            frag = entry['fragments'][0]   # get first monomer only
                            s = entry['fragments'][0].get('start')
                            e = entry['fragments'][0].get('end')
                            starts_ends = [[s,e]]
                        else:
                            if nfrag >= 2 and ( entry['fragments'][0].get('start') - 50 <= start_PS <= entry['fragments'][0].get('start')+50):
                                # if truncated domain AND correspond to the prosite domain
                                print(f"splitting {domain}-{uniprot_acc}")
                                queries = []
                                starts_ends = []
                                for frag in entry['fragments']:
                                    s = int(frag.get("start"))
                                    e = int(frag.get("end"))
                                    starts_ends.append([s,e])
                                if use_uniprot_boundaries == True:
                                    starts_ends[0][0] = start_PS
                                    starts_ends[-1][-1] = end_PS
                            else:
                                # use prosite fragment
                                starts_ends = [[start_PS, end_PS]]

                    QueryString = " or ".join([f"({x} <= residue_number <= {y})" for x,y in starts_ends])
        else:
            QueryString = " or ".join([f"({x} <= residue_number <= {y})" for x,y in starts_ends])
        
        return QueryString

    def _get_json(self, uniprot_acc, domain, source='ssf'):
        jsonfolder = f"{self.settings.ALPHAFOLDFOLDER}/{domain}/json"
        if not os.path.exists(jsonfolder):
            os.makedirs(jsonfolder)
        
        jsonfile = f"{jsonfolder}/{uniprot_acc}.json"
        if os.path.isfile(jsonfile):
            f = open(jsonfile)
            interpro = json.load(f)
        else:
            #make the query on ebi/interpro
            url_prefix = self.settings.config_file['ALPHAFOLD_UTILS']['AF_interpro_url_prefix']
            url_middle = self.settings.config_file['ALPHAFOLD_UTILS']['AF_interpro_url_middle']
            url_suffix = self.settings.config_file['ALPHAFOLD_UTILS']['AF_interpro_url_suffix']
            url = url_prefix + source + url_middle + uniprot_acc + url_suffix

            response = self.__request_URL(url)
            if response == None:
                return None
            try:
                interpro = json.loads(response)
            except:
                print(f"no data for {uniprot_acc}")
                return None
            with open(jsonfile,'w') as out:
                json.dump(interpro, out, indent=2)
                
        return(interpro)

    def __request_URL(self, link, trial=1):
        try:
            response = requests.get(link).text
            return response
        except URLError as e:
            print(e, link)
            if trial > 3 :
                print('3rd fail, skipping this one')
                return None
            else:
                print(f"Trial {trial}, waiting 10s and trying again")
                sleep(10)
                return self.__request_URL(link, trial=trial+1)

    def align_fasta_files(self):
        pass

    def save_PH_domains_from_alphafold(self):
        pass
