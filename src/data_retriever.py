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

The only difference here is that we moved notebook-related imports and initial 
settings to NotebookHandle. When using jupyter notebooks, it suffices to call 
the using_notebook() method from Settings before running everything else, e.g. 
as currently done in main().

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

    def run(self):
        self.retrieve_cath_domains()
        self.retrieve_uniprot_to_pdb_correspondence()
        self.retrieve_prosite()
        self.retrieve_cath_pdb_files()

    def retrieve_cath_domains(self):
        # TO DO: make these options instead of hardcoding, e.g. previous ftp address not working anymore (replaced by http)
        # TO DO: release 4_2_0 or latest? at least use CATHVERSION on Settings?
        self.dom_file = 'cath-domain-list.txt'
        url = "http://download.cathdb.info/cath/releases/all-releases/v4_2_0/cath-classification-data/cath-domain-list-v4_2_0.txt"

        destination = self.settings.CATHFOLDER + self.dom_file
        if not os.path.isfile(destination) or self.UPDATE: 
            urllib.request.urlretrieve(url, destination)

        self.column_dom_file = [ 'Domain',
                                 'Class',
                                 'Architecture',
                                 'Topology',
                                 'Homologous',
                                 'S35',
                                 'S60',
                                 'S95',
                                 'S100',
                                 'S100Count',
                                 'DomSize',
                                 'resolution', ]

    def retrieve_uniprot_to_pdb_correspondence(self):
        # TO DO: make this an option instead of hardcoding
        url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz"

        destination = self.settings.CATHFOLDER + "pdb_chain_uniprot.csv"

        # NB! made a small fix here (open/write operations were outside the if)
        if not os.path.isfile(destination) or self.UPDATE:
            tmp_destination = destination + ".gz"
            
            urllib.request.urlretrieve(url, tmp_destination)
            
            with gzip.open(tmp_destination, 'rb') as f:
                file = f.read()
                with open(destination, 'wb') as output:
                    output.write(file)
            os.remove(tmp_destination)

    def retrieve_prosite(self):
        # TO DO: make this an option instead of hardcoding
        url = "ftp://ftp.expasy.org/databases/prosite/prosite_alignments.tar.gz"

        destination = self.settings.PROSITEFOLDER + "prosite_alignments.tar.gz" 

        if self.UPDATE:
            if os.path.exists(self.settings.PROSITEFOLDER + "msa"):
                shutil.rmtree(self.settings.PROSITEFOLDER + "msa/")
            else:
                os.makedirs(self.settings.PROSITEFOLDER + "msa/")

        if not os.path.exists(self.settings.PROSITEFOLDER+"msa") or self.UPDATE:
            urllib.request.urlretrieve(url, destination)
            tf = tarfile.open(destination)
            tf.extractall(self.settings.PROSITEFOLDER)
            tf.close()
            os.rename(self.settings.PROSITEFOLDER + "prosite_alignments",
                      self.settings.PROSITEFOLDER + "msa")
            os.remove(destination)

    def retrieve_cath_pdb_files(self):
        # reading cath domain list
        self.cath_domains = pd.read_csv(self.settings.CATHFOLDER + self.dom_file,
                                        comment = '#',
                                        sep = r"\s+",
                                        header = None)
        self.cath_domains.columns = self.column_dom_file

        if self.PARALLEL:
            self.cath_domains['Superfamily'] = self.cath_domains.parallel_apply(
                lambda x : f"{x.Class}.{x.Architecture}.{x.Topology}.{x.Homologous}",
                axis=1)
        else:
            self.cath_domains['Superfamily'] = self.cath_domains.progress_apply(
                lambda x : f"{x.Class}.{x.Architecture}.{x.Topology}.{x.Homologous}",
                axis=1)

        # creating the superfamily
        self.cath_superfamily = pd.DataFrame()
        self.cath_superfamily['Superfamily'] = self.cath_domains.Superfamily
        self.cath_superfamily['Domain'] = self.cath_domains.Domain

        # creating a dictionary superfamily -> list of cathdomain (pdb format)
        # NB! do not parallelize this one
        self.cath_domains_per_superfamily = defaultdict(list)
        _ = self.cath_superfamily.progress_apply(
                lambda x : self.cath_domains_per_superfamily[x.Superfamily].append(x.Domain),
                axis = 1)

        # gather data for all domains
        for superfamily, domain in self.settings.SUPERFAMILY.items():
            print(f"> Fetching domain {domain} files")
            self._fetch_dom_for_superfamily(superfamily, domain)

    def _fetch_dom_for_superfamily(self, superfamily, domName):
        prefix = self.settings.CATHFOLDER
        folder = prefix + 'domains/' + domName + '/raw/'
        if not os.path.exists(folder):
            os.makedirs(folder)
        if not os.path.exists(prefix + 'domains/' + domName + '/cleaned/'):
            os.makedirs(prefix + 'domains/' + domName + '/cleaned/')

        dom_list = self.cath_domains_per_superfamily[superfamily]
        
        if self.PARALLEL:
            pd.Series(dom_list).parallel_apply(
                lambda x : self._fetch_pdb_from_cath_dom(x, folder) )
        else:
            print(dom_list)
            pd.Series(dom_list).progress_apply(
                lambda x : self._fetch_pdb_from_cath_dom(x, folder) )

    def _fetch_pdb_from_cath_dom(self, dom, folder):
        # TO DO: make this an option instead of hardcoding
        url = "http://www.cathdb.info/version/" + self.settings.CATHVERSION + "/api/rest/id/" + dom + ".pdb"
        destination = folder + dom + '.pdb'
        if not os.path.isfile(destination): 
            urllib.request.urlretrieve(url, destination)
