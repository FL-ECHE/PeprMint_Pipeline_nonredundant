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

There are two main (proxy) methods for client classes: fetch() downloads most of
the data necessary, and superpose() runs cath-superpose to update downloaded pdb
files.

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
import stat
from pathlib import Path
import glob
import subprocess

import shutil
import urllib
import tarfile
import gzip

from src.settings import Settings

class DataRetriever:

    def __init__(self, global_settings: Settings):
        self.settings = global_settings
        self.UPDATE = False

    def fetch(self) -> bool:
        """ Public method #1
        Fetches the main data needed to construct the dataset, from CATH, EBI, 
        and UniProt. Not included in this module are AlphaFold PDBs (that were
        ported to a separate module, following Notebook #3) and the UniProt xml 
        sequence files (that are retrieved as part of the earlier architecture 
        of pepr2ds (called through DatasetManager).
        """

        if self.settings.FORMER_WORKING_DIR:
            print('Error: working directory already exists; remove it if you want to fetch all data from scratch with DataRetriever.fetch()')
            return False
        else:
            print(f'> Data retriever (domains selected in .config file: {self.settings.active_superfamilies})')
            self._retrieve_cath_domains()
            self._retrieve_uniprot_to_pdb_correspondence()
            self._retrieve_prosite()
            self._retrieve_cath_pdb_files()
            return True
    
    def superpose(self) -> bool:
        """ Public method #2
        Running the superposition over downloaded PDBs (cath-superpose 
        external tool, with built-in SSAP for all-vs-all pairwise structure 
        alignment)
        """
        print(f"> Superposing downloaded PDBs with cath-superpose and SSAP alignments")

        if not self._fetch_cath_superpose_binary():
            return False

        # TO DO: consider adding an option to avoid overwriting downloaded PDBs with superimposed ones
        if not self._run_cath_superpose():
            return False

        print(f"> Aligning superposed PDBs along z axis with respect to the domain representative")
        return self._align_on_z_axis()

    ############################################################################

    def _retrieve_cath_domains(self):
        # TO DO: release 4_2_0 or latest?
        self.dom_file = 'cath-domain-list.txt'
        url = self.settings.config_file['CATH']['domain_list_url']

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

    def _retrieve_uniprot_to_pdb_correspondence(self):
        url = self.settings.config_file['UNIPROT']['url']

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

    def _retrieve_prosite(self):
        url = self.settings.config_file['PROSITE']['url']

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

    def _retrieve_cath_pdb_files(self):
        # cath_domains: the complete list from cath with all domains
        self.cath_domains = pd.read_csv(self.settings.CATHFOLDER + self.dom_file,
                                        comment = '#',
                                        sep = r"\s+",
                                        header = None)
        self.cath_domains.columns = self.column_dom_file

        if self.settings.PARALLEL:
            self.cath_domains['Superfamily'] = self.cath_domains.parallel_apply(
                lambda x : f"{x.Class}.{x.Architecture}.{x.Topology}.{x.Homologous}",
                axis=1)
        else:
            self.cath_domains['Superfamily'] = self.cath_domains.progress_apply(
                lambda x : f"{x.Class}.{x.Architecture}.{x.Topology}.{x.Homologous}",
                axis=1)

        # cath_superfamily: just superfamily-domain correspondence
        self.cath_superfamily = pd.DataFrame()
        self.cath_superfamily['Superfamily'] = self.cath_domains.Superfamily
        self.cath_superfamily['Domain'] = self.cath_domains.Domain

        # cath_domains_per_superfamily: a map superfamily -> list of domains (from cath_superfamily) 
        # NB! do not parallelize this one
        self.cath_domains_per_superfamily = defaultdict(list)
        _ = self.cath_superfamily.progress_apply(
                lambda x : self.cath_domains_per_superfamily[x.Superfamily].append(x.Domain),
                axis = 1)

        # gather data for domains enabled in config file
        for superfamily, domain in self.settings.SUPERFAMILY.items():
            if domain in self.settings.active_superfamilies:
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

        # limit the dataset size when experimenting only
        if self.settings.XP_MODE and len(dom_list) > self.settings.xp_domain_limit:
                dom_list = dom_list[0:self.settings.xp_domain_limit]

                # make sure the reference pdb (to "align on z" later) is kept
                ref_pdb = self.settings.config_file['ALIGNMENT_ON_Z_AXIS']["ref_"+domName+"_pdb"]
                if ref_pdb is not None and ref_pdb not in dom_list:
                    dom_list[0] = ref_pdb

        if self.settings.PARALLEL:
            pd.Series(dom_list).parallel_apply(
                lambda x : self._fetch_pdb_from_cath_dom(x, folder) )
        else:
            print(dom_list)
            pd.Series(dom_list).progress_apply(
                lambda x : self._fetch_pdb_from_cath_dom(x, folder) )

    def _fetch_pdb_from_cath_dom(self, dom, folder):
        url = self.settings.config_file['CATH']['fetch_pdb_url'] +  dom + ".pdb"
        destination = folder + dom + '.pdb'
        if not os.path.isfile(destination): 
            urllib.request.urlretrieve(url, destination)

    ############################################################################

    def _fetch_cath_superpose_binary(self) -> bool:
        if self.settings.OS == "linux":
            tool_url = self.settings.config_file['CATH']['superpose_url_linux']
        elif self.settings.OS == "macos":
            tool_url = self.settings.config_file['CATH']['superpose_url_macos']
        else:
            print('Error: current superposition method (cath-superpose) only available on GNU/Linux and MacOS')
            return False

        self.superpose_tool = self.settings.CATHFOLDER + 'cath-superpose'
        if not os.path.isfile(self.superpose_tool): 
            urllib.request.urlretrieve(tool_url, self.superpose_tool)

        # make executable: 'current permission bits' BINARY OR 'executable bit'
        st = os.stat(self.superpose_tool)
        os.chmod(self.superpose_tool, st.st_mode | stat.S_IEXEC)

        return True

    def _run_cath_superpose(self, verbose=False) -> bool:
        success_only = True
        error_message = ""

        for domain in self.settings.active_superfamilies:
            # temporary output folders
            output_dir = self.settings.CATHFOLDER + 'domains/' + domain + '/super/'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            tmp_dir = self.settings.CATHFOLDER + 'domains/' + domain + '/cath_tools_tmp/'
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)

            # domain files as arguments to cath-superpose
            filelist = Path(self.settings.CATHFOLDER).glob(f"domains/{domain}/raw/*.pdb")
            run_input = ""
            for f in filelist:
                run_input += "--pdb-infile " + str(f) + " "

            # need to export environment var for cath-ssap
            raw_dir = self.settings.CATHFOLDER + 'domains/' + domain + "/raw"
            env_setting = f"export CATH_TOOLS_PDB_PATH={raw_dir} ; "

            # run cath-superpose
            run_options = "--do-the-ssaps " + tmp_dir
            run_output = "--sup-to-pdb-files-dir " + output_dir
            trigger = f"{self.superpose_tool} {run_options} {run_output} {run_input}"
            about = subprocess.run(env_setting+trigger, shell=True, capture_output=True)
            if verbose:
                print("done executing:")
                x = str(about.args)
                print(x)
                print("retuncode is:")
                x = int(about.returncode)
                print(str(x))
                print("stdout is:")
                x = str(about.stdout, "utf-8")
                print(x)
                print("stderr is")
                x = str(about.stderr, "utf-8")
                print(x)

            # TO DO: maybe call about.check_returncode() instead (raise exception)
            if int(about.returncode) != 0:
                success_only = False
                out = str(about.stderr, "utf-8")
                error_message += out

            # delete extra folder with alignments, replace raw with superposed
            shutil.rmtree(tmp_dir)
            raw_files = len([name for name in os.listdir(raw_dir) if os.path.isfile(os.path.join(raw_dir, name))])
            super_files = len([name for name in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, name))])
            if raw_files == super_files:
                shutil.rmtree(raw_dir)
                os.rename(output_dir, raw_dir)
            else:
                success_only = False
                out = f"Error: superposed file count for {domain} does not match raw file count\n"
                out += "Check directory: " + output_dir + "\n"
                error_message += out
        
        if not success_only:
            print("Errors found while executing cath-superpose:")
            print(error_message)

        return success_only

    def _align_on_z_axis(self):
        # TO DO: main() on old script becomes this method
        # TO DO: other functions on old script become accessory methods below
        # TO DO: replace current pdbs under 'raw' with reoriented ones

        """
        align_on_z.py -d PH      -ref 2da0A00 -res1 19   -res2 42   -res3 50   -i raw        -o zaligned
        align_on_z.py -d C2      -ref 1rsyA00 -res1 169  -res2 178  -res3 237  -i raw        -o zaligned
        align_on_z.py -d START   -ref 2e3mA00 -res1 412  -res2 448  -res3 515  -i raw        -o orientationA
        align_on_z.py -d START   -ref 2e3mA00 -res1 567  -res2 470  -res3 509  -i raw        -o orientationB
        align_on_z.py -d C1      -ref 1ptrA00 -res1 243  -res2 257  -res3 237  -i raw        -o zaligned
        align_on_z.py -d C2DIS   -ref 1czsA00 -res1 23   -res2 76   -res3 45   -i superposed -o zaligned
        align_on_z.py -d PX      -ref 1h6hA00 -res1 33   -res2 74   -res3 100  -i superposed -o zaligned
        align_on_z.py -d PLD     -ref 3rlhA00 -res1 59   -res2 205  -res3 198  -i superposed -o orientationOPM
        align_on_z.py -d PLD     -ref 3rlhA00 -res1 53   -res2 41   -res3 99   -i superposed -o orientationCAGE
        align_on_z.py -d ANNEXIN -ref 1a8aA01 -res1 25   -res2 68   -res3 77   -i superposed -o zaligned
        align_on_z.py -d PLA     -ref 1pocA00 -res1 7    -res2 92   -res3 76   -i superposed -o zaligned
        align_on_z.py -d SH2     -ref 2oq1A03 -res1 180  -res2 209  -res3 243  -i superposed -o zaligned
        align_on_z.py -d FYVE    -ref 1jocA02 -res1 1373 -res2 1392 -res3 1382 -i superposed -o zaligned
        align_on_z.py -d ENTH    -ref 1h0aA00 -res1 17   -res2 70   -res3 116  -i superposed -o zaligned
        """
        pass
