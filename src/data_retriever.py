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
import numpy as np
import vg
from collections import defaultdict
from tqdm import trange, tqdm

from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import os
import stat
from pathlib import Path
import glob
import subprocess
import warnings

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

        print(f"> Orienting superposed PDBs along z axis with respect to the domain representative")
        return self._orient_on_z_axis()

    ############################################################################
    # fetch-related methods below

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
                ref_pdb = self.settings.config_file['REORIENT_ALONG_Z']["ref_"+domName+"_pdb"]
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
    # superpose-related methods below

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
    
    ############################################################################
    # reorient-related methods below

    def _orient_on_z_axis(self) -> bool:
        success_only = True
        error_message = ""

        for domain in self.settings.active_superfamilies:
            ref_pdb = self.settings.config_file['REORIENT_ALONG_Z']["ref_"+domain+"_pdb"]
            if ref_pdb is None:
                print(f"  Warning: no reference protein for '{domain}' defined on .config file - skipping reorientation along z")
                continue
            
            # read pdbs in raw, save in temporary output folder, then replace them
            raw_dir = self.settings.CATHFOLDER + 'domains/' + domain + "/raw"
            output_dir = self.settings.CATHFOLDER + 'domains/' + domain + '/z_oriented/'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            print(f"  '{domain}' orientation:")
            rotation, translation = self._get_transformation_from_reference(domain)
            print(f"  rotation matrix: {rotation}")
            print(f"  translation vector: {translation}")
            self._transform_pdbs(domain, rotation, translation)

            # replace raw with reoriented
            raw_files = len([name for name in os.listdir(raw_dir) if os.path.isfile(os.path.join(raw_dir, name))])
            oriented_files = len([name for name in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, name))])
            if raw_files == oriented_files:
                shutil.rmtree(raw_dir)
                os.rename(output_dir, raw_dir)
            else:
                success_only = False
                out = f"Error: reoriented file count for '{domain}' does not match raw file count\n"
                out += "Check directory: " + output_dir + "\n"
                error_message += out
        
        if not success_only:
            print("Errors found while orienting along z axis:")
            print(error_message)

        return success_only

    def _get_transformation_from_reference(self, domain):
        raw_dir = self.settings.CATHFOLDER + 'domains/' + domain + "/raw"
        ref_pdb = self.settings.config_file['REORIENT_ALONG_Z']["ref_"+domain+"_pdb"]
        structure =  self._get_structure(f"{raw_dir}/{ref_pdb}.pdb")

        res1 = self.settings.config_file.getint('REORIENT_ALONG_Z', 'ref_'+domain+'_res1')
        res2 = self.settings.config_file.getint('REORIENT_ALONG_Z', 'ref_'+domain+'_res2')
        res3 = self.settings.config_file.getint('REORIENT_ALONG_Z', 'ref_'+domain+'_res3')

        translation = self._get_translation_vector(structure, res1, res2, res3)
        
        structure = self._apply_translation_vector(structure, translation)
        rotation = self._get_rotation_matrix(structure, res1, res2, res3, domain)

        return (rotation,translation)

    def _transform_pdbs(self,
                        domain,
                        rotation: np.array,       # rotation matrix (3x3)
                        translation: np.array):   # translation vector (1x3)

        raw_dir = self.settings.CATHFOLDER + 'domains/' + domain + "/raw"
        output_dir = self.settings.CATHFOLDER + 'domains/' + domain + '/z_oriented/'
        pdb_list = glob.glob(f"{raw_dir}/*.pdb")

        for pdb_path in tqdm(pdb_list):
            structure = self._get_structure(pdb_path)
            structure = self._apply_translation_vector(structure,translation)
            structure = self._apply_rotation_matrix(structure, rotation)

            # write new file
            pdb_io = PDBIO()
            pdb_io.set_structure(structure)
            pdb_io.save(f"{output_dir}/{structure.id}.pdb")

    def _get_structure(self, pdb_path):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', PDBConstructionWarning)

                parser = PDB.PDBParser()
                pdb_code = Path(pdb_path).stem   # filename without extension
                return parser.get_structure(id = pdb_code, file = pdb_path)
        except:
            print(f"Reading Error: {pdb_path}")

    def _get_translation_vector(self,
                                structure: PDB.Structure.Structure,   # Biopython PDB Structure
                                res1: int,
                                res2: int,
                                res3: int) -> np.array:
        # Get the translation vector between the origin and the centroid of the
        # triangle formed by the CA of three aminoacid residues
        chain = structure[0].child_list[0].id
        p1 = structure[0][chain][res1]['CA'].get_coord()
        p2 = structure[0][chain][res2]['CA'].get_coord()
        p3 = structure[0][chain][res3]['CA'].get_coord()

        # translation vector (1x3)
        return list(-np.mean(np.array([p1, p2, p3]), axis=0))

    def _apply_translation_vector(self,
                                  structure: PDB.Structure.Structure,   # Biopython PDB Structure
                                  translation: np.array):   # (1x3) translation vector
        # apply translation without rotation
        null_rotation = np.identity(3).tolist()

        structure[0].transform(null_rotation, translation)
        return structure

    def _apply_rotation_matrix(self,
                               structure: PDB.Structure.Structure,   # Biopython PDB Structure
                               rotation:np.array):   # (3x3) rotation matrix
        # apply rotation without translation
        structure[0].transform(rotation, [0, 0, 0])
        return structure

    def _get_rotation_matrix(self,
                             structure: PDB.Structure.Structure,   # Biopython PDB Structure
                             res1: int,
                             res2: int,
                             res3: int,
                             domain: str,
                             orientation = 'z') -> np.array:       # axis used for alignment
        # Get the rotation matrix between the normal of a triangle formed by the CA 
        # of three aminoacid residues and a plane (x,y,z)

        def get_normal_COM(structure, res1, res2, res3):
            # normal vector and and the geom. center of a structure

            # compute vectors (get new coordinates - rotation with null translation)
            chain = structure[0].child_list[0].id
            p1 = structure[0][chain][res1]['CA'].get_coord()
            p2 = structure[0][chain][res2]['CA'].get_coord()
            p3 = structure[0][chain][res3]['CA'].get_coord()

            # compute normal to the triangle: the cross product of two vectors
            A = p2 - p1
            B = p3 - p1
            N = np.cross(A, B)

            coords = np.array([x.get_coord() for x in structure.get_atoms()])
            COM = coords.mean(axis=0)

            return N, COM

        def test_rotation(structure, res1, res2, res3):
            # Recalculate angle etc...
            N, COM = get_normal_COM(structure, res1, res2, res3)

            # Recalculate angle
            angle = vg.angle(N, np.array([0, 0, -1]))
            return (angle == 0 and COM[2] > 0) or (angle == 180 and COM[2] > 0)

        N,COM = get_normal_COM(structure, res1, res2, res3)

        # This norm will be our translation vector to all our atoms
        axis = {'x':[-1,0,0],
                'y':[0,-1,0],
                'z':[0,0,-1]}
        #Create the reference vector, per default we want to align on the z axis so it will be [0,0,1]
        refVector = PDB.vectors.Vector(axis[orientation])
        normal = PDB.vectors.Vector(N) #The normal should be a biopython object

        # Transformation 1
        temp = structure.copy()

        #case1 : normal and rotation
        rotation = PDB.vectors.rotmat(normal, refVector)
        temp.transform(rotation, [0, 0, 0])

        #If it doesn't work, case2: -normal and rotation
        if not test_rotation(temp, res1, res2, res3):
            temp = structure.copy()
            rotation = PDB.vectors.rotmat(-normal, refVector)
            temp.transform(rotation, [0, 0, 0])
            # If it doesn't work, case3: normal and rotation.T
            if not test_rotation(temp, res1, res2, res3):
                temp = structure.copy()
                rotation = PDB.vectors.rotmat(normal, refVector).T
                temp.transform(rotation, [0, 0, 0])
                # If it doesn't work, case2: -normal and rotation.T
                if not test_rotation(temp, res1, res2, res3):
                    temp = structure.copy()
                    rotation = PDB.vectors.rotmat(-normal, refVector).T

        # (3x3) rotation matrix
        return rotation
