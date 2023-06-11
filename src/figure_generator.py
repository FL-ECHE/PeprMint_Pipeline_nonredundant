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
from decimal import Decimal
import scipy.stats as stats

from tqdm import tnrange, tqdm

import urllib
import glob
from urllib.error import HTTPError
from pathlib import Path

from src.settings import Settings
from src.notebook_handle import NotebookHandle
from pepr2ds.dataset.tagibs import Dataset

class FigureGenerator:

    def __init__(self, global_settings: Settings, tagged_dataset: Dataset):
        self.settings = global_settings
        self.pepr2ds = tagged_dataset

        self.INCLUDE_AF_FROM_START = self.settings.config_file.getboolean(
            'FIGURE_GENERATION', 'include_alphafold_from_the_beginning')

        self.FILESUFFIX = "-AF" if self.INCLUDE_AF_FROM_START else ""
        self._silent_eq_test = False

        self._libs_setup()
        self._palette_setup()
        self._data_setup()

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
        sns.set_style("whitegrid")

        # Matplotlib
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams["font.family"] = "Arial"

    def _palette_setup(self):
        self.dict_palette_IBS = {"IBS":"#A5CCE2",'nonIBS':"#DFA3AE"}
        self.palette_IBS = [self.dict_palette_IBS["nonIBS"],self.dict_palette_IBS["IBS"]]
        self.FORMAT="tif"
        
        self.colorsPerType = {"Positive":"tab:blue",
                 "Negative":"tab:red",
                 "Non-polar":"tab:gray",
                 "Hydrophobic,H-non-aromatic":"tab:brown",
                 "Hydrophobic,H-aromatic":"tab:pink",
                 "Polar":"tab:green"}

        #From http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm#shapelycolours
        self.COLORS_taylor = {
            "LEU": "#33FF00",
            "ILE": "#66FF00",
            "CYS": "#FFFF00",
            "MET": "#00FF00",
            "TYR": "#00FFCC",
            "TRP": "#00CCFF",
            "PHE": "#00FF66",
            "HIS": "#0066FF",
            "LYS": "#6600FF",
            "ARG": "#0000FF",
            "ASP": "#FF0000",
            "GLU": "#FF0066",
            "VAL": "#99FF00",
            "ALA": "#CCFF00",
            "GLY": "#FF9900",
            "PRO": "#FFCC00",
            "SER": "#FF3300",
            "ASN": "#CC00FF",
            "GLN": "#FF00CC",
            "THR": "#FF6600",
            "UNK": "#000000"
        }

    def _data_setup(self):
        # Thibault: "Temporary fix: redefine type to have HIS as polar"
        self._AATYPE = {
            "LEU": "Hydrophobic,H-non-aromatic",
            "ILE": "Hydrophobic,H-non-aromatic",
            "CYS": "Hydrophobic,H-non-aromatic",
            "MET": "Hydrophobic,H-non-aromatic",
            "TYR": "Hydrophobic,H-aromatic",
            "TRP": "Hydrophobic,H-aromatic",
            "PHE": "Hydrophobic,H-aromatic",
            "HIS": "Polar",
            "LYS": "Positive",
            "ARG": "Positive",
            "ASP": "Negative",
            "GLU": "Negative",
            "VAL": "Non-polar",
            "ALA": "Non-polar",
            "SER": "Polar",
            "ASN": "Polar",
            "GLY": "Non-polar",
            "PRO": "Non-polar",
            "GLN": "Polar",
            "THR": "Polar",
            "UNK": "none"
        }
        self.pepr2ds.domainDf["type"] = self.pepr2ds.domainDf.residue_name.apply(
            lambda x: self._AATYPE[x])

        # backup info
        self._backup_all = self.pepr2ds.domainDf.copy()

        number_of_HP_per_pdbs = self.pepr2ds.domainDf.query('IBS == True and atom_name == "CB"').groupby("cathpdb").progress_apply(lambda group: self._count_hydr_protr_perpdb(group))
        nohydrprotrpdbs = list(number_of_HP_per_pdbs[number_of_HP_per_pdbs==0].index)
        self._backup_prot_HPIBS = self.pepr2ds.domainDf.query("cathpdb not in @nohydrprotrpdbs")
        self._backup_prot_noHPIBS = self.pepr2ds.domainDf.query("cathpdb in @nohydrprotrpdbs")
    
    """
    ### AUXILIARY METHODS
    """

    def _count_hydr_protr_perpdb(self, group):
        g = group.query("protrusion == True and atom_name == 'CB' and type.str.contains('Hydrophobic')",
                        engine="python")
        return (len(g))

    def _save_fig(self,
                  figname, 
                  folder="article", 
                  format="png", 
                  dpi=300, 
                  bbox_inches='tight', 
                  transparent=False, ):

        Path(f"{self.settings.FIGURESFOLDER}/{folder}").mkdir(parents=True, exist_ok=True)

        if format in ["tiff",'tif']:
            plt.savefig(f"{self.settings.FIGURESFOLDER}/{folder}/{figname}.{format}", dpi=dpi, bbox_inches=bbox_inches, transparent=transparent, pil_kwargs={"compression": "tiff_lzw"})
        else:
            plt.savefig(f"{self.settings.FIGURESFOLDER}/{folder}/{figname}.{format}", dpi=dpi, bbox_inches=bbox_inches, transparent=transparent)

    def _get_palette_OR(self, oddsratio):
        palette = []
        for index, row in oddsratio.iterrows():
            if row["pvalue"] < 0.05:
                if row.oddsratio > 0:
                    palette.append(self.dict_palette_IBS["IBS"])
                else:
                    palette.append(self.dict_palette_IBS["nonIBS"])
            else:
                palette.append("gray")
        return(palette)

    def _significance(self, ortable):
        def __add_symbol(row):
            if row.pvalue > 0.05:
                return "ns"
            elif row.pvalue <= 0.0001:
                return '****'
            elif row.pvalue <= 0.001:
                return '***'
            elif row.pvalue <= 0.01:
                return '**'
            elif row.pvalue <= 0.05:
                return '*'
            else:
                return 'NA'
            
        ortable["significance"] = ortable.apply(lambda x: __add_symbol(x), axis=1)
        ortable["labels"] ="("+ortable["significance"]+") " + ortable.iloc[:,0].astype(str)
        return ortable
    
    def _printmd(self, xstring, color=None):
        if self._silent_eq_test:
            return
        else:
            if self.settings.USING_NOTEBOOK:
                colorstr = "<span style='color:{}'>{}</span>".format(color, xstring)
                display(Markdown(colorstr))
            else:
                print(xstring)
                print("")

    def _equality_test(self,
                       pop1,                         # array with continuous or discrete data
                       pop2,                         # array with continuous or discrete data
                       ALTERNATIVE = "two-sided",    # 'two-sided', 'less' or 'greater'
                       pairwised = False,
                       silent = False,
                       returnPval = False):

        self._silent_eq_test = silent

        print("\n----------------------------------------")
        self._printmd("**STATISTICAL TEST BETWEEN TWO SAMPLES**")
        self._printmd(f" - ALTERNATIVE HYPOTHESIS = {ALTERNATIVE}")
        
        sign = {"two-sided":"≠",
           "less":"<",
           "greater":">"}

        self._printmd("**NORMALITY TEST (shapiro)**")
        normality = True
        self._printmd("*The two samples should follow a normal law to use a standard t.test*")

        normpop1 = stats.shapiro(pop1).pvalue
        normpop2 = stats.shapiro(pop2).pvalue
        if normpop1 < 0.05:
            self._printmd(f"---- Sample 1 shapiro test pvalue = {normpop1:.2E}, <= 0.05. This sample does NOT follow a normal law", color='red')
            normality = False
        else: 
            self._printmd(f"---- Sample 1 shapiro test pvalue = {normpop1:.2E}, > 0.05. This sample follows a normal law", color='blue')
        if normpop2 < 0.05:
            self._printmd(f"---- Sample 1 shapiro test pvalue = {normpop2:.2E}, <= 0.05. This sample does NOT follow a normal law", color='red')
            normality = False
        else: 
            self._printmd(f"---- Sample 1 shapiro test pvalue = {normpop2:.2E}, > 0.05. This sample follows a normal law", color='blue')

        if normality == True:
            self._printmd("Both samples follow a normal law")

            if pairwised == True:
                self._printmd("**TTest_REL Pairwise test **")
                equalstat, equalpval = stats.ttest_rel(pop1,pop2)
                if returnPval:
                    print("-- DONE\n\n")
                    return equalpval
            else: 
                print("Performing variance equality test")
                varstat, varpval = stats.levene(pop1,pop2)
                #Levene, pval < 0.05 --> Variances not equal.
                #Levene, pval > 0.05 --> Not significative and the hypothesis H0 is not rejected

                self._printmd("-- Null hypothesis : the variance of both samples are equal")
                print(f"---- Variance test --> stat={varstat:.2E}, p-value={varpval:.3E}")
                if varpval < 0.05:
                    self._printmd("P value <= 0.05, H0 rejected. The variances are not equal. Performing Welch’s t-test", color="red")
                    equal_var = False
                else:
                    self._printmd("Pvalue > 0.05, the variances are not equal. Performing a standard independent 2-sample test", color="blue")
                    equal_var = True
                equalstat, equalpval = stats.ttest_ind(pop1,
                                           pop2,
                                           equal_var=equal_var,)
                
            print(f"t-test --> stat={equalstat:.2E}, p-value={equalpval:.3E}")
            print(f"  Null hypothesis: the averages of both samples are equal")
            print(f"  Alternative hypothesis: average(sample2) {sign[ALTERNATIVE]} average(sample2)")
            if equalpval > 0.05:
                self._printmd("pvalue > 0.05, we cannot reject the null hypothesis of identical averages between both populations", color="blue")
            else:
                self._printmd("pvalue <= 0.05, the null hypothesis is rejected, the two samples are different", color="red")
        else:
            self._printmd("At least one sample does not follow a normal law")
            if pairwised==True:
                self._printmd("**WILCOXON SIGNED-RANK TEST**")
                self._printmd(f"  Null hypothesis: the two distributions are equal")
                self._printmd(f"  Alternative hypothesis: pop1 {sign[ALTERNATIVE]} pop2")
                stat, pval = stats.wilcoxon(pop1,pop2, alternative=ALTERNATIVE)
            else:
                self._printmd("Performing a Wilcoxon rank sum test with continuity correction")
                self._printmd("**WILCOXON RANK SUM TEST WITH CONTINUITY CORRECTION (or Mann-Whitney test)**")
                self._printmd(f"  Null hypothesis: the two distributions are equal")
                self._printmd(f"  Alternative hypothesis: pop1 {sign[ALTERNATIVE]} pop2")
                stat, pval = stats.mannwhitneyu(pop1,pop2, alternative=ALTERNATIVE)
            if pval < 0.05:
                self._printmd(f"  pvalue = {pval:.2E} which is <= 0.05. The null hypothesis is rejected. The alternative hypothesis (f{ALTERNATIVE}) is valid and the two distributions are statistically different", color="red")
            else:
                self._printmd(f"  pvalue = {pval:.2E} which is > 0.05. Null hypothesis not rejected. Both distributions are NOT statistically different", color="blue")
            if returnPval:
                print("-- DONE\n\n")
                return pval
        
        print("-- DONE\n")

    def _move_seaborn_legend(self, ax, new_loc, title=None, invert=True, order=None, **kws):
        # from https://github.com/mwaskom/seaborn/issues/2280#issuecomment-692350136
        old_legend = ax.legend_
        handles = old_legend.legendHandles
        labels = [t.get_text() for t in old_legend.get_texts()]

        labels_handles = dict(zip(labels,handles))

        if invert:
            handles=handles[::-1]
            labels=labels[::-1]
        if title == None:
            title = old_legend.get_title().get_text()

        if order:
            handles = [labels_handles[x] for x in order]
            labels = order

        ax.legend(handles, labels, loc=new_loc, title=title, **kws)

    """
    ### All methods below just encapsulate the steps in Notebook #4
    """

    #############################################################################
    # Figure 2 in the 2022 paper: composition of the exposed IBS for all proteins
    def make_figure_composition_of_exposed_IBS(self):
        dataCathIBS = self.pepr2ds.domainDf.query("IBS == True and exposed == True").drop_duplicates(['residue_name', 'residue_number', 'cathpdb']).residue_name.value_counts(normalize=True)*100
        dataCathnoIBS = self.pepr2ds.domainDf.query("IBS == False and exposed == True").drop_duplicates(['residue_name', 'residue_number', 'cathpdb']).residue_name.value_counts(normalize=True)*100

        plt.subplots(figsize=(14, 4), )
        gs = gridspec.GridSpec(ncols=2, nrows=2)
        gs.update(hspace=0.3,)
        ax0 = plt.subplot(gs[0,0])
        ax1 = plt.subplot(gs[1,0])
        ax2 = plt.subplot(gs[0,1])
        ax3 = plt.subplot(gs[1,1])
        self._composition_of_exposed_IBS(dataCathIBS, ax=ax0, PERTYPE=True, putLegend=False, xlabel="")
        self._composition_of_exposed_IBS(dataCathnoIBS, ax=ax1, PERTYPE=True, putLegend=True, )

        self._composition_of_exposed_IBS(dataCathIBS, ax=ax2, PERTYPE=False, putLegend=False, xlabel="")
        self._composition_of_exposed_IBS(dataCathnoIBS, ax=ax3, PERTYPE=False, putLegend=True)
        _ = ax0.text(-0.07,0.15, "IBS Exposed", transform=ax0.transAxes, fontsize=10, rotation=90,weight="bold")
        _ = ax1.text(-0.07,0, "non-IBS Exposed", transform=ax1.transAxes, fontsize=10, rotation=90, weight="bold")
        _ = ax1.text(0.5,1.1,"Physicochemical properties", transform=ax0.transAxes, fontsize=12,  ha="center", weight="bold")
        _ = ax1.text(0.5,1.1,"Amino acid", transform=ax2.transAxes, fontsize=12, ha="center", weight="bold")

        _= ax0.text(-0.1,1.02, "A",transform=ax0.transAxes, fontsize=20)
        _= ax2.text(-0.1,1.02, "B",transform=ax2.transAxes, fontsize=20)
        _= ax1.text(-0.1,1.02, "C",transform=ax1.transAxes, fontsize=20)
        _= ax3.text(-0.1,1.02, "D",transform=ax3.transAxes, fontsize=20)

        self._save_fig("Fig 2",format=self.FORMAT)

    def _composition_of_exposed_IBS(self, data, ax=None, PERTYPE=False, putLegend=True, xlabel="Composition (%)"):
        graph_data = data.to_frame()
        if PERTYPE:
            order_legend=["Positive",
                          "Negative",
                          "Polar",
                          "Non-polar",
                          "Hydrophobic,H-aromatic",
                          "Hydrophobic,H-non-aromatic",
                         ]
            color_palette=self.colorsPerType
            hue="type"
            weights="Percentage_Type"
            # graph_res_data=graph_res_data.drop_duplicates(["domain","type","Percentage_Type"])
        else:
            order_legend=["LYS","ARG",
                          "ASP","GLU",
                          "HIS","ASN","GLN","THR","SER",
                          "PRO","ALA","VAL","GLY",
                           "TYR","TRP","PHE",
                           "LEU","ILE","CYS","MET",
                         ]
            color_palette = {x:self.COLORS_taylor[x] for x in list(self.COLORS_taylor.keys())}
            hue='residue_name'
            weights='Percentage'
        
        graph_data.reset_index(inplace=True)
        graph_data.columns = ["residue_name","Percentage"]

        graph_data["type"] = graph_data["residue_name"].apply(lambda x: self._AATYPE[x])

        graph_data = graph_data.set_index(["type"])
        graph_data["Percentage_Type"] = graph_data.groupby("type").Percentage.sum()
        graph_data.reset_index(inplace=True)

        graph_data["data"] = "data"

        if PERTYPE:
            graph_data = graph_data.drop_duplicates("type")

        graph = sns.histplot(graph_data,
                             y="data",
                             hue=hue,
                             weights=weights,
                             multiple='stack',
                             hue_order=order_legend[::-1],
                             edgecolor='k',
                             linewidth=0.1,
                             palette=color_palette,
                             legend=putLegend,
                             ax=ax)
        graph.set(ylabel="", xlabel=xlabel)
        graph.set_yticklabels("")
        graph.set(xlim=[-1,101])

        if PERTYPE:
            for rec, label in zip(graph.patches,graph_data['Percentage_Type'].round(1).astype(str)):
                        height = rec.get_height()
                        width = rec.get_width()
                        val = f"{rec.get_width():.1f} "
                        size=12
                        """
                        if PERTYPE:
                            size = 8
                        else:
                            size = 4
                        """
                        ax.text( (rec.xy[0]+rec.get_width()/2), (rec.xy[1]+rec.get_height()/2), val, size=size, color="#383838", ha = 'center', va='center',)

        if putLegend==True:
            self._move_seaborn_legend(graph, 
                                      new_loc='center', 
                                      title="",
                                      order=order_legend,
                                      ncol=math.ceil(len(order_legend)/4), 
                                      bbox_to_anchor=(0.5,-0.7))


    ####################################################################################
    # Figure 3 in the 2022 paper: protrusions and hydrophobic protrusions in the dataset
    def make_figure_protrusions(self, show_equality_test=True, show_percentage_without_HP=True):
        count_protr_whole = self.pepr2ds.domainDf.query('atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_protrusion_per_pdb(group)).to_frame("count") 
        count_hydr_protr_whole = self.pepr2ds.domainDf.query('atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_hydr_protrusion_per_pdb(group)).to_frame("count") 
        # count_polar_protr_whole = self.pepr2ds.domainDf.query('atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_polar_protrusion_per_pdb(group)).to_frame("count") 
        # count_nonpolar_protr_whole = self.pepr2ds.domainDf.query('atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_nonpolar_protrusion_per_pdb(group)).to_frame("count") 
        # count_positive_protr_whole = self.pepr2ds.domainDf.query('atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_positive_protrusion_per_pdb(group)).to_frame("count") 
        # count_negative_protr_whole = self.pepr2ds.domainDf.query('atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_negative_protrusion_per_pdb(group)).to_frame("count") 

        count_protr_IBS = self.pepr2ds.domainDf.query('IBS == True and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_protrusion_per_pdb(group)).to_frame("count") 
        count_hydr_protr_IBS = self.pepr2ds.domainDf.query('IBS == True and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_hydr_protrusion_per_pdb(group)).to_frame("count") 
        # count_polar_protr_IBS = self.pepr2ds.domainDf.query('IBS == True and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_polar_protrusion_per_pdb(group)).to_frame("count") 
        # count_nonpolar_protr_IBS = self.pepr2ds.domainDf.query('IBS == True and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_nonpolar_protrusion_per_pdb(group)).to_frame("count") 
        # count_positive_protr_IBS = self.pepr2ds.domainDf.query('IBS == True and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_positive_protrusion_per_pdb(group)).to_frame("count") 
        # count_negative_protr_IBS = self.pepr2ds.domainDf.query('IBS == True and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_negative_protrusion_per_pdb(group)).to_frame("count") 

        count_protr_nonIBS = self.pepr2ds.domainDf.query('IBS == False and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_protrusion_per_pdb(group)).to_frame("count") 
        count_hydr_protr_nonIBS = self.pepr2ds.domainDf.query('IBS == False and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_hydr_protrusion_per_pdb(group)).to_frame("count") 
        # count_polar_protr_nonIBS = self.pepr2ds.domainDf.query('IBS == False and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_polar_protrusion_per_pdb(group)).to_frame("count") 
        # count_nonpolar_protr_nonIBS = self.pepr2ds.domainDf.query('IBS == False and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_nonpolar_protrusion_per_pdb(group)).to_frame("count") 
        # count_positive_protr_nonIBS = self.pepr2ds.domainDf.query('IBS == False and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_positive_protrusion_per_pdb(group)).to_frame("count") 
        # count_negative_protr_nonIBS = self.pepr2ds.domainDf.query('IBS == False and atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._count_negative_protrusion_per_pdb(group)).to_frame("count") 

        freqs_hydro_protrusion_whole = self.pepr2ds.domainDf.query('atom_name == "CB"').groupby('cathpdb').progress_apply(lambda group: self._calc_fract(group,"type.str.contains('Hydrophobic')", "protrusion == True")).to_frame("count") 
        freqs_hydro_protrusion_IBS = self.pepr2ds.domainDf.query('atom_name == "CB" and IBS == True').groupby('cathpdb').progress_apply(lambda group: self._calc_fract(group,"type.str.contains('Hydrophobic')", "protrusion == True")).to_frame("count") 
        freqs_hydro_protrusion_nonIBS = self.pepr2ds.domainDf.query('atom_name == "CB" and IBS == False').groupby('cathpdb').progress_apply(lambda group: self._calc_fract(group,"type.str.contains('Hydrophobic')", "protrusion == True")).to_frame("count")

        ###
        sns.set(font_scale=1.2)
        sns.set_style("whitegrid") #Seaborn style

        count_protr_whole["surface"] = "whole"
        count_hydr_protr_whole["surface"] = "whole"
        count_hydr_protr_IBS["surface"] = "IBS"
        count_protr_IBS["surface"] = "IBS"
        count_hydr_protr_nonIBS["surface"] = "nonIBS"
        count_protr_nonIBS["surface"] = "nonIBS"

        #WHOLE
        count_protr_whole["Type"] = "All"
        count_hydr_protr_whole["Type"] = "Hydrophobic"
        # count_positive_protr_whole["Type"] = "positive"
        # count_polar_protr_whole["Type"] = "polar"
        # count_nonpolar_protr_whole["Type"] = "nonpolar"
        # count_negative_protr_whole["Type"] = "negative"

        #IBS
        count_protr_IBS["Type"] = "All"
        count_hydr_protr_IBS["Type"] = "Hydrophobic"
        # count_positive_protr_IBS["Type"] = "positive"
        # count_polar_protr_IBS["Type"] = "polar"
        # count_nonpolar_protr_IBS["Type"] = "nonpolar"
        # count_negative_protr_IBS["Type"] = "negative"

        #non IBS
        count_protr_nonIBS["Type"] = "All"
        count_hydr_protr_nonIBS["Type"] = "Hydrophobic"
        # count_positive_protr_nonIBS["Type"] = "positive"
        # count_polar_protr_nonIBS["Type"] = "polar"
        # count_nonpolar_protr_nonIBS["Type"] = "nonpolar"
        # count_negative_protr_nonIBS["Type"] = "negative"

        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        
        if self.INCLUDE_AF_FROM_START:
            ylim=[0,900]
        else:
            ylim=[0,575]
        
        #GRAPH A
        count_whole = pd.concat([count_protr_whole,
                                 count_hydr_protr_whole,
                                #count_positive_protr_whole,
                                #count_polar_protr_whole,
                                #count_nonpolar_protr_whole,
                                #count_negative_protr_whole,
                                ], axis=0).reset_index()
        max_x_whole = count_whole["count"].max()

        graph_whole = sns.histplot(count_whole, 
                                   x="count",
                                   hue="Type",
                                   #stat="probability",
                                   #kind="kde",
                                   palette=["green","#ff7F0E"],
                                   alpha=0.5,
                                   bins=list(range(0,max_x_whole)),
                                   edgecolor="gray", linewidth=0.2,
                                   common_bins=False,
                                   common_norm=False,
                                   ax=axs[0],)

        _ = sns.despine(left=False, bottom=False, top=False, right=False) #All 4 borders

        graph_whole.set(xlabel="Number of protrusions",
                        ylabel="Number of structures",
                        #ylim=[0,0.42],
                        #ylim=[0,575]
                        ylim=ylim,
                        xlim=[0,48],)
        # _ = graph_whole.set_title("Whole surface", fontsize=11)

        #Graph B
        count_IBS = pd.concat([count_protr_IBS,
                               count_hydr_protr_IBS,
                               #count_positive_protr_IBS,
                               #count_polar_protr_IBS,
                               #count_nonpolar_protr_IBS,
                               #count_negative_protr_IBS,
                              ], axis=0).reset_index()
        max_x_IBS = count_IBS["count"].max()

        graph_IBS = sns.histplot(count_IBS, 
                                 x="count",
                                 hue="Type",
                                 #stat="probability",
                                 bins=list(range(0,max_x_IBS)),
                                 #kind="kde",
                                 alpha=0.5,
                                 palette=["green","#ff7F0E"],
                                 edgecolor="gray", linewidth=0.2,
                                 common_norm=False,
                                 ax=axs[1])

        _ = graph_IBS.set(xlabel="Number of protrusions",
                          ylabel="",
                          ylim=ylim,
                          #ylim=[0,575],
                          xlim=[0,48],)
        #_ = graph_IBS.set_title("IBS surface", fontsize=11)
        _ = sns.despine(left=False, bottom=False, top=False, right=False) #All 4 borders

        #Write letters
        freqs_hydro_protrusion_whole["Surface"] = "whole"
        freqs_hydro_protrusion_nonIBS["Surface"] = "nonIBS"
        freqs_hydro_protrusion_IBS["Surface"] = "IBS"
        ratio_graph_data = pd.concat([freqs_hydro_protrusion_nonIBS,freqs_hydro_protrusion_IBS], axis=0).reset_index()
        max_x = ratio_graph_data["count"].max()

        #graph1 = sns.histplot(test.query("Surface == 'whole'"), x="count", bins=np.arange(0,0.8,0.05), alpha=0.6,color="#8de5A1", stat="probability",  ax=axs[2], linewidth=0.2,edgecolor="gray")
        #graph2 = sns.histplot(test.query("Surface == 'IBS'"), x="count", bins=np.arange(0,0.8,0.05), alpha=0.6, color="#cb6679",stat="probability",ax=axs[2], linewidth=0.2,edgecolor="gray")
        graph_ratio = sns.histplot(ratio_graph_data, 
                                  hue="Surface", 
                                  x="count", 
                                  bins=np.arange(0,0.8,0.05), 
                                  alpha=0.6,
                                  #stat="probability", 
                                  palette=["#cb6679","#69aacf"],  
                                  ax=axs[2],
                                  linewidth=0.2,
                                  kde=True,
                                  common_norm=False,
                                  edgecolor="gray")
        #old color = #cb6679 / #69aacf
        _ = sns.despine(left=False, bottom=False, top=False, right=False) #All 4 borders

        #graph.set(xlabel="Number of protrusions", ylabel="Number of protein",title="Number of hydrophobic protrusions per protein (Whole)")

        # ---- Legend ----
        #top_bar = mpatches.Patch(color='#8de5a1', label='Whole surface', alpha=0.6, linewidth=0)
        #bottom_bar = mpatches.Patch(color='#cb6679', label='IBS surface', alpha=0.6, linewidth=0)
        #axs[2].legend(handles=[top_bar, bottom_bar])

        _ = axs[2].set(xlabel="Ratio of protrusions being hydropobic",
                       ylabel="",
                       xlim=[0,0.7],
                       #ylim=[0,575],
                       ylim=ylim,)
        #_ = axs[2].text(-0.1,1.1, "Ratio of protrusions being hydrophobic", transform=axs[2].transAxes, fontsize=15)
        _ = axs[2].tick_params(labelsize=10)

        _ = axs[0].text(-0.1,1.02, "A", transform=axs[0].transAxes, fontsize=20)
        #_ = axs[0].text(0.5,1.1, "Number of protrusions per structure", transform=axs[0].transAxes, fontsize=15)
        _ = axs[1].text(-0.1,1.02, "B", transform=axs[1].transAxes, fontsize=20)
        _ = axs[2].text(-0.1,1.02, "C", transform=axs[2].transAxes, fontsize=20)

        self._save_fig(f"Fig 3{self.FILESUFFIX}", transparent=False, format=self.FORMAT)

        if show_equality_test:
            pop1 = ratio_graph_data.query("Surface == 'IBS'")["count"]
            pop2 = ratio_graph_data.query("Surface == 'nonIBS'")["count"]
            self._equality_test(pop1,pop2, "greater")

        if show_percentage_without_HP:
            print("----------------------------------------")
            print("**PERCENTAGE OF THE DATASET WITH/WITHOUT HP**")
            print("-- nonIBS SURFACE --")
            print(f" protrusions .............. {count_protr_nonIBS['count'].mean():.2f} ± {count_protr_nonIBS['count'].std():.2f}") 
            print(f" hydrophobic protrusions .. {count_hydr_protr_nonIBS['count'].mean():.2f} ± {count_hydr_protr_nonIBS['count'].std():.2f}") 

            print("-- IBS SURFACE --")
            print(f" protrusions .............. {count_protr_IBS['count'].mean():.2f} ± {count_protr_IBS['count'].std():.2f}") 
            print(f" hydrophobic protrusions .. {count_hydr_protr_IBS['count'].mean():.2f} ± {count_hydr_protr_IBS['count'].std():.2f}") 

            print("-- whole surface --")
            print(f" protrusions .............. {count_protr_whole['count'].mean():.2f} ± {count_protr_whole['count'].std():.2f}") 
            print(f" hydrophobic protrusions .. {count_hydr_protr_whole['count'].mean():.2f} ± {count_hydr_protr_whole['count'].std():.2f}")

            ###
            print("\nPercentage of the dataset without hydrophobic protrusions")
            perc = ratio_graph_data.groupby("Surface").apply(lambda x: len(x.query("count <0.05")) / len(x) * 100)
            vals = ratio_graph_data.groupby("Surface").apply(lambda x: len(x.query("count <0.05")))
            print("      percentage")
            print(perc)
            print("      values")
            print(vals)

            ratio_graph_data.groupby("Surface").apply(lambda x: len(x.query("count >=0.05")))

            ###
            res = self.pepr2ds.domainDf.groupby("cathpdb", as_index=False).progress_apply(lambda x: self._count_percentage_of_the_convexhull(x))
            print(f"Percentage of the convexhull being part of the IBS: {res.loc[(slice(None), 'fraction_IBS'), :].mean()[0]:.2f}% ± {res.loc[(slice(None), 'fraction_IBS'), :].std()[0]:.2f}%")
            print(f"Percentage of the convexhull NOT being part of the IBS: {res.loc[(slice(None), 'fraction_nonIBS'), :].mean()[0]:.2f}% ± {res.loc[(slice(None), 'fraction_nonIBS'), :].std()[0]:.2f}%")

            ###
            res = self._backup_prot_noHPIBS.groupby("cathpdb", as_index=False).progress_apply(lambda x: self._count_fraction_of_IBS(x))
            print(f"Percentage of the IBS covered by protrusion for protreins without hydrophobic protrusions at their IBS: {res.iloc[:,1].mean():.2f}±{res.iloc[:,1].std():.2f}%")
            print("\n-- DONE\n\n")

    def _count_protrusion_per_pdb(self, group):
        N = group.query("protrusion == True")
        return len(N)

    def _count_hydr_protrusion_per_pdb(self, group):
        N = group.query("protrusion == True and type.str.contains('Hydrophobic')", engine='python')
        return len(N)

    def _count_polar_protrusion_per_pdb(self, group):
        N = group.query("protrusion == True and type.str.contains('Polar')", engine='python')
        return len(N)

    def _count_nonpolar_protrusion_per_pdb(self, group):
        N = group.query("protrusion == True and type.str.contains('Non-polar')", engine='python')
        return len(N)

    def _count_negative_protrusion_per_pdb(self, group):
        N = group.query("protrusion == True and type.str.contains('Negative')", engine='python')
        return len(N)

    def _count_positive_protrusion_per_pdb(self, group):
        N = group.query("protrusion == True and type.str.contains('Positive')", engine='python')
        return len(N)

    def _calc_fract(self, group, property1, property2):
        # determine the fraction to get property1 IN property2
        total = len(group.query(property2, engine='python'))
        pr1 = len(group.query(f"{property1} and {property2}", engine='python'))
        if total == 0:
            return 0
        return pr1/total

    def _count_percentage_of_the_convexhull(self, pdb):
        # evaluate the percentage of the convexhull to be part of the IBS and not IBS
        df = pdb.query("atom_name == 'CB' and convhull_vertex == True")
        n_vertices = len(df)
        n_vertices_IBS = len(df.query("IBS == True"))
        n_vertices_nonIBS = len(df.query("IBS == False"))

        ret = [n_vertices_IBS/n_vertices*100,n_vertices_nonIBS/n_vertices*100]

        return(pd.DataFrame(ret, index=["fraction_IBS","fraction_nonIBS"]))
    
    def _count_fraction_of_IBS(self, pdb):
        #evaluate the percentage of the convexhull to be part of the IBS and not IBS.
        df = pdb.query("atom_name == 'CB' and convhull_vertex == True")
        n_vertices = len(df)
        n_vertices_IBS = len(df.query("IBS == True"))
        n_protrusion = len(df.query("protrusion == True"))

        ret = [n_vertices_IBS/n_vertices*100]

        return(ret)

    ######################################################################################
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

