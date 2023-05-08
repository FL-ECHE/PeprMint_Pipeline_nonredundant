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


""" Small test comparing CATH and AF for the PH superfamily

This example was previously on the tag_IBS_domains_withAlphaFold notebook under 
"tools".

__author__ = ["Thibault Tubiana", "Phillippe Samer"]
__organization__ = "Computational Biology Unit, Universitetet i Bergen"
__copyright__ = "Copyright (c) 2022 Reuter Group"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Phillippe Samer"
__email__ = "samer@uib.no"
__status__ = "Prototype"
"""

# NB! Modifying sys.path is not a recommended practice in general! Exceptionally
# doing it here to make examples without modifing the repository structure
import sys
sys.path.append('..')

from src.settings import Settings
from src.dataset_manager import DatasetManager
from src.ibs_tagging import IBSTagging

from IPython.display import Markdown, display
import random
import numpy

class ExampleComparingAFCATH():

    def __init__(self):
        seeds = [1092593, 1984337]
        random.seed(seeds[0])
        numpy.random.seed(seeds[1])

        self.USING_NOTEBOOK = self._is_notebook()

    def main(self):
        self._setup()
        #self._make_reports()
        self._compare()

    def _setup(self):
        self.settings = Settings("./peprmint_default.config")

        self.cath_manager = DatasetManager(self.settings)
        self.cath_manager.load_light_dataset()
        self.af_manager = DatasetManager(self.settings)
        self.af_manager.load_light_dataset()

        self.cath_manager.add_IBS_data(db="cath")
        self.af_manager.add_IBS_data(db="alphafold")

    def _make_reports(self):
        self.cath_manager.make_IBS_analysis_report()
        self.af_manager.make_IBS_analysis_report()

    def _compare(self):
        protrusions_per_structure_CATH = self.cath_manager.get_protusion_count_after_IBS()
        protrusions_per_structure_AF = self.af_manager.get_protusion_count_after_IBS()
        uniprot_in_common = list(set(protrusions_per_structure_CATH.index).intersection(protrusions_per_structure_AF.index))
        self._equality_test(protrusions_per_structure_CATH[uniprot_in_common],
                            protrusions_per_structure_AF[uniprot_in_common],
                            pairwised=True)

    def test_equality_test1(self):
        n = 100
        sample1 = random.sample(range(1, n+1), n)               # all different
        sample2 = [random.randint(1, n) for u in range(n)]
        self._equality_test(sample1, sample2)

    def test_equality_test2(self):
        n = 1000
        sample1 = random.sample(range(1, n+1), n)               # all different
        sample2 = [random.gauss(n^2, n/4) for u in range(n)]
        self._equality_test(sample1, sample2)
    
    def _is_notebook(self) -> bool:
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell':
                return True   # jupyter notebook
            elif shell == 'TerminalInteractiveShell':
                return False  # terminal running IPython
            else:
                return False  # other type (?)
        except NameError:
            return False      # probably standard Python interpreter

    def _printmd(self, xstring, color=None):
        if self.USING_NOTEBOOK:
            colorstr = "<span style='color:{}'>{}</span>".format(color, xstring)
            display(Markdown(colorstr))
        else:
            print(xstring)
            print("")

    def _equality_test(self,
                       pop1,                         # array with continuous or discrete data
                       pop2,                         # array with continuous or discrete data
                       ALTERNATIVE = "two-sided",    # 'two-sided', 'less' or 'greater'
                       pairwised = False):
        
        print("----------------------------------------")
        self._printmd("**STATISTICAL TEST BETWEEN TWO SAMPLES**")
        
        self._printmd(f" - ALTERNATIVE HYPOTHESIS = {ALTERNATIVE}")
        from decimal import Decimal
        import scipy.stats as stats
        
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
        
        print("-- DONE\n\n")


if __name__ == "__main__":
    ex = ExampleComparingAFCATH()
    ex.main()
    #ex.test_equality_test1()
    #ex.test_equality_test2()
