#! /usr/bin/env python
from collections import defaultdict
import csv
import os

input_directory = "/afs/cern.ch/user/j/jrobinso/HGamma/PyAnalysis/text/"

scale_factors = defaultdict(dict)
with open(os.path.join(input_directory, "scale_factors.csv"), "rb") as f_input:
    for row in csv.reader(f_input, delimiter="\t"):
        scale_factors[row[0]][row[1]] = row[2]

for mass_category in ["low", "high"]:
    for tag_category in [0, 1, 2]:
        with open(os.path.join("input", "m_yyjj_SM_bkg_{}Mass_{}tag_tightIsolated.csv".format(mass_category, tag_category)), "wb") as f_output:
            # for bkg in ["SM_bkg_yy", "SM_bkg_yj", "SM_bkg_jy", "SM_bkg_jj", "bbH", "ggH", "ttH", "VBFH", "WH", "ZH"]:
            for bkg in ["SM_bkg_yy", "SM_bkg_yj", "SM_bkg_jy", "SM_bkg_jj", "ggH", "ttH", "VBFH", "WH", "ZH"]:
                scale_factor = float(scale_factors["m_yyjj_{}Mass_{}tag_tightIsolated".format(mass_category, tag_category)][bkg])
                with open(os.path.join(input_directory, "m_yyjj_{}_{}Mass_{}tag_tightIsolated.csv".format(bkg, mass_category, tag_category)), "rb") as f_input:
                    for row in csv.reader(f_input, delimiter="\t"):
                        f_output.write("\t".join((row[0], str(float(row[1]) * scale_factor))) + "\n")
