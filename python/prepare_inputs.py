#! /usr/bin/env python
from collections import defaultdict
import csv
import os
import glob
import shutil

input_directory = "/afs/cern.ch/user/j/jrobinso/HGamma/PyAnalysis/text/blinded"

# Background samples
for mass_category in ["low", "high"]:
    for tag_category in [0, 1, 2]:
        print "Now considering {} mass, {}-tag category...".format(mass_category, tag_category)
        with open(os.path.join("input", "m_yyjj_SM_bkg_{}Mass_{}tag_tightIsolated.csv".format(mass_category, tag_category)), "wb") as f_output:
            nTotal, sumWTotal = 0, 0.0
            for bkg in ["SM_bkg_yy", "SM_bkg_yj", "SM_bkg_jy", "SM_bkg_jj", "ggH", "ttH", "VBFH", "WH", "ZH", "bbH"]:
                nBkg, sumWBkg = 0, 0.0
                with open(os.path.join(input_directory, "m_yyjj_{}_{}Mass_{}tag_tightIsolated.csv".format(bkg, mass_category, tag_category)), "rb") as f_input:
                    for row in csv.reader(f_input, delimiter="\t"):
                        output_weight = float(row[1])
                        f_output.write("\t".join((row[0], "{0:.20f}".format(output_weight))) + "\n")
                        nBkg += 1
                        sumWBkg += output_weight
                print ".... {} had {} events, equivalent to {} in data".format(bkg, nBkg, sumWBkg)
                nTotal += nBkg
                sumWTotal += sumWBkg
        print "=> Summary: {} events, equivalent to {} in data".format(nTotal, sumWTotal)

# BSM Higgs
for csv_file in glob.glob(os.path.join(input_directory, "m_yyjj_Xhh_m*.csv")):
    with open(os.path.join("input", os.path.basename(csv_file)), "wb") as f_output:
        with open(csv_file, "rb") as f_input:
            for row in csv.reader(f_input, delimiter="\t"):
                f_output.write("\t".join((row[0], "{0:.20f}".format(float(row[1])))) + "\n")
