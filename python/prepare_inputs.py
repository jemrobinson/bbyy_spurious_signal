#! /usr/bin/env python
from collections import defaultdict
import csv
import os
import glob
import shutil

input_directory = "/afs/cern.ch/user/j/jrobinso/HGamma/PyAnalysis/text/"

# Sherpa + fakes + SM Higgs
scale_factors = defaultdict(dict)
with open(os.path.join(input_directory, "scale_factors.csv"), "rb") as f_input:
    for row in csv.reader(f_input, delimiter="\t"):
        scale_factors[row[0]][row[1]] = row[2]


mean_bkg_weight = {("low", 0): 0.0042164922200, ("low", 1): 0.00090105694345, ("low", 2): 0.00040002332151, ("high", 0): 0.0033953944451, ("high", 1): 0.00035164075374, ("high", 2): 0.000285774005317}

# Background samples
for mass_category in ["low", "high"]:
    for tag_category in [0, 1, 2]:
        weight_threshold = 0.0 #0.001 * mean_bkg_weight[(mass_category, tag_category)]
        with open(os.path.join("input", "m_yyjj_SM_bkg_{}Mass_{}tag_tightIsolated.csv".format(mass_category, tag_category)), "wb") as f_output:
            nAccepted, nRejected = 0, 0
            for bkg in ["SM_bkg_yy", "SM_bkg_yj", "SM_bkg_jy", "SM_bkg_jj", "ggH", "ttH", "VBFH", "WH", "ZH", "bbH"]:
                scale_factor = float(scale_factors["m_yyjj_{}Mass_{}tag_tightIsolated".format(mass_category, tag_category)][bkg])
                with open(os.path.join(input_directory, "m_yyjj_{}_{}Mass_{}tag_tightIsolated.csv".format(bkg, mass_category, tag_category)), "rb") as f_input:
                    for row in csv.reader(f_input, delimiter="\t"):
                        output_weight = float(row[1]) * scale_factor
                        if abs(output_weight) > weight_threshold:
                            f_output.write("\t".join((row[0], "{0:.20f}".format(output_weight))) + "\n")
                            nAccepted += 1
                        else:
                            nRejected += 1
            print "{} mass, {}-tag: accepted {}/{} [{}%] of events ".format(mass_category, tag_category, nAccepted, nAccepted + nRejected, (100 * nAccepted) / (nAccepted + nRejected))

# BSM Higgs
for csv_file in glob.glob(os.path.join(input_directory, "m_yyjj_Xhh_m*.csv")):
    with open(os.path.join("input", os.path.basename(csv_file)), "wb") as f_output:
        with open(csv_file, "rb") as f_input:
            for row in csv.reader(f_input, delimiter="\t"):
                f_output.write("\t".join((row[0], "{0:.20f}".format(float(row[1])))) + "\n")
