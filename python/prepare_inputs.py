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
            with open(os.path.join("input", "m_yyjj_SM_bkg_{}Mass_{}tag_tightIsolated_positive_weights.csv".format(mass_category, tag_category)), "wb") as f_output_positive_weights:
                nTotal, sumWTotal = 0, 0.0
                for bkg in ["SM_bkg_yy", "SM_bkg_yj", "SM_bkg_jy", "SM_bkg_jj", "ggH", "ttH", "VBFH", "WH", "ZH", "bbH"]:
                    _masses, _weights = [], []
                    with open(os.path.join(input_directory, "m_yyjj_{}_{}Mass_{}tag_tightIsolated.csv".format(bkg, mass_category, tag_category)), "rb") as f_input:
                        for row in csv.reader(f_input, delimiter="\t"):
                            _masses.append(row[0])
                            _weights.append(float(row[1]))
                    nBkg, sum_weights = len(_weights), sum(_weights)
                    sum_positive_weights = sum([_w for _w in _weights if _w > 0.])
                    _scale_factor = sum_weights / sum_positive_weights if sum_weights > 0.0 else 0.0
                    for _mass, _weight in zip(_masses, _weights):
                        f_output.write("\t".join([_mass, "{0:.20f}".format(_weight)]) + "\n")
                        if _weight > 0.0:
                            f_output_positive_weights.write("\t".join([_mass, "{0:.20f}".format(_weight * _scale_factor)]) + "\n")
                    print ".... {} had {} events, equivalent to {:.3f} in data ({:.3f} if negative weights are ignored): scale factor of {:.3f} was applied".format(bkg, nBkg, sum_weights, sum_positive_weights, _scale_factor)
                    nTotal += nBkg
                    sumWTotal += sum_weights
        print "=> SM bkg summary: {} events, equivalent to {} in data".format(nTotal, sumWTotal)

        # BSM Higgs
        print "... and BSM Higgs"
        for csv_file in sorted(glob.glob(os.path.join(input_directory, "m_yyjj_Xhh_m*{}Mass*{}tag_tightIsolated.csv".format(mass_category, tag_category)))):
            mX = os.path.basename(csv_file).split("_")[3]
            with open(csv_file, "rb") as f_input:
                _masses, _weights = [], []
                for row in csv.reader(f_input, delimiter="\t"):
                    _masses.append(row[0])
                    _weights.append(float(row[1]))
            sum_weights = sum(_weights)
            sum_positive_weights = sum([_w for _w in _weights if _w > 0.])
            _scale_factor = sum_weights / sum_positive_weights if sum_weights > 0.0 else 0.0
            nSig = len(_weights)
            # Write outputs
            with open(os.path.join("input", os.path.basename(csv_file)), "wb") as f_output:
                with open(os.path.join("input", os.path.basename(csv_file)).replace(".csv", "_positive_weights.csv"), "wb") as f_output_positive_weights:
                    for _mass, _weight in zip(_masses, _weights):
                        f_output.write("\t".join([_mass, "{0:.20f}".format(_weight)]) + "\n")
                        if _weight > 0.0:
                            f_output_positive_weights.write("\t".join([_mass, "{0:.20f}".format(_weight * _scale_factor)]) + "\n")
            print ".... {} had {} events, equivalent to {:.3f} in data ({:.3f} if negative weights are ignored): scale factor of {:.3f} was applied".format(mX, nSig, sum_weights, sum_positive_weights, _scale_factor)
