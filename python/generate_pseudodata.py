#! /usr/bin/env python
import csv
import os
from mATLASplotlib import canvases
import numpy as np
import random

# datasets = ["bkg_only", "bkg_only", "bkg_only", "m260:5000", "m260:2100", "m300:2000", "m350:1500", "m500:600", "m750:100", "m1000:40"]
datasets = ["m300:5000",  "m750:5000", "bkg_only", "bkg_only", "bkg_only", "m260:350", "m275:1550", "m300:1200", "m325:700", "m350:700", "m400:400", "m450:600", "m500:350", "m750:150", "m1000:40"]
np.random.seed(20170711)

for idx, dataset in enumerate(np.random.permutation(datasets), start=1):
    print "** Generating sample {}: {} **".format(idx, dataset)
    for mass_category in ["high", "low"]:
        for tag_category in ["0", "1", "2"]:
            f_output_name = os.path.join("output", "pseudodata", "pseudodata_sample_{}_{}Mass_{}tag.txt".format(idx, mass_category, tag_category))
            if os.path.isfile(f_output_name):
                print "File {} exists! Skipping".format(f_output_name)
                continue
            with open(f_output_name, "wb") as f_output:
                nBkg, nSig = 0, 0
                # Background
                with open(os.path.join("input", "m_yyjj_SM_bkg_{}Mass_{}tag_tightIsolated.csv".format(mass_category, tag_category)), "rb") as f_input:
                    for row in csv.reader(f_input, delimiter="\t"):
                        # if random.random() < float(row[1]):
                        if np.random.random() < float(row[1]):
                            f_output.write(row[0] + "\n")
                            nBkg += 1
                # Signal
                if ":" in dataset:
                    mass, xs_fb = dataset.split(":")
                    weight_scale = float(xs_fb) / 5000.
                    with open(os.path.join("input", "m_yyjj_Xhh_{}_{}Mass_{}tag_tightIsolated.csv".format(mass, mass_category, tag_category)), "rb") as f_input:
                        for row in csv.reader(f_input, delimiter="\t"):
                            # if random.random() < (float(row[1]) * weight_scale):
                            if np.random.random() < (float(row[1]) * weight_scale):
                                f_output.write(row[0] + "\n")
                                nSig += 1
                print "Generated {} mass {}-tag pseudodata with {} background events and {} signal events".format(mass_category, tag_category, nBkg, nSig)

# # Plot distributions
# for idx, _ in enumerate(datasets, start=1):
#     for mass_category in ["high", "low"]:
#         x_range = {"high":(335, 1140), "low":(245, 485)}[mass_category]
#         for tag_category in ["0", "1", "2"]:
#             data = np.loadtxt(os.path.join("output", "pseudodata", "pseudodata_sample_{}_{}Mass_{}tag.txt".format(idx, mass_category, tag_category)))
#             hist, _edges = np.histogram(data, "auto", range=x_range)
#             # hist, _edges = np.histogram(data, int(data.size / 10 + 1), range=x_range)
#             centres = [0.5 * (e1 + e2) for e1, e2 in zip(_edges[:-1], _edges[1:])]
#             canvas = canvases.Simple()
#             canvas.plot_dataset((centres, None, hist, np.sqrt(hist)), style="scatter yerror", colour="black")
#             canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
#             canvas.set_axis_label("y", "Events / bin")
#             canvas.set_axis_range("x", x_range)
#             canvas.save_to_file(os.path.join("plots", "pseudodata_sample_{}_{}Mass_{}tag".format(idx, mass_category, tag_category)))

