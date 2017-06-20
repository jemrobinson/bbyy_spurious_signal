#! /usr/bin/env python
import csv
import math
import os
from collections import defaultdict

from mATLASplotlib import canvases

colours = {"novosibirsk": "purple", "modified_gamma": "green", "modified_landau": "blue"}
base_path = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))


for mass_category in ["low", "high"]:
    # Check that input file exists
    input_path = os.path.join(base_path, "output", "spurious_signal_{}.csv".format(mass_category))
    if not os.path.exists(input_path):
        continue

    n_spurious = defaultdict(list)
    Z_spurious = defaultdict(list)

    with open(input_path, "rb") as f_input:
        for row in csv.reader(f_input, delimiter=" "):
            # function, mass, n_sig, n_sig_err, n_bkg, n_bkg_err, fit_result, Z = row
            function, mass, n_sig, Z, Z_uncertainty, chi2, ndof = row
            if int(mass) > 0:
                n_spurious[function].append((int(mass), float(n_sig)))
                Z_spurious[function].append((int(mass), float(Z), float(Z_uncertainty)))
            else:
                print row, "did not converge!"

    # Ensure that output directory exists
    if not os.path.exists(os.path.join(base_path, "output", "plots")):
        os.makedirs(os.path.join(base_path, "output", "plots"))

    # Plot n_spurious
    canvas = canvases.Simple()
    for function in n_spurious.keys():
        x, y = zip(*n_spurious[function])
        print x
        canvas.plot_dataset((x, y), style="line join centres", label=function.replace("_", " ").title(), colour=colours[function])
    canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
    canvas.set_axis_label("y", "$N_{spurious}$")
    canvas.add_legend(0.97, 0.96, fontsize="medium", anchor_to="upper right")
    canvas.save_to_file(os.path.join(base_path, "output", "plots", "n_spurious_{}".format(mass_category)))

    # Plot Z_spurious
    canvas = canvases.Simple()
    for function in Z_spurious.keys():
        print x
        x, y, y_err= zip(*Z_spurious[function])
        canvas.plot_dataset((x, None, y, y_err), style="binned band", colour=colours[function], alpha=0.4)
        canvas.plot_dataset((x, y), style="line join centres", label=function.replace("_", " ").title(), colour=colours[function])
    canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
    canvas.set_axis_label("y", "$S / \Delta S$")
    canvas.add_legend(0.97, 0.96, fontsize="medium", anchor_to="upper right")
    canvas.save_to_file(os.path.join(base_path, "output", "plots", "Z_spurious_{}".format(mass_category)))
