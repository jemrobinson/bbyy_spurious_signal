#! /usr/bin/env python
from mATLASplotlib import canvases
from collections import defaultdict
import csv
import os
import math

colours = {"novosibirsk": "purple", "modified_gamma": "green", "modified_cauchy": "blue"}
base_path = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))


for mass_category in ["low", "high"]:
    n_spurious = defaultdict(list)
    Z_spurious = defaultdict(list)

    with open(os.path.join(base_path, "output", "spurious_signal_{}.csv".format(mass_category)), "rb") as f_input:
        for row in csv.reader(f_input, delimiter=" "):
            function, mass, n_sig, n_sig_err, n_bkg, n_bkg_err, fit_result, Z = row
            if True:
            # if fit_result == "0":
                n_spurious[function].append((int(mass), float(n_sig)))
                Z_spurious[function].append((int(mass), float(Z)))
                if float(Z) != 0:
                    print "n_sig_err = {}, sqrt(n_bkg) = {}, n_sig_err[no sum(w^2)] = {}".format(n_sig_err, math.sqrt(float(n_bkg)), float(n_sig) / float(Z))
            else:
                print row, "did not converge!"

    # Plot n_spurious
    canvas = canvases.Simple()
    for function in n_spurious.keys():
        x, y = zip(*n_spurious[function])
        canvas.plot_dataset((x, y), style="line join centres", label=function, colour=colours[function])
    canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
    canvas.set_axis_label("y", "$N_{spurious}$")
    canvas.add_legend(0.97, 0.96, fontsize="medium", anchor_to="upper right")
    canvas.save_to_file(os.path.join(base_path, "output", "n_spurious_{}".format(mass_category)))

    # Plot Z_spurious
    canvas = canvases.Simple()
    for function in Z_spurious.keys():
        x, y = zip(*Z_spurious[function])
        canvas.plot_dataset((x, y), style="line join centres", label=function, colour=colours[function])
    canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
    canvas.set_axis_label("y", "$S / \Delta S$")
    canvas.add_legend(0.97, 0.96, fontsize="medium", anchor_to="upper right")
    canvas.save_to_file(os.path.join(base_path, "output", "Z_spurious_{}".format(mass_category)))
