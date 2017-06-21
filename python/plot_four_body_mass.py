#! /usr/bin/env python
import csv
import math
import os
from collections import defaultdict, OrderedDict

from mATLASplotlib import canvases

colours = OrderedDict([("novosibirsk", "purple"), ("modified_gamma", "green"), ("modified_landau", "blue")])
base_path = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))


for mass_category in ["high", "low"]:
    for tag_category in ["0", "1", "2"]:
        # Check that input file exists
        input_path = os.path.join(base_path, "output", "spurious_signal_{}Mass_{}tag.csv".format(mass_category, tag_category))
        if not os.path.exists(input_path):
            continue

        n_spurious = defaultdict(list)
        Z_spurious = defaultdict(list)
        mass_points = []

        with open(input_path, "rb") as f_input:
            for row in csv.reader(f_input, delimiter=" "):
                function, mass, n_sig, Z, Z_uncertainty, chi2, ndof = row
                if int(mass) > 0:
                    n_spurious[function].append((int(mass), float(n_sig)))
                    Z_spurious[function].append((int(mass), float(Z), float(Z_uncertainty)))
                    mass_points.append(int(mass))
        x_min, x_max = min(mass_points), max(mass_points)

        # Ensure that output directory exists
        if not os.path.exists(os.path.join(base_path, "output", "plots")):
            os.makedirs(os.path.join(base_path, "output", "plots"))

        # Plot n_spurious
        canvas = canvases.Simple()
        for function in colours.keys():
            x, y = zip(*n_spurious[function])
            canvas.plot_dataset((x, y), style="line join centres", label=function.replace("_", " ").title(), colour=colours[function])
        canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
        canvas.set_axis_label("y", "$N_{spurious}$")
        canvas.set_axis_range("x", (x_min, x_max))
        canvas.add_legend(0.97, 0.96, fontsize="medium", anchor_to="upper right")
        canvas.save_to_file(os.path.join(base_path, "output", "plots", "n_spurious_{}Mass_{}tag".format(mass_category, tag_category)))

        # Plot Z_spurious
        canvas = canvases.Simple()
        for function in colours.keys():
            x, y, y_err = zip(*Z_spurious[function])
            canvas.plot_dataset((x, None, y, y_err), style="binned band", colour=colours[function], alpha=0.2)
            canvas.plot_dataset((x, y), style="line join centres", label=function.replace("_", " ").title(), colour=colours[function])
        canvas.plot_dataset(([x_min, x_max], [0.2, 0.2]), style="line join centres", colour="red", linestyle="dashed")
        canvas.plot_dataset(([x_min, x_max], [-0.2, -0.2]), style="line join centres", colour="red", linestyle="dashed")
        canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
        canvas.set_axis_label("y", "$S / \Delta S$")
        canvas.set_axis_range("x", (x_min, x_max))
        canvas.set_axis_range("y", (-0.5, 0.6))
        canvas.add_legend(0.97, 0.96, fontsize="medium", anchor_to="upper right")
        canvas.save_to_file(os.path.join(base_path, "output", "plots", "Z_spurious_{}Mass_{}tag".format(mass_category, tag_category)))
