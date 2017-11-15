#! /usr/bin/env python
import csv
import os
from mATLASplotlib import canvases

base_path = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
# colours = ["#66a61e", "#e6ab02", "#1b9e77", "#d95f02", "#7570b3", "#e7298a"]
colours = ["#1b9e77", "#d95f02", "#7570b3"]
signal_source = "Expected_bbyy"

input_rows = []
with open(os.path.join(base_path, "output", "csv", "signal_bias.csv"), "rb") as f_input:
    for row in csv.reader(f_input, delimiter=" "):
        input_rows.append(row)

for mass_category in ["low", "high"]:
    for tag_category in ["0", "1", "2"]:
        category_rows = [[_x if "_" in _x else float(_x) for _x in row[2:]] for row in input_rows if row[0] == mass_category and row[1] == tag_category]
        if len(category_rows) == 0: continue
        # Initialise bias plot canvas
        canvas = canvases.Simple()
        maximum_bias = 0
        # Restrict to only rows with this signal source
        rows = [row for row in category_rows if row[1] == signal_source]
        mass_points = [row[0] for row in rows]
        # Get nEvents
        nEvtBiasAsimov = [100 * (row[3] - row[2]) / row[2] for row in rows]
        nEvtBiasAsimovErr = [100 * row[4] / row[2] for row in rows] # f = (S-T) / T => f_err = s_err / T
        nEvtBiasToys = [100 * (row[5] - row[2]) / row[2] for row in rows]
        nEvtBiasToysErr = [100 * row[6] / row[2] for row in rows]
        nToys = int(rows[0][7])
        # Construct boundaries
        maximum_bias = sorted(nEvtBiasToys, key=lambda x: abs(x))[-1]
        x_min, x_max = min(mass_points), max(mass_points)
        # Add Asimov points
        canvas.plot_dataset(mass_points, nEvtBiasAsimov, style="line join centres", colour=colours[0], label="Asimov bkg + signal MC")
        # Add toys points
        canvas.plot_dataset(mass_points, nEvtBiasToys, style="line join centres", colour=colours[1], label="Toy signal + bkg (N = {})".format(nToys))
        canvas.plot_dataset(mass_points, None, nEvtBiasToys, nEvtBiasToysErr, style="binned band", colour=colours[1], alpha=0.4)
        # Add dashed line
        canvas.plot_dataset((x_min, x_max), (0.0, 0.0), style="line join centres", colour="red", linestyle="dashed")
        canvas.set_axis_label("x", "$m_{X}$")
        canvas.set_axis_label("y", "$(N_{fit} - N_{inj}) / N_{inj}$ [%]")
        canvas.set_axis_range("x", (x_min, x_max))
        canvas.set_axis_range("y", (-10, 10))
        canvas.add_ATLAS_label(0.05, 0.96, plot_type="Simulation Internal", anchor_to="upper left")
        canvas.add_luminosity_label(0.05, 0.90, sqrts_TeV=13, luminosity=36.1, units="fb-1", anchor_to="upper left")
        canvas.add_text(0.05, 0.85, "{}-tag {} mass".format(tag_category, mass_category), anchor_to="upper left")
        canvas.add_text(0.05, 0.78, "Signal parameterisation bias", anchor_to="upper left")
        canvas.add_legend(0.05, 0.72, anchor_to="upper left")
        # canvas.add_text(0.05, 0.65, "Largest bias: {:.1f}%".format(maximum_bias), anchor_to="upper left")
        canvas.internal_header_fraction = 0.3
        canvas.save_to_file(os.path.join(base_path, "plots", "signal_bias", "signal_bias_{}Mass_{}tag".format(mass_category, tag_category)))
