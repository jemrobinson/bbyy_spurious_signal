#! /usr/bin/env python
import csv
import os
from mATLASplotlib import canvases

base_path = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
colours = ["#66a61e", "#e6ab02", "#1b9e77", "#d95f02", "#7570b3", "#e7298a"]

input_rows = []
with open(os.path.join(base_path, "output", "signal_bias.csv"), "rb") as f_input:
    for row in csv.reader(f_input, delimiter=" "):
        input_rows.append(row)

for mass_category in ["low", "high"]:
    for tag_category in ["0", "1", "2"]:
        category_rows = [[_x if "_" in _x else float(_x) for _x in row[2:]] for row in input_rows if row[0] == mass_category and row[1] == tag_category]
        if len(category_rows) == 0: continue
        # Initialise bias plot canvas
        canvas = canvases.Simple()
        maximum_bias = 0
        for idx, signal_source in enumerate(sorted(set([row[1] for row in category_rows]))):
            if signal_source != "Expected_bbyy": continue
            # if signal_source == "ATLAS_Run_2_bbbb" and mass_category == "low": continue
            # if signal_source == "ATLAS_Run_2_bbyy" and mass_category == "high": continue
            rows = [row for row in category_rows if row[1] == signal_source]
            mass_points = [row[0] for row in rows]
            nEvents_bias_percent = [100 * (row[3] - row[2]) / row[2] for row in rows]
            nEvents_bias_err_percent = [100 * row[4] / row[3] for row in rows]
            # label = "{} limits".format(signal_source.replace("_", " ").replace("yy", "$\gamma\gamma$"))
            label = "Cross-section from expected limits"
            canvas.plot_dataset((mass_points, nEvents_bias_percent), style="line join centres", colour=colours[idx], label=label)
            # canvas.plot_dataset((mass_points, None, nEvents_bias_percent, nEvents_bias_err_percent), style="binned band", colour=colours[idx], alpha=0.4)
            canvas.plot_dataset(([min(mass_points), max(mass_points)], [0.0, 0.0]), style="line join centres", colour="red", linestyle="dashed")
            maximum_bias = max([abs(x) for x in nEvents_bias_percent])
            maximum_bias = sorted(nEvents_bias_percent, key=lambda x: abs(x))[-1]
        x_min, x_max = min(mass_points), max(mass_points)
        canvas.set_axis_label("x", "$m_{X}$")
        canvas.set_axis_label("y", "$(N_{fit} - N_{inj}) / N_{inj}$ [%]")
        canvas.set_axis_range("x", (x_min, x_max))
        canvas.set_axis_range("y", (-10, 10))
        canvas.add_ATLAS_label(0.05, 0.96, plot_type="Simulation Internal", anchor_to="upper left")
        canvas.add_luminosity_label(0.05, 0.90, sqrts_TeV=13, luminosity=36.1, units="fb-1", anchor_to="upper left")
        canvas.add_text(0.05, 0.85, "{}-tag {} mass".format(tag_category, mass_category), anchor_to="upper left")
        canvas.add_legend(0.05, 0.78, anchor_to="upper left")
        canvas.add_text(0.05, 0.71, "Largest bias: {:.1f}%".format(maximum_bias), anchor_to="upper left")
        canvas.internal_header_fraction = 0.3
        canvas.save_to_file(os.path.join(base_path, "plots", "signal_bias_{}Mass_{}tag".format(mass_category, tag_category)))
