#! /usr/bin/env python
import csv
import math
import numpy as np
import os
from collections import defaultdict, OrderedDict
from mATLASplotlib import canvases


colours = {"novosibirsk": "#e7298a", "modified_gamma": "#1b9e77", "modified_landau": "#7570b3",
           "exppoly1": "#e7298a", "exppoly2": "#1b9e77", "invpoly2": "#7570b3", "invpoly3": "#d95f02"}
base_path = os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
resonance_range = {"low": (260, 400), "high": (400, 1000)}
RATIO_UNCERTAINTY = 10
SPIKE_CUT_OFF = 0.2

for mass_category in ["high", "low"]:
    if mass_category == "low":
        functions = ["novosibirsk", "modified_gamma", "modified_landau"]
    else:
        functions = ["exppoly1", "exppoly2", "invpoly2", "invpoly3"]

    for tag_category in ["0", "1", "2"]:
        # Check that input file exists
        input_path = os.path.join(base_path, "output", "csv", "spurious_signal_{}Mass_{}tag.csv".format(mass_category, tag_category))
        if not os.path.exists(input_path):
            print "{} does not exist".format(input_path)
            continue

        n_spurious = defaultdict(list)
        Z_spurious = defaultdict(list)
        chi_squared = {}
        mass_points = []

        # Read inputs from csv
        with open(input_path, "rb") as f_input:
            for row in csv.reader(f_input, delimiter=" "):
                try:
                    function, mass, n_sig, Z, Z_uncertainty, chi2, ndof = row
                    mass, n_sig, Z, Z_uncertainty, chi2, ndof = float(mass), float(n_sig), float(Z), float(Z_uncertainty), float(chi2), int(ndof)
                except ValueError:
                    print row
                    raise
                if resonance_range[mass_category][0] <= mass <= resonance_range[mass_category][1]:
                    if not np.isinf(Z) and not np.isinf(Z_uncertainty):
                        if Z_uncertainty == 0.0: Z_uncertainty = 0.5 * Z
                        n_spurious[function].append([mass, n_sig, 0.0])
                        Z_spurious[function].append([mass, Z, Z_uncertainty])
                    else:
                        print "Rejecting", function, "at", mass, "GeV because fit did not converge"
                    mass_points.append(float(mass))
                else:
                    chi_squared[function] = (chi2, ndof)

        # Remove points where Z is more than SPIKE_CUT_OFF away from the points around it
        for function in functions:
            Z_values = [abs(point[1]) for point in Z_spurious[function]]
            pathological_filter = [abs((Z_value - 0.5 * (low + high))) < SPIKE_CUT_OFF for low, Z_value, high in zip(Z_values[0:1] + Z_values[:-1], Z_values, Z_values[1:] + Z_values[-2:-1])]
            Z_spurious[function] = [_elem for _elem, _filter in zip(Z_spurious[function], pathological_filter) if _filter]
            n_spurious[function] = [_elem for _elem, _filter in zip(n_spurious[function], pathological_filter) if _filter]

        # Cap Z uncertainties at +/- 5 to avoid matplotlib issue
        for function in functions:
            for idx_point in range(len(Z_spurious[function])):
                Z_spurious[function][idx_point][2] = min(Z_spurious[function][idx_point][2], 5.0)

        # Sort by mass
        [tup.sort(key=lambda x:x[0]) for tup in n_spurious.values()]
        [tup.sort(key=lambda x:x[0]) for tup in Z_spurious.values()]

        # Get x-axis min and max for plotting
        if len(mass_points) == 0:
            continue
        x_min, x_max = min(mass_points), max(mass_points)

        # Ensure that output directory exists
        if not os.path.exists(os.path.join(base_path, "output", "plots")):
            os.makedirs(os.path.join(base_path, "output", "plots"))

        # Plot n_spurious
        canvas = canvases.Simple()
        for function in functions:
            x, y, y_err = zip(*n_spurious[function])
            canvas.plot_dataset((x, y), style="line join centres", label=function.replace("_", " ").title() + ": $N_{{spur}}^{{max}}$ = {:.2f} ".format(max([abs(_y) for _y in y])), colour=colours[function])
        canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
        canvas.set_axis_label("y", "$N_{spurious}$")
        canvas.set_axis_range("x", (x_min, x_max))
        canvas.add_ATLAS_label(0.05, 0.96, plot_type="Simulation Internal", anchor_to="upper left")
        canvas.add_luminosity_label(0.05, 0.90, sqrts_TeV=13, luminosity=36.1, units="fb-1", anchor_to="upper left")
        canvas.add_text(0.05, 0.85, "{}-tag {} mass".format(tag_category, mass_category), anchor_to="upper left")
        canvas.add_legend(0.97, 0.85, fontsize="medium", anchor_to="upper right")
        canvas.internal_header_fraction = 0.3
        canvas.save_to_file(os.path.join(base_path, "plots", "n_spurious_{}Mass_{}tag".format(mass_category, tag_category)))

        # Plot Z_spurious
        canvas = canvases.Simple()
        y_min, y_max = -0.3, 0.3
        for function in functions:
            x, y, y_err = zip(*Z_spurious[function])
            canvas.plot_dataset((x, None, y, y_err), style="binned band", colour=colours[function], alpha=0.2)
            canvas.plot_dataset((x, y), style="line join centres", label=function.replace("_", " ").title(), colour=colours[function])
            y_restricted = [_y for _y in y if -1.5 < float(_y) < 1.0]
            y_min, y_max = min(y_min, min(y_restricted)), max(y_max, max(y_restricted))
        canvas.plot_dataset(([x_min, x_max], [0.2, 0.2]), style="line join centres", colour="red", linestyle="dashed")
        canvas.plot_dataset(([x_min, x_max], [-0.2, -0.2]), style="line join centres", colour="red", linestyle="dashed")
        canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
        canvas.set_axis_label("y", "$S / \Delta S$")
        canvas.set_axis_range("x", (x_min, x_max))
        canvas.set_axis_range("y", (1.1 * y_min, 1.1 * y_max))
        canvas.add_ATLAS_label(0.05, 0.96, plot_type="Simulation Internal", anchor_to="upper left")
        canvas.add_luminosity_label(0.05, 0.90, sqrts_TeV=13, luminosity=36.1, units="fb-1", anchor_to="upper left")
        canvas.add_text(0.05, 0.85, "{}-tag {} mass".format(tag_category, mass_category), anchor_to="upper left")
        canvas.internal_header_fraction = 0.3
        canvas.add_legend(0.97, 0.85, fontsize="medium", anchor_to="upper right")
        canvas.save_to_file(os.path.join(base_path, "plots", "Z_spurious_{}Mass_{}tag".format(mass_category, tag_category)))

        if not os.path.exists("tex"):
            os.makedirs("tex")
        with open(os.path.join("tex", "spurious-myyjj-{}tag-{}.tex".format(tag_category, mass_category)), "wb") as f_output:
            f_output.write("\\begin{table}\\footnotesize\n")
            f_output.write("\\begin{center}\n")
            f_output.write("\\caption{{Performance of different functions that are considered for the background modeling of the {} b-tag category ({} mass selection).\n".format(tag_category, mass_category))
            f_output.write("The associated systematic uncertainty on the signal amplitude in terms of spurious signal $N_{{spur}}$ and its ratio to the statistical uncertainty on the fitted number of signal events, $Z_{{spur}}$ computed using the background estimate described in Section~\ref{sec:background}.\n")
            f_output.write("The $\chi^2$/ndof obtained from performing a background-only fit is also shown.\n")
            f_output.write("The functions used are defined above and the number of free parameters for each function is given by nPars.}\n")
            f_output.write("\\label{{spurious-myyjj-{}tag-{}}}\n".format(tag_category, mass_category))
            f_output.write("\\begin{tabular}{|c|c|c|c|c|c|}\n")
            f_output.write("\\hline\n")
            f_output.write("\\textbf{Model}             & \\textbf{max(Z\_spur) {[}\%{]}} & \\textbf{max(N\_spur)} & \\textbf{nPars} & \\textbf{$\chi^2$/ndof}\\\\\n")
            f_output.write("\\hline\n")
            if "novosibirsk" in functions:
                f_output.write("\\textbf{{Novosibirsk}}     & {:.2f}                     & {:.2f}            & 3              & {:.1f} / {}                 \\\\\n".format(100 * max(abs(x[1]) for x in Z_spurious["novosibirsk"]), max(abs(x[1]) for x in n_spurious["novosibirsk"]), chi_squared["novosibirsk"][0], chi_squared["novosibirsk"][1]))
            if "modified_gamma" in functions:
                f_output.write("\\textbf{{Modified Gamma}}  & {:.2f}                     & {:.2f}            & 5              & {:.1f} / {}                 \\\\\n".format(100 * max(abs(x[1]) for x in Z_spurious["modified_gamma"]), max(abs(x[1]) for x in n_spurious["modified_gamma"]), chi_squared["modified_gamma"][0], chi_squared["modified_gamma"][1]))
            if "modified_landau" in functions:
                f_output.write("\\textbf{{Modified Landau}} & {:.2f}                     & {:.2f}            & 3              & {:.1f} / {}                 \\\\\n".format(100 * max(abs(x[1]) for x in Z_spurious["modified_landau"]), max(abs(x[1]) for x in n_spurious["modified_landau"]), chi_squared["modified_landau"][0], chi_squared["modified_landau"][1]))
            if "exppoly1" in functions:
                f_output.write("\\textbf{{Exppoly1}}        & {:.2f}                     & {:.2f}            & 1              & {:.1f} / {}                 \\\\\n".format(100 * max(abs(x[1]) for x in Z_spurious["exppoly1"]), max(abs(x[1]) for x in n_spurious["exppoly1"]), chi_squared["exppoly1"][0], chi_squared["exppoly1"][1]))
            if "exppoly2" in functions:
                f_output.write("\\textbf{{Exppoly2}}        & {:.2f}                     & {:.2f}            & 1              & {:.1f} / {}                 \\\\\n".format(100 * max(abs(x[1]) for x in Z_spurious["exppoly2"]), max(abs(x[1]) for x in n_spurious["exppoly2"]), chi_squared["exppoly2"][0], chi_squared["exppoly2"][1]))
            if "invpoly2" in functions:
                f_output.write("\\textbf{{Invpoly2}}        & {:.2f}                     & {:.2f}            & 1              & {:.1f} / {}                 \\\\\n".format(100 * max(abs(x[1]) for x in Z_spurious["invpoly2"]), max(abs(x[1]) for x in n_spurious["invpoly2"]), chi_squared["invpoly2"][0], chi_squared["invpoly2"][1]))
            if "invpoly3" in functions:
                f_output.write("\\textbf{{Invpoly3}}        & {:.2f}                     & {:.2f}            & 1              & {:.1f} / {}                 \\\\\n".format(100 * max(abs(x[1]) for x in Z_spurious["invpoly3"]), max(abs(x[1]) for x in n_spurious["invpoly3"]), chi_squared["invpoly3"][0], chi_squared["invpoly3"][1]))
            f_output.write("\\hline\n")
            f_output.write("\\end{tabular}\n")
            f_output.write("\\end{center}\n")
            f_output.write("\\end{table}\n")
