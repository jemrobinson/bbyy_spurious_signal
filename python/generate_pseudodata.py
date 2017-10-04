#! /usr/bin/env python
import csv
import logging
import numpy as np
import os
import tqdm

def ensure_path(*args):
    path = os.path.join(*args)
    if not os.path.exists(path):
        os.makedirs(path)
    return path


# Set up logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)8s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
logger.setLevel(logging.INFO)

# Input settings
expected_hh_limits_pb = {260: 1.15, 275: 1.0, 300: 0.9, 325: 0.8, 350: 0.7, 400: 0.55, 450: 0.5, 500: 0.38,  750: 0.18,  1000: 0.13}
masses = [260, 275, 300, 400, 750, None]
mass_categories = ["high", "low"]
tag_categories = ["1", "2"]
nPseudodataSamples = 100
np.random.seed(20170711)

# Create pseudodata
for idx_mass, mass in enumerate(masses, start=1):
    logger.info("Now working on mass mX = {} GeV".format(mass))
    output_dir = "bkg_only" if mass is None else "mX_{}".format(mass)
    mass_dir = ensure_path("/afs/cern.ch/work/j/jrobinso", "public", "bbyy", "pseudodata", output_dir)
    total_nBkg = dict((m, dict((t, []) for t in tag_categories)) for m in mass_categories)
    total_nSig = dict((m, dict((t, []) for t in tag_categories)) for m in mass_categories)
    for idx_sample in tqdm.trange(1, nPseudodataSamples + 1):
        sample_dir = ensure_path(mass_dir, "sample_{}".format(idx_sample))
        for mass_category in mass_categories:
            for tag_category in tag_categories:
                f_output_name = os.path.join(sample_dir, "pseudodata_myyjj_{}_sample_{}_{}Mass_{}tag.txt".format(output_dir, idx_sample, mass_category, tag_category))
                if os.path.isfile(f_output_name):
                    logger.warning("File {} exists! Skipping".format(f_output_name))
                    continue
                with open(f_output_name, "wb") as f_output:
                    nBkg, nSig = 0, 0
                    # Background
                    with open(os.path.join("input", "m_yyjj_SM_bkg_{}Mass_{}tag_tightIsolated.csv".format(mass_category, tag_category)), "rb") as f_input:
                        for row in csv.reader(f_input, delimiter="\t"):
                            if np.random.random() < float(row[1]):
                                f_output.write(row[0] + "\n")
                                nBkg += 1
                    # Signal
                    if mass is not None:
                        hh_xs_pb = (expected_hh_limits_pb[mass]) * 2
                        weight_scale = float(hh_xs_pb) / 5. # default weights are for a 5 pb hh cross-section
                        with open(os.path.join("input", "m_yyjj_Xhh_m{}_{}Mass_{}tag_tightIsolated.csv".format(mass, mass_category, tag_category)), "rb") as f_input:
                            for row in csv.reader(f_input, delimiter="\t"):
                                if np.random.random() < (float(row[1]) * weight_scale):
                                    f_output.write(row[0] + "\n")
                                    nSig += 1
                    logger.debug("... => {} mass {}-tag pseudodata with {} background events and {} signal events".format(mass_category, tag_category, nBkg, nSig))
                    total_nBkg[mass_category][tag_category].append(nBkg)
                    total_nSig[mass_category][tag_category].append(nSig)
    for mass_category in mass_categories:
        for tag_category in tag_categories:
            nBkg = np.mean(total_nBkg[mass_category][tag_category])
            nSig = np.mean(total_nSig[mass_category][tag_category])
            logger.info("... overall average {} mass {}-tag pseudodata is {} background events and {} signal events".format(mass_category, tag_category, nBkg, nSig))







# datasets = ["bkg_only", "bkg_only", "bkg_only", "m260:5000", "m260:2100", "m300:2000", "m350:1500", "m500:600", "m750:100", "m1000:40"]
# datasets = ["m300:5000",  "m750:5000", "bkg_only", "bkg_only", "bkg_only", "m260:350", "m275:1550", "m300:1200", "m325:700", "m350:700", "m400:400", "m450:600", "m500:350", "m750:150", "m1000:40"]
# datasets = ["m300:5000",  "m750:5000", "bkg_only", "bkg_only", "bkg_only", "m260:350", "m275:1550", "m300:1200", "m325:700", "m350:700", "m400:400", "m450:600", "m500:350", "m750:150", "m1000:40"]

# datasets = ["260:2.3", "275:2.0", "300:1.8", "400:1.1", "750:0.36"]
# np.random.seed(20170711)

# for idx, dataset in enumerate(np.random.permutation(datasets), start=1):
#     print "** Generating sample {}: {} **".format(idx, dataset)
#     for mass_category in ["high", "low"]:
#         for tag_category in ["0", "1", "2"]:
#             f_output_name = os.path.join("output", "pseudodata", "pseudodata_sample_{}_{}Mass_{}tag.txt".format(idx, mass_category, tag_category))
#             if os.path.isfile(f_output_name):
#                 print "File {} exists! Skipping".format(f_output_name)
#                 continue
#             with open(f_output_name, "wb") as f_output:
#                 nBkg, nSig = 0, 0
#                 # Background
#                 with open(os.path.join("input", "m_yyjj_SM_bkg_{}Mass_{}tag_tightIsolated.csv".format(mass_category, tag_category)), "rb") as f_input:
#                     for row in csv.reader(f_input, delimiter="\t"):
#                         # if random.random() < float(row[1]):
#                         if np.random.random() < float(row[1]):
#                             f_output.write(row[0] + "\n")
#                             nBkg += 1
#                 # Signal
#                 if ":" in dataset:
#                     mass, hh_xs_fb = dataset.split(":")
#                     weight_scale = float(xs_fb) / 5000.
#                     with open(os.path.join("input", "m_yyjj_Xhh_{}_{}Mass_{}tag_tightIsolated.csv".format(mass, mass_category, tag_category)), "rb") as f_input:
#                         for row in csv.reader(f_input, delimiter="\t"):
#                             # if random.random() < (float(row[1]) * weight_scale):
#                             if np.random.random() < (float(row[1]) * weight_scale):
#                                 f_output.write(row[0] + "\n")
#                                 nSig += 1
#                 print "Generated {} mass {}-tag pseudodata with {} background events and {} signal events".format(mass_category, tag_category, nBkg, nSig)

# # # Plot distributions
# # for idx, _ in enumerate(datasets, start=1):
# #     for mass_category in ["high", "low"]:
# #         x_range = {"high":(335, 1140), "low":(245, 485)}[mass_category]
# #         for tag_category in ["0", "1", "2"]:
# #             data = np.loadtxt(os.path.join("output", "pseudodata", "pseudodata_sample_{}_{}Mass_{}tag.txt".format(idx, mass_category, tag_category)))
# #             hist, _edges = np.histogram(data, "auto", range=x_range)
# #             # hist, _edges = np.histogram(data, int(data.size / 10 + 1), range=x_range)
# #             centres = [0.5 * (e1 + e2) for e1, e2 in zip(_edges[:-1], _edges[1:])]
# #             canvas = canvases.Simple()
# #             canvas.plot_dataset((centres, None, hist, np.sqrt(hist)), style="scatter yerror", colour="black")
# #             canvas.set_axis_label("x", "$m_{\gamma\gamma jj}$")
# #             canvas.set_axis_label("y", "Events / bin")
# #             canvas.set_axis_range("x", x_range)
# #             canvas.save_to_file(os.path.join("plots", "pseudodata_sample_{}_{}Mass_{}tag".format(idx, mass_category, tag_category)))
