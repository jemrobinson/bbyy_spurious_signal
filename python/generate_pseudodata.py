#! /usr/bin/env python
import csv
import logging
import numpy as np
import os
import sys
import tqdm

def ensure_path(*args):
    path = os.path.join(*args)
    if not os.path.exists(path):
        os.makedirs(path)
    return path

# hh_injected_xs_pb_hypo_0 = {260: 1.15, 275: 1.0, 300: 0.9, 325: 0.8, 350: 0.7, 400: 0.55, 450: 0.5, 500: 0.38,  750: 0.18,  1000: 0.13}
# hh_injected_xs_pb_hypo_1 = {260: 1.725, 275: 1.5, 300: 1.35, 325: 1.2, 350: 1.05, 400: 0.825, 450: 0.75, 500: 0.57,  750: 0.27,  1000: 0.195}
# hh_injected_xs_pb_hypo_2 = {260: 2.3, 275: 2.0, 300: 1.8, 325: 1.6, 350: 1.4, 400: 1.1, 450: 1.0, 500: 0.76,  750: 0.36,  1000: 0.26}

# Set up logging
logger = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, format="%(asctime)s %(levelname)8s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
logger.setLevel(logging.INFO)

# Input settings
version_directory = "20171212"
expected_hh_limits_pb = {260: 1.15, 275: 1.0, 300: 0.9, 325: 0.8, 350: 0.7, 400: 0.55, 450: 0.5, 500: 0.38,  750: 0.18,  1000: 0.13}
masses = [260, 275, 300, 325, 350, 400, 450, 500, 750, 1000, None]
mass_categories = ["high", "low"]
tag_categories = ["1", "2"]
signal_hypos = [1, 1.5, 2]
nPseudodataSamples = 1000
np.random.seed(20171212)

# Create pseudodata
for idx_mass, mass in enumerate(masses, start=1):
    logger.info("Now working on mass mX = {} GeV".format(mass))
    output_dir = "bkg_only" if mass is None else "mX_{}".format(mass)
    mass_dir = ensure_path("/afs/cern.ch/work/j/jrobinso", "public", "bbyy", version_directory, "blinded", "pseudodata", output_dir)
    for idx_signal_hypo, signal_hypo in enumerate(signal_hypos):
        total_nBkg = dict((m, dict((t, []) for t in tag_categories)) for m in mass_categories)
        total_nSig = dict((m, dict((t, []) for t in tag_categories)) for m in mass_categories)
        signal_hypo_dir = ensure_path(mass_dir, "signal_hypothesis_{}".format(idx_signal_hypo))
        logger.info("... signal hypothesis {}/{}".format(idx_signal_hypo + 1, len(signal_hypos)))
        # Iterate over tag and mass categories for a given signal hypothesis
        sum_weights_bkg, sum_weights_sig = {}, {}
        for mass_category in mass_categories:
            for tag_category in tag_categories:
                # Load background events
                masses_bkg, weights_bkg = [], []
                with open(os.path.join("input", "m_yyjj_SM_bkg_{}Mass_{}tag_tightIsolated_positive_weights.csv".format(mass_category, tag_category)), "rb") as f_input:
                    for row in csv.reader(f_input, delimiter="\t"):
                        masses_bkg.append(float(row[0]))
                        weights_bkg.append(float(row[1]))
                sum_weights_bkg[(mass_category, tag_category)] = sum(weights_bkg)
                logger.debug("... => {} mass {}-tag pseudodata has {} background entries, corresponding to {:.2f} events".format(mass_category, tag_category, len(masses_bkg), sum_weights_bkg[(mass_category, tag_category)]))
                nEventsPoisson_bkg = np.random.poisson(sum_weights_bkg[(mass_category, tag_category)], nPseudodataSamples)
                normalised_weights_bkg = [_w / sum_weights_bkg[(mass_category, tag_category)] for _w in weights_bkg]
                # Load signal events
                masses_sig, weights_sig = [], []
                if mass is not None:
                    hh_xs_pb = (expected_hh_limits_pb[mass]) * signal_hypo
                    weight_scale = float(hh_xs_pb) / 5. # default weights are for a 5 pb hh cross-section
                    with open(os.path.join("input", "m_yyjj_Xhh_m{}_{}Mass_{}tag_tightIsolated_positive_weights.csv".format(mass, mass_category, tag_category)), "rb") as f_input:
                        for row in csv.reader(f_input, delimiter="\t"):
                            masses_sig.append(float(row[0]))
                            weights_sig.append(float(row[1]))
                sum_weights_sig[(mass_category, tag_category)] = sum(weights_sig)
                logger.debug("... => {} mass {}-tag pseudodata has {} signal entries, corresponding to {:.2f} events".format(mass_category, tag_category, len(masses_sig), sum_weights_sig[(mass_category, tag_category)]))
                nEventsPoisson_sig = np.random.poisson(sum_weights_sig[(mass_category, tag_category)], nPseudodataSamples)
                normalised_weights_sig = [_w / sum_weights_sig[(mass_category, tag_category)] for _w in weights_sig]
                # Write out pseudodata
                for idx_sample in tqdm.trange(1, nPseudodataSamples + 1):
                    sample_dir = ensure_path(signal_hypo_dir, "sample_{}".format(idx_sample))
                    f_output_name = os.path.join(sample_dir, "pseudodata_myyjj_{}_sample_{}_{}Mass_{}tag.txt".format(output_dir, idx_sample, mass_category, tag_category))
                    if os.path.isfile(f_output_name):
                        logger.warning("File {} exists! Skipping".format(f_output_name))
                        continue
                    with open(f_output_name, "wb") as f_output:
                        nBkg, nSig = nEventsPoisson_bkg[idx_sample - 1], nEventsPoisson_sig[idx_sample - 1]
                        for _myyjj in np.random.choice(masses_bkg, nBkg, replace=False, p=normalised_weights_bkg):
                            f_output.write("{}\n".format(_myyjj))
                        if mass is not None:
                            for _myyjj in np.random.choice(masses_sig, nBkg, replace=False, p=normalised_weights_sig):
                                f_output.write("{}\n".format(_myyjj))
                        for _myyjj in events_bkg + events_sig:
                        logger.debug("... => {} mass {}-tag pseudodata with {} background events and {} signal events".format(mass_category, tag_category, nBkg, nSig))
                        total_nBkg[mass_category][tag_category].append(nBkg)
                        total_nSig[mass_category][tag_category].append(nSig)
        for mass_category in mass_categories:
            for tag_category in tag_categories:
                nBkg = np.mean(total_nBkg[mass_category][tag_category])
                nSig = np.mean(total_nSig[mass_category][tag_category])
                logger.info("... overall average {} mass {}-tag pseudodata is {:.3f} background events (expected {:.3f}) and {:.3f} signal events (expected {:.3f})".format(mass_category, tag_category, nBkg, sum_weights_bkg[(mass_category, tag_category)], nSig, sum_weights_sig[(mass_category, tag_category)]))
