#! /usr/bin/env python
import math
import numpy as np
import operator
from itertools import izip_longest


class DataPoint:
    """A class to hold a point with error bars"""
    def __init__(self, x_low, x_high, y, y_err):
        self.x_low, self.x_high = x_low, x_high
        self.y, self.y_err = y, y_err

    @property
    def range(self):
        return [self.x_low, self.x_high]

    @property
    def x(self):
        return 0.5 * (self.range[0] + self.range[1])

    @property
    def value(self):
        return [self.y, self.y_err]

    def significance(self):
        difference, error = abs(self.y), abs(self.y_err)
        if error == 0.0 : return 1e99
        return difference / error

    def __add__(self, other):
        if self.value[1] == 0: self.y_err = 1e-99
        if other.value[1] == 0: other.y_err = 1e-99
        x_low = min(self.range[0], other.range[0])
        x_high = max(self.range[1], other.range[1])
        y = (self.value[0] / self.value[1] + other.value[0] / other.value[1]) / (1.0 / self.value[1] + 1.0 / other.value[1])
        y_err = math.sqrt(2.0 / (1.0 / self.value[1] + 1.0 / other.value[1])**2)
        return DataPoint(x_low, x_high, y, y_err)

    def __str__(self) :
        return str([self.range, self.value])


def simple(datapoints, rebin=2):
    def grouper(iterable, n, fillvalue=None):
        args = [iter(iterable)] * n
        return izip_longest(*args, fillvalue=fillvalue)
    x, y, y_err = zip(*datapoints)
    _x, _y, _y_err = [], [], []
    if len(x) % rebin != 0:
        print "Rebinning {} points by factor {} will leave ungrouped points!".format(len(x), rebin)
    for group in grouper(zip(x, y, y_err), rebin):
        try:
            group = [elem for elem in group if elem is not None]
            _reciprocal_sqrd_err = sum([1.0/ (elem[2]**2) for elem in group])
            _x.append(sum([elem[0] for elem in group]) / len(group))
            _y.append(sum([elem[1] / (elem[2]**2) for elem in group]) / _reciprocal_sqrd_err)
            _y_err.append(1.0 / np.sqrt(_reciprocal_sqrd_err))
        except:
            print group
            raise
    return zip(_x, _y, _y_err)

def recursive(datapoints, min_sigma=0.6, bins_to_merge=None):
    # Assume that datapoints contains (position, value, err_value)
    bin_edges = [0.5 * (x_low[0] + x_high[0]) for x_low, x_high in zip(datapoints[:-1], datapoints[1:])]
    bin_edges = [2 * bin_edges[0] - bin_edges[1]] + bin_edges + [2 * bin_edges[-1] - bin_edges[-2]]
    input_data = [DataPoint(max([b for b in bin_edges if b < datapoint[0]]), min([b for b in bin_edges if b > datapoint[0]]), datapoint[1], datapoint[2]) for datapoint in datapoints]
    input_range = (min([d[0] for d in datapoints]), max([d[0] for d in datapoints]))

    # Set up initial bin list
    decorated_bins = [[input_point, idx, input_point.significance()] for idx, input_point in enumerate(input_data)]

    print input_range
    print "starting with", len(decorated_bins), "bins from", min([d[0].x for d in decorated_bins]), "to", max([d[0].x for d in decorated_bins])

    # Iterate over values merging until all bins are above min_sigma
    _bins_to_merge = []
    while True :
        if bins_to_merge is None:
            # Get lowest significance bin
            bin_lowest_sig = sorted(decorated_bins, key=operator.itemgetter(2))[0]
            lowest_sig, idx_lowest_sig = bin_lowest_sig[2], bin_lowest_sig[1]
            if lowest_sig > min_sigma: break

            # Get bins to merge
            if len(decorated_bins) == 1 : break
            adjacent_bins = [_bin for _bin in decorated_bins if abs(_bin[1] - idx_lowest_sig) == 1]
            idx_merge_bin = sorted(adjacent_bins, key=operator.itemgetter(2))[0][1]
            _bins_to_merge.append((idx_lowest_sig, idx_merge_bin))
        else:
            if len(bins_to_merge) == 0: break
            idx_lowest_sig, idx_merge_bin = bins_to_merge.pop(0)
        # print "   merging...", idx_lowest_sig, idx_merge_bin

        # Create merged data point and replace the existing bins
        new_data_point = decorated_bins[idx_lowest_sig][0] + decorated_bins[idx_merge_bin][0]
        decorated_bins[idx_lowest_sig] = [new_data_point, idx_lowest_sig, new_data_point.significance()]
        del decorated_bins[idx_merge_bin]
        for idx, bin in enumerate(decorated_bins): bin[1] = idx

    # Construct output
    output = [[_bin[0].x, _bin[0].value[0], _bin[0].value[1]] for _bin in decorated_bins]
    # Move extremal points to range ends
    output[0][0] = input_range[0]
    output[-1][0] = input_range[1]
    print "ending with", len(decorated_bins), "bins from", min([d[0].x for d in decorated_bins]), "to", max([d[0].x for d in decorated_bins])
    return output, _bins_to_merge
