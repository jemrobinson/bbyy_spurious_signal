#! /usr/bin/env python
import numpy as np

def _apply_window(_list, window_half_width):
    _list_extended = [_list[0]] * window_half_width + list(_list) + [_list[-1]] * window_half_width
    _list_groups = zip(*[_list_extended[idx:] for idx in range(2 * window_half_width + 1)])
    return [sum(_list_group) / len(_list_group) for _list_group in _list_groups]


def kNN(x, y, y_err, window_half_width=5):
    y_merged = _apply_window(y, window_half_width)
    y_err_merged = _apply_window(y_err, window_half_width)
    return x, y_merged, y_err_merged

def _gaussian_kernel(x_out, x_in, radius):
    return np.exp(-(x_out - x_in) ** 2 / (2 * radius ** 2))

def kernel(x, y, y_err, radius=2):
    y_smoothed, y_err_smoothed = [], []
    for _x_out in x:
        kernels = [_gaussian_kernel(_x_out, _x_in, radius) for _x_in in x]
        y_smoothed.append(sum([_kernel * y_in for _kernel, y_in in zip(kernels, y)]) / sum(kernels))
        y_err_smoothed.append(sum([_kernel * y_err_in for _kernel, y_err_in in zip(kernels, y_err)]) / sum(kernels))
    return x, y_smoothed, y_err_smoothed