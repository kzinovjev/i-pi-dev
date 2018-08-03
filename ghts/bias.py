import numpy as np


def harmonic_bias(value, force_constant, reference):
    return harmonic_bias_force(value, force_constant, reference)


def side_harmonic_bias(value, force_constant, reference, side=1):
    if side is not None and (value.value - reference) * side < 0:
        return np.zeros(value.gradient.shape, value.gradient.dtype)
    return harmonic_bias_force(value, force_constant, reference)


def harmonic_bias_force(value, force_constant, reference):
    return - value.gradient * (value.value - reference) * force_constant
