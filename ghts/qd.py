"""
Contains functions to calculate q and d
as well as related quantities for sigma optimization
"""
import numpy as np
import cv
from cv_geometry import orthogonal, mod, mod2


def get_q(cv_set, ghts):
    return cv.get_cv_lincomb(cv_set, ghts['Minvn'], -np.matmul(ghts['z'], ghts['Minvn']))


def get_d(cv_set, ghts):
    ort = ort_n(cv_set.value - ghts['z'], ghts)
    d = mod(ort, ghts['Minv'])
    coef = np.matmul(ghts['Minv'], ort) / d if d > 0 else np.zeros(ort.shape, ort.dtype)
    grad = np.matmul(coef, cv_set.jacobian)
    return cv.CV(d, grad)


def get_dq_dz(ghts):
    return - ghts['Minvn']


def get_dq_dn(cv_set, ghts):
    return np.matmul(ghts['Minv'], cv_set.value - ghts['z'])


def get_grad_grad_dq_dz(ghts):
    return np.zeros(len(ghts['z']))


def get_grad_grad_dq_dn(cv_set, ghts):
    return np.matmul(ghts['Minv'], np.matmul(cv_set.m, ghts['Minvn']))


def get_gradq_mod2(cv_set, ghts):
    return mod2(ghts['Minvn'], cv_set.m)


def ort_n(vec, ghts):
    return orthogonal(vec, ghts['n'], ghts['Minv'])