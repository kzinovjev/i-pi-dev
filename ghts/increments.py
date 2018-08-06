import numpy as np
from cv_geometry import orthogonal
import rp

def b_increment(sigma, gradsigma_mod, a, bk):
    return - bk * sigma.value * a * gradsigma_mod


def a_increment(sigma, r, gradsigma_mod, gradq_mod2, b, bk):
    return bk * sigma.value * (r.value - b) * gradsigma_mod - \
           np.matmul(np.matmul(r.gradient, np.diag(gradq_mod2)), sigma.gradient) / gradsigma_mod


def delta_increment(r, a, b):
    return (r.value - b) - np.dot(r.value - b, a) * a


def z_increment(ghts, cv_set, sigma, gradsigma_mod, gradq_mod2, bk):
    return c_increment(ghts, sigma, gradsigma_mod, gradq_mod2, rp.get_dq_dz(ghts, cv_set),
                       rp.get_grad_grad_dq_dz(ghts, cv_set), bk)


def n_increment(ghts, cv_set, sigma, gradsigma_mod, gradq_mod2, bk):
    return ort_n(c_increment(ghts, sigma, gradsigma_mod, gradq_mod2, rp.get_dq_dn(ghts, cv_set),
                       rp.get_grad_grad_dq_dn(ghts, cv_set), bk), ghts)


def d_increment(ghts, cv_set):
    mean_cv_value = rp.mean_beads(getattr, cv_set, 'value')
    return ort_n(mean_cv_value - ghts['z'], ghts)


def m_increment(ghts, cv_set):
    return rp.mean_beads(getattr, cv_set, 'm') - ghts['M']


def c_increment(ghts, sigma, gradsigma_mod, gradq_mod2, dq_dc, grad_grad_dq_dc, bk):
    return np.matmul(ghts['M'], bk * sigma.value * gradsigma_mod * np.matmul(sigma.gradient, dq_dc) - \
           grad_grad_dsigma_dc(sigma, gradq_mod2, dq_dc, grad_grad_dq_dc) / gradsigma_mod)


def grad_grad_dsigma_dc(sigma, gradq_mod2, dq_dc, grad_grad_dq_dc):
    return np.matmul(sigma.gradient * gradq_mod2, np.matmul(sigma.hessian, dq_dc)) + \
           np.matmul(sigma.gradient**2, grad_grad_dq_dc)


def ort_n(vec, ghts):
    return orthogonal(vec, ghts['n'], ghts['Minv'])