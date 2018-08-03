import numpy as np
from ipi.utils.nmtransform import mk_nm_matrix
from collections import namedtuple

CV2 = namedtuple('CV2', ['value', 'gradient', 'hessian'])


def get_gradsigma_mod(sigma, gradq_mod2):
    return np.sqrt(np.dot(sigma.gradient**2, gradq_mod2))


def get_sigma(modes, r):
    return CV2(np.dot(modes['a'], r.value - modes['b']), np.matmul(modes['a'], r.gradient),
               sum(d2r_k * a_k for d2r_k, a_k in zip(r.hessian, modes['a'])))


def get_r(q, nmodes):
    eta = get_eta(q)
    nbeads = len(eta)
    r = _get_r(eta, nmodes)
    return CV2(r / np.sqrt(nbeads),
               _get_dr(r, eta) / np.sqrt(nbeads),
               _get_d2r(r, eta) / np.sqrt(nbeads))


def _get_r(eta, nmodes):
    nbeads = len(eta)
    r = np.zeros(nbeads/2 if nbeads > 1 else 1)
    r[0] = eta[0]
    for i in range(1, (nbeads - 1) / 2):
        r[i] = np.sqrt(eta[k(i)] ** 2 + eta[_k(i, nbeads)] ** 2)
    if nbeads % 2 == 0:
        r[nbeads / 2 - 1] = eta[nbeads - 1]
    return r[0:nmodes]


def _get_dr(r, eta):
    nr = len(r)
    nbeads = len(eta)
    dr = np.zeros((nr, nbeads))
    U = mk_nm_matrix(nbeads)
    dr[0] = U[0]
    for i in range(1, min(nr, (nbeads - 1) / 2)):
        dr[i] = (U[k(i)] * eta[k(i)] + U[_k(i, nbeads)] * eta[_k(i, nbeads)])/r[i]
    if nr == (nbeads + 1) / 2 + 1:
        dr[nbeads / 2] = U[nbeads]
    return dr


def _get_d2r(r, eta):
    nr = len(r)
    nbeads = len(eta)
    d2r = np.zeros((nr, nbeads, nbeads))
    U = mk_nm_matrix(nbeads)
    for i in range(1, nr):
        d2r[i] = (
                     out(U[k(i)]) + out(U[_k(i, nbeads)]) -
                     out(U[k(i)] * eta[k(i)] - U[_k(i, nbeads)] * eta[_k(i, nbeads)])
                 ) / r[i]
    return d2r


def get_eta(q):
    return np.matmul(mk_nm_matrix(len(q.value)), q.value)


def out(v):
    return np.outer(v, v)


def k(i):
    return i


def _k(i, nbeads):
    return nbeads - i
