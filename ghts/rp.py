"""
RP dependent functions of cartesian coordinates
"""
import numpy as np
import cv
import qd


def get_sigma(sigma, q):
    return cv.CV(
        value=sigma.value,
        gradient=np.array([dq*dsigma for dsigma, dq
                            in zip(sigma.gradient, q.gradient)])
    )


def get_q(ghts, cv_set):
    # Here q.gradient is not formally dq/dx, since q is a vector and x is a matrix, so
    # dq/dx would be a 3D tensor. Instead q.gradient[i,j] = dq[i]/dx[i][j]
    q_beads = map_beads(qd.get_q, cv_set, ghts)
    return cv.CV(
        value=np.array([q_bead.value for q_bead in q_beads]),
        gradient=np.array([q_bead.gradient for q_bead in q_beads])
    )


def get_d(ghts, cv_set):
    d_beads = map_beads(qd.get_d, cv_set, ghts)
    return cv.CV(
        value=np.average([d_bead.value for d_bead in d_beads], 0),
        gradient=np.array([d_bead.gradient for d_bead in d_beads]) / len(cv_set)
    )


def get_dq_dz(ghts, cv_set):
    return np.array([qd.get_dq_dz(ghts) for cv_set_bead in cv_set])


def get_dq_dn(ghts, cv_set):
    return np.array(map_beads(qd.get_dq_dn, cv_set, ghts))


def get_grad_grad_dq_dz(ghts, cv_set):
    return np.array([qd.get_grad_grad_dq_dz(ghts) for cv_set_bead in cv_set])


def get_grad_grad_dq_dn(ghts, cv_set):
    return np.array(map_beads(qd.get_grad_grad_dq_dn, cv_set, ghts))


def get_gradq_mod2(ghts, cv_set):
    return np.array(map_beads(qd.get_gradq_mod2, cv_set, ghts))


def mean_beads(callback, cv_set, *args):
    return np.average(map_beads(callback, cv_set, *args), 0)


def map_beads(callback, cv_set, *args):
    return [callback(cv_set_bead, *args) for cv_set_bead in cv_set]