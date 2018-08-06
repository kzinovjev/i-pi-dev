"""
Responsible for RC optimization (changes b, a, z, n and M)
"""
from mpi4py import MPI
import nm
from increments import *
import rp


def move(modes, ghts, cv_set, r, sigma, params, bk, dt):
    gradq_mod2 = rp.get_gradq_mod2(ghts, cv_set)
    gradsigma_mod = nm.get_gradsigma_mod(sigma, gradq_mod2)

    if params['move_modes']:
        move_modes(modes, sigma, r, gradsigma_mod, gradq_mod2, params, bk, dt)

    if params['move_ghts']:
        move_ghts(ghts, cv_set, sigma, gradsigma_mod, gradq_mod2, params, bk, dt)

    if not params['fix_M']:
        ghts['M'] += mpi_reduce(m_increment(ghts, cv_set)) / params['gamma_M'] * dt


def mpi_reduce(inc):
    return MPI.COMM_WORLD.allreduce(inc) / MPI.COMM_WORLD.Get_size()


def move_modes(modes, sigma, r, gradsigma_mod, gradq_mod2, params, bk, dt):
    modes['b'] += mpi_reduce(b_increment(sigma, gradsigma_mod,
                                         modes['a'], bk)) / params['gamma_b'] * dt
    modes['b'] += mpi_reduce(delta_increment(r, modes['a'],
                                             modes['b'])) / params['gamma_delta'] * dt
    modes['a'] += mpi_reduce(a_increment(sigma, r, gradsigma_mod,
                                         gradq_mod2, modes['b'], bk)) / params['gamma_a'] * dt


def move_ghts(ghts, cv_set, sigma, gradsigma_mod, gradq_mod2, params, bk, dt):
    ghts['z'] += mpi_reduce(z_increment(ghts, cv_set, sigma,
                                        gradsigma_mod, gradq_mod2, bk)) / params['gamma_z'] * dt
    ghts['z'] += mpi_reduce(d_increment(ghts, cv_set)) / params['gamma_d'] * dt
    ghts['n'] += mpi_reduce(n_increment(ghts, cv_set, sigma,
                                        gradsigma_mod, gradq_mod2, bk)) / params['gamma_n'] * dt