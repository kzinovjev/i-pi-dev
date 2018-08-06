import numpy as np
import cv
import rp
import qd
import nm
from cv_geometry import normal
import optimizer
from bias import harmonic_bias, side_harmonic_bias
from ipi.utils.units import UnitMap
from math import sqrt
from ipi.utils.io import print_file_path
from ipi.utils.softexit import softexit

out_file = open('CV.out', 'w')
bead_out_files = []

PS = 1E12 / UnitMap['time']['second']
AA = 1. / UnitMap['length']['angstrom']
KCAL_MOL = 0.001 / UnitMap['energy']['cal/mol']
AMU = 1. / UnitMap['mass']['dalton']
SQAMU = sqrt(AMU)
MWAA = SQAMU * AA


def force(beads, cell, masses, temp, dt, state):
    global bead_out_files

    nbeads = len(beads.q)
    if len(bead_out_files) == 0:
        bead_out_files = [open('CV_' + str(i+1) + '.out', 'w') for i in range(nbeads)]

    if len(state['modes']['a']) > (nbeads + 1) / 2:
        state['modes']['a'] = np.resize(state['modes']['a'], (nbeads + 1) / 2)

    cv_set = [cv.get_cv_set(q, state['CV'], masses) for q in beads.q]
    if state['ghts'].get('M') is None:
        state['ghts']['M'] = np.average([bead_cv_set.m for bead_cv_set in cv_set], 0) / AMU

    params = convert_params(state['params'])
    modes = convert_modes(state['modes'])
    nmodes = len(modes['a'])
    ghts = convert_ghts(state['ghts'])

    q = rp.get_q(ghts, cv_set)
    r = nm.get_r(q, nmodes)
    sigma = nm.get_sigma(modes, r)
    d = rp.get_d(ghts, cv_set)

    stage = state['stage']
    stage['step'] = stage.get('step', 0) + 1

    write_centroid_data(cv_set, sigma, q, d, r, out_file)
    write_bead_data(cv_set, ghts, bead_out_files)

    if stage['name'] == 'optimize':
        optimizer.move(modes, ghts, cv_set, r, sigma, params, params['K'] / temp, dt)

    if stage['name'] == 'sample':
        stage['last_save_step'] = stage.get('last_save_step', stage['step'])
        stage['last_saved'] = stage.get('last_saved', (stage['walker']-1) * stage['structures'])
        if abs(sigma.value * SQAMU) < stage['q_threshold'] and \
           stage['step'] >= stage['last_save_step'] + stage['offset']:
            idx = stage['last_saved'] + 1
            with open(str(idx) + ".xyz", 'w') as f:
                print_file_path('xyz', beads, cell, f, units='angstrom')
            stage['last_saved'] = idx
            stage['last_save_step'] = stage['step']
            if idx == stage['walker']*stage['structures']:
                softexit.trigger('exit ' + str(stage['walker']))

    if stage['name'] == 'prepare' and stage['step'] < stage['prepare_steps']:
        params['K'] *= stage['step'] / stage['prepare_steps']
        params['K_d'] *= stage['step'] / stage['prepare_steps']
        if stage['step'] == state['prepare_steps']:
            stage['name'] = 'optimize'
            stage['step'] = 0

    sigma_bias = harmonic_bias(rp.get_sigma(sigma, q), params['K']*nbeads, 0)
    d_bias = side_harmonic_bias(d, params['K_d']*nbeads, params['d_max'])

    restraint_biases = np.zeros(beads.q.shape, beads.q.dtype)
    for restraint in state['restraints']:
            restraint_biases += np.array([restraint_bias(bead, restraint) for bead in beads.q])

    if stage['name'] == 'committor':
        if abs(sigma.value * SQAMU) > stage['q_threshold']:
            softexit.trigger('q_threshold reached')
        return sigma_bias * 0

    recover_ghts(state, ghts)
    recover_modes(state, modes)

    return sigma_bias + d_bias + restraint_biases


def restraint_bias(bead, bias_def):
    return side_harmonic_bias(cv.get_cv(bead, bias_def),
                              bias_def['K'] / KCAL_MOL,
                              bias_def['reference'],
                              bias_def.get('side', None))


def write_centroid_data(cv_set, sigma, q, d, r, out_file):
    mean_cv_value = rp.mean_beads(getattr, cv_set, 'value')
    nitems = len(mean_cv_value) + len(r.value) + 3
    sigma_value = sigma.value * SQAMU
    q_value = np.mean(q.value) * SQAMU
    d_value = d.value * SQAMU
    r_value = r.value * SQAMU
    data = [sigma_value, q_value, d_value] + list(r_value) + list(mean_cv_value)
    out_file.write(('{:>12.4e}' * nitems + '\n').format(*data))
    out_file.flush()


def write_bead_data(cv_set, ghts, bead_out_files):
    for bead_file, bead in zip(bead_out_files, cv_set):
        q = qd.get_q(bead, ghts)
        d = qd.get_d(bead, ghts)
        nitems = len(bead.value) + 2
        data = [q.value * SQAMU, d.value * SQAMU] + list(bead.value)
        bead_file.write(('{:>12.4e}' * nitems + '\n').format(*data))
        bead_file.flush()


def convert_modes(modes):
    return {
        'b': modes['b'] / SQAMU,
        'a': modes['a'] / np.linalg.norm(modes['a'])
    }


def recover_modes(state, modes):
    state['modes'].update({
        'b': modes['b'] * SQAMU,
        'a': modes['a'] / np.linalg.norm(modes['a'])
    })


def convert_ghts(ghts):
    M = ghts['M'] * AMU
    Minv = np.linalg.inv(M)
    n = normal(ghts['n'] * SQAMU, Minv)
    return {
        'z': ghts['z'],
        'n': n,
        'M': M,
        'Minv': Minv,
        'Minvn': np.matmul(Minv, n)
    }


def recover_ghts (state, ghts):
    state['ghts'].update({
        'z': ghts['z'],
        'n': ghts['n'] / SQAMU,
        'M': ghts['M'] / AMU,
    })


def convert_params (params):
    return {
        'K': params['K'] / KCAL_MOL * AMU,
        'K_d': params['K_d'] / KCAL_MOL * AMU,
        'd_max': params['d_max'] / SQAMU,
        'gamma_b': params['gamma_b'] * AMU / PS,
        'gamma_a': params['gamma_a'] / SQAMU / PS,
        'gamma_delta': params['gamma_delta'] / PS,
        'gamma_z': params['gamma_z'] * AMU / PS,
        'gamma_d': params['gamma_d'] / PS,
        'gamma_n': params['gamma_n'] / PS,
        'gamma_M': params['gamma_M'] / PS,
        'move_modes': params['move_modes'],
        'move_ghts': params['move_ghts'],
        'fix_M': params['fix_M']
    }