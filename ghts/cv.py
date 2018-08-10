import numpy as np
import molmod as mm
from collections import namedtuple
from cv_geometry import m_tensor

CV = namedtuple('CV', ['value', 'gradient'])
CVSet = namedtuple('CVSet', ['value', 'jacobian', 'm'])

def get_xyz(q, index):
    return q[index*3-3:index*3]


def set_xyz(array, index, xyz):
    array[index*3-3:index*3] = xyz


def get_full_gradient(q, atoms, gradient):
    result = np.zeros(q.shape, q.dtype)
    for atom, atom_gradient in zip(atoms, gradient):
        set_xyz(result, atom, atom_gradient)
    return result


def eval_cv(q, atoms, callback):
    value, gradient = callback([get_xyz(q, atom) for atom in atoms], 1)
    return CV(value, get_full_gradient(q, atoms, gradient))


def eval_dist(q, atoms): return eval_cv(q, atoms, mm.bond_length)


def eval_transfer(q, atoms):
    dist1 = eval_dist(q, atoms[:2])
    dist2 = eval_dist(q, atoms[1:])
    return CV(dist1.value-dist2.value, dist1.gradient-dist2.gradient)


def eval_angle(q, atoms): return eval_cv(q, atoms, mm.bend_angle)


def eval_dihedral(q, atoms): return eval_cv(q, atoms, mm.dihed_angle)


def eval_pplane(q, atoms):
    value, gradient = mm.opbend_dist([get_xyz(q, atom) for atom in atoms], 1)
    sign = np.sign(mm.opbend_angle([get_xyz(q, atom) for atom in atoms]))[0]
    return CV(value*sign, get_full_gradient(q, atoms, gradient*sign))


def eval_x(q, atoms):
    return CV(get_xyz(q, atoms[0])[0],
              get_full_gradient(q, atoms, np.array([[1, 0, 0],])))


CV_KIND_DICT = {
    'distance': eval_dist,
    'transfer': eval_transfer,
    'angle': eval_angle,
    'dihedral': eval_dihedral,
    'pplane': eval_pplane,
    'x': eval_x
}


def get_cv(q, cv_def):
    result = CV_KIND_DICT.get(cv_def['kind'])(q, cv_def['atoms'])
    power = cv_def.get('power')
    if power is None:
        return result
    return CV(
        value=np.power(result.value, power),
        gradient=power * np.power(result.value, power - 1) * result.gradient
    )


def get_cv_set(q, cv_def_list, masses):
    values, gradients = zip(*[get_cv(q, cv) for cv in cv_def_list])
    jacobian = np.array(gradients)
    return CVSet(
        value=np.array(values),
        jacobian=jacobian,
        m=m_tensor(jacobian, masses)
    )


def get_cv_lincomb(cv_set, coef, coef0):
    value = np.dot(cv_set.value, coef) + coef0
    gradient = np.matmul(coef, cv_set.jacobian)
    return CV(value, gradient)