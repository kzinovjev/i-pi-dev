"""Defines dot products, distances, angles etc. for a CV space with a given a metric tensor"""
import numpy as np


def m_tensor(jacobian, masses):
    return np.matmul(jacobian * np.reciprocal(masses), jacobian.transpose())

class CVVector(object):
    """A wrapper class to facilitate working with geometry routines. For all binary operations
    with other CVVector instances the metric tensor of the first (left) CVVector is used"""
    def __init__(self, value, m):
        self.value = value
        self.m = m

    def __add__(self, other):
        if isinstance(other, CVVector):
            return CVVector(self.value + other.value, self.m)
        else:
            return CVVector(self.value + other, self.m)

    def __sub__(self, other):
        if isinstance(other, CVVector):
            return CVVector(self.value - other.value, self.m)
        else:
            return CVVector(self.value - other, self.m)

    def __mul__(self, other):
        if isinstance(other, CVVector):
            return CVVector(dot_product(self.value, other.value, self.m), self.m)
        elif isinstance(other, np.ndarray):
            return CVVector(np.matmul(self.value, other), self.m)
        else:
            return CVVector(self.value * other, self.m)

    def __len__(self):
        return len(self.value)

    def angle(self, a): return angle(self.value, a, self.m)

    def cosine(self, a): return cosine(self.value, a, self.m)

    def orthogonal(self, a): return orthogonal(self.value, a, self.m)

    def parallel(self, a): return parallel(self.value, a, self.m)

    def normal(self): return normal(self.value, self.m)

    def mod(self): return mod(self.value, self.m)


def angle(a, b, m):
    """Angle between vectors a and b"""
    return np.arccos(cosine(a, b, m))


def cosine(a, b, m):
    """Cosine of the angle between vectors a and b"""
    return dot_product(normal(a, m), normal(b, m), m)


def orthogonal(a, b, m):
    """Component of vector a orthogonal to vector b"""
    return a - parallel(a, b, m)


def parallel(a, b, m):
    """Component of vector a parallel to vector b"""
    return b * dot_product(a, normal(b, m), m)


def normal(a, m):
    """Unit vector collinear to vector a"""
    return a/mod(a, m)


def dist(a, b, m):
    """Distance between points a and b"""
    return mod(a - b, m)


def mod(a, m):
    """Module of vector a"""
    return np.sqrt(mod2(a, m))


def mod2(a, m):
    """Squared module of vector a"""
    return dot_product(a, a, m)


def dot_product(a, b, m):
    """Dot product between vectors a and b"""
    return a.dot(np.matmul(m, b))