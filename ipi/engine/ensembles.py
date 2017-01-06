"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time

import numpy as np

from ipi.utils.depend import *
from ipi.utils.softexit import softexit
from ipi.utils.io import read_file
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import unit_to_internal, Constants
from ipi.engine.thermostats import *
from ipi.engine.barostats import *


__all__ = ['Ensemble', 'ensemble_swap']

### IMPORTANT - THIS MUST BE KEPT UP-TO-DATE WHEN THE ENSEMBLE CLASS IS CHANGED
def ensemble_swap(ens1, ens2):
    """ Swaps the definitions of the two ensembles, by
    exchanging all of the inner properties. """

    if ens1.temp != ens2.temp :
        swp = ens1.temp; ens1.temp = ens2.temp; ens2.temp = swp
        

class Ensemble(dobject):
    """Base ensemble class.

    Defines the thermodynamic state of the system.

    Depend objects:
        temp: The system's temperature.
        pext: The systems's pressure
        stressext: The system's stress tensor
    """

    def __init__(self, eens=0.0, econs=0.0, temp=None, pext=None, stressext=None):
        """Initialises Ensemble.

        Args:
            temp: The temperature.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """

        dset(self, "temp", depend_value(name='temp'))
        if temp is not None:
            self.temp = temp
        else:
            self.temp = 0.0

        dset(self, "stressext", depend_array(name='stressext', value=np.zeros((3,3), float)))
        if stressext is not None:
            self.stressext = np.reshape(np.asarray(stressext), (3,3))
        else:
            self.stressext = 0.0

        dset(self, "pext", depend_value(name='pext'))
        if pext is not None:
            self.pext = pext
        else:
            self.pext = 0.0

        dset(self, "eens", depend_value(name='eens'))
        if eens is not None:
            self.eens = eens
        else:
            self.eens = 0.0
                    

    def bind(self, beads, nm, cell, bforce, bbias, elist=[]):
        self.beads = beads
        self.cell = cell
        self.forces = bforce
        self.bias = bbias
        self.nm = nm
        dset(self, "econs", depend_value(name='econs', func=self.get_econs))        
        # dependencies of the conserved quantity
        dget(self, "econs").add_dependency(dget(self.nm, "kin"))
        dget(self, "econs").add_dependency(dget(self.forces, "pot"))
        dget(self, "econs").add_dependency(dget(self.bias, "pot"))
        dget(self, "econs").add_dependency(dget(self.beads, "vpath"))
        dget(self, "econs").add_dependency(dget(self, "eens"))

        self._elist = []

        for e in elist:
            self.add_econs(e)
            
        dset(self, "lpens", depend_value(name='lpens', func=self.get_lpens, 
                dependencies=[ dget(self,"temp"), dget(self,"pext"), dget(self,"stressext") ] ))        
        dget(self, "lpens").add_dependency(dget(self.nm, "kin"))
        dget(self, "lpens").add_dependency(dget(self.forces, "pot"))
        dget(self, "lpens").add_dependency(dget(self.bias, "pot"))
        dget(self, "lpens").add_dependency(dget(self.beads, "vpath"))

    def add_econs(self, e):
        self._elist.append(e)
        dget(self, "econs").add_dependency(e)

    def get_econs(self):
        """Calculates the conserved energy quantity for constant energy
        ensembles.
        """
        eham = self.beads.vpath*self.nm.omegan2 + self.nm.kin + self.forces.pot
        eham += self.bias.pot   # bias
        for e in self._elist:
            eham += e.get()

        return eham + self.eens
        
    def get_lpens(self):
        """Returns the ensemble probability (modulo the partition function) 
        for the ensemble. 
        """
        
        # TODO include the stress 
        lpens = -(self.pext*self.cell.V+self.forces.pot+self.bias.pot+self.nm.kin+self.beads.vpath)/(Constants.kb*self.temp)
        print "Computing ensemble stuff ", lpens
        return lpens
