import ghts.io as io
import numpy as np

state = io.load_state("ghts.json")

if state is not None and state['stage']['name'] != 'committor':
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.Get_rank()
else:
    import mpi4py.rc
    mpi4py.rc.initialize = False
    mpi4py.rc.finalize = False
    rank = 0

from ghts.ghts import force

step = 0


def ipi_force(integrator):
    global step
    if state is None:
        return 0

    result = force(integrator.beads, integrator.cell, integrator.beads.m3[0],
                   integrator.ensemble.temp, integrator.dt, state)

    if rank == 0:
        step += 1
        if state['stage']['name'] != 'optimize':
            return result
        if step % state['output']['print_every'] == 0:
            print_ghts(state['modes'], state['ghts'], step)
        if step % state['output']['save_every'] == 0:
            save_state(state, step)
    return result


def save_state(state, step):
    io.dump_state(state, str(step) + ".json")


def print_ghts(modes, ghts, step):
    with open('ghts.out', 'a') as f:
        nitems = len(ghts['z']) * 2 + len(modes['a']) * 2 + 2
        sigma0 = - np.dot(modes['a'], modes['b'])
        data = [step, sigma0] + \
               list(modes['b']) + list(modes['a']) + list(ghts['z']) + list(ghts['n'])
        f.write(('{:>12.4e}' * nitems + '\n').format(*data))
    if state['ghts'].get('M') is not None:
        with open('M.out', 'w') as f:
            for row in ghts['M']:
                f.write(('{:>12.4e}' * len(ghts['z']) + '\n').format(*row))


#Initial print
if state is not None and rank == 0 and state['stage']['name'] == 'optimize':
    print_ghts(state['modes'], state['ghts'], step)
    save_state(state, step)