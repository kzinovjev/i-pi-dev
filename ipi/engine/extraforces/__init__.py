import pkgutil
import importlib

modules = [name for _, name, _ in pkgutil.iter_modules(__path__)]
for module in modules:
    importlib.import_module(__package__+'.'+module)

forces = [eval(module).ipi_force for module in modules if hasattr(eval(module), 'ipi_force')]

force_cache = None
call_calc_forces = False


def calc(integrator):
    global force_cache, call_calc_forces

    #calc is called twice per velocity verlet step, so cache the forces for the next call
    call_calc_forces = not call_calc_forces
    if call_calc_forces:
        force_cache = calc_forces(integrator)
    return force_cache


def calc_forces(integrator):
    return reduce(lambda f, force: f + force(integrator), forces, 0)
