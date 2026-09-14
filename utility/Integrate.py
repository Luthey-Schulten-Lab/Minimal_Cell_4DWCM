"""
Integrate the metabolic ODE system.

Authors
-------
Alfia Parvez — opt-space param capture before prepareFunctor; ``setSolverCached``
Zane Thornburg, David Biancho — original odecell / Cython solver setup
"""

from pycvodes import integrate_predefined
from pycvodes import integrate_adaptive
from scipy import integrate
import odecell
import numpy as np


### Constants
step = 0.1 # s
atol = 1e-6 # tolerance
rtol = 1e-6 # tolerance

#########################################################################################
def setSolver(model):
    """
    Set the solver for the model

    Parameters:
         
        (odecell.model) - the model object

    Returns:

        solvFunctor - a functor for the solver

    """

    ## We are NOT building for odeint (gives us more room to chose between CVODES and SciPy-ODE).
    ## We are NOT using a jacobian, since we do not have the partial derivatives for all rate forms.
    ## We are building with Cython for speed, this is a big model.

    # Builds the solver using a *Functor* interface
    solvFunctor = odecell.solver.ModelSolver(model)

    # Numeric opt-space params before prepareFunctor() (which rewrites .val to
    # "self.params[k]" strings). Do not use getInitVals() here — that sizes
    # self.params to metabolite count and causes IndexError in calcFlux_c.
    _, initParamVals = model.getOptSpace()

    solvFunctor.prepareFunctor() #OG uncomment

    # Set verbosity to 0 for now, below uncomment OG
    rxnIdList = solvFunctor.buildCall(odeint=False, useJac=False, cythonBuild=True, functor=True, verbose=0)
    #rxnIdList = solvFunctor.buildCall(odeint=True, useJac=False, cythonBuild=False, functor=False, transpJac=False, verbose=0, noBuild=True)

    solvFunctor = solvFunctor.functor( np.asarray(initParamVals, dtype=np.double) )

    return solvFunctor
#########################################################################################


#########################################################################################
def setSolverCached(model, cache):
    """
    Cython-compile the ODE flux functor once; later calls only refill params.

    Requires Enzyme (and other per-step quantities) in the opt-space so the
    generated code is stable across hooks. ``cache`` is a Hook-owned dict;
    recompiles if metabolite/opt-space topology (``sig``) changes.
    """
    # Numeric opt-space values before any prepareFunctor() mutation.
    optList, initParamVals = model.getOptSpace()

    # Topology only: metabolite count + ordered opt-space identities.
    nMet = len(model.getInitVals())
    sig = (nMet, tuple((o.type, str(o.indx)) for o in optList))

    if cache.get('builder') is None or cache.get('sig') != sig:
        builder = odecell.solver.ModelSolver(model)
        builder.prepareFunctor()
        builder.buildCall(odeint=False, useJac=False, cythonBuild=True,
                          functor=True, verbose=0)
        cache['builder'] = builder
        cache['sig'] = sig

    solver = cache['builder'].functor( np.asarray(initParamVals, dtype=np.double) )

    return solver
#########################################################################################


#########################################################################################
### NOTE: Have to get a callable f(y,t) for scipy.ode without creating the functor
def noCythonSetSolver(model):
    """
    Set the solver without compiling via Cython

    Parameters:

    model (odecell Model object): The model object

    Returns:

    solver (odecell Solver object): The Solver object, to solve the system of ODEs representing metabolic reactions
    """

    # Construct a Model Solver Object
    solver = odecell.solver.ModelSolver(model)

    rxnIdList =solver.buildCall(verbose=0, useJac=False, transpJac=0, nocheck=False, odeint=False, cythonBuild=False, functor=False, noBuild=True)

    return solver 
#########################################################################################


#########################################################################################
def f_wrap(solv, t, y, dydt):

    #solv = setSolver(model)
    dydt[:] = solv(0,np.asarray(y))[:]
#########################################################################################


#########################################################################################
def getInitVals(model):
    y0=model.getInitVals()
#     print(y0)
    return y0
#########################################################################################

# def noCythonRunODE():
#     """
#     Run the ODE Model without compiling via Cython

#     Parameters:

#     Returns:
#     None
#     """
#     return 0

#########################################################################################
def runODE(y0, solv, model):
    """
    Run the ODE Model after getting initial conditions

    Parameters:

    y0 (seems non-necessary) - can remove
    time (float): the current hybrid simulation time
    delt (float): the communication timestep between stochastic and deterministic simulation
    ts (float): the timestep for the adaptive ODE Solver
    solv (odecell Solver object): The solver object, with call built
    model (odecell Model object): The model object

    Returns:

    results (np.array): the array containing ODE Simulation Results (Maybe only the last time should be passed?)
    """

    integrator = integrate.ode(solv)#, solv.calcJac)

    integrator.set_integrator("lsoda")

    integrator.set_initial_value(model.getInitVals())

    ### With fixed timestepping
    step = 0.01
    totalTime = 1.0
    results = np.empty((0,len(model.getInitVals())), float)

    while integrator.successful() and integrator.t < totalTime:
        currConcentration = integrator.integrate(integrator.t + step)
#         print(integrator.t)
        # Silence integrator output for now
        #print(integrator.t, currConcentration)
        results = np.append(results, [np.asarray(currConcentration)], axis=0 )

#     results = results[-1,:]

    # Return only the last timestep for the results, this is all thats really needed
    return results
#########################################################################################