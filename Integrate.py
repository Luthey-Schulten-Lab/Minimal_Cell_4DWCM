"""
Authors: Zane Thornburg, David Biancho

Integrate system of ODE reactions
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

    # Capture the numeric optimization-space parameter values BEFORE
    # prepareFunctor() runs. prepareFunctor() calls applyOpt(), which overwrites
    # every parameter's .val with the *string* "self.params[k]" used for code
    # generation -- after that point getOptSpace() no longer returns numbers.
    # The compiled functor's self.params[] array must be filled with these
    # parameter values (length == number of optimization-space parameters, i.e.
    # the highest self.params[k] index in calcFlux_c + 1), NOT the metabolite
    # initial concentrations from getInitVals() (length == number of tracked
    # metabolites). Using getInitVals() here sizes self.params to the wrong
    # length and makes calcFlux_c index off the end (the documented IndexError).
    _, initParamVals = model.getOptSpace()

    solvFunctor.prepareFunctor() #OG uncomment

    # Set verbosity to 0 for now, below uncomment OG
    rxnIdList = solvFunctor.buildCall(odeint=False, useJac=False, cythonBuild=True, functor=True, verbose=0)
    #rxnIdList = solvFunctor.buildCall(odeint=True, useJac=False, cythonBuild=False, functor=False, transpJac=False, verbose=0, noBuild=True)

    # Sets up the actual solver, with the optimization-space parameter values.
    solvFunctor = solvFunctor.functor( np.asarray(initParamVals, dtype=np.double) )

    return solvFunctor
#########################################################################################


#########################################################################################
def setSolverCached(model, cache):
    """
    Build & Cython-compile the ODE flux functor ONCE, then on later calls only
    re-instantiate the already-compiled module with the current parameter values.

    Why: odecell's buildCall(cythonBuild=True) regenerates the .pyx and shells out
    to gcc every single call (~95 s). It only had to recompile because the per-step
    enzyme concentrations were *baked* into the generated code. With "Enzyme" moved
    into the optimization space (Rxns_ODE: addParameter("Enzyme", ..., lb!=ub)),
    every per-step-varying quantity now flows through the compiled stepClass'
    self.params[] array, so the generated/compiled code is byte-for-byte identical
    across steps. We therefore compile once and afterwards just call
    stepClass(currentParams) -- O(milliseconds) -- instead of recompiling.

    Correctness: metabolite concentrations are ODE state (passed as y0 via
    getInitVals each step, not params); kcat/KM/medium are constants baked at the
    first compile and unchanged thereafter. The only thing that changes per step is
    the opt-space parameter vector, which we refill on every call. If the opt-space
    or metabolite *structure* changes (e.g. replication adds/removes reactions) the
    cache signature changes and we transparently recompile.

    Parameters:
        model  - a freshly built odecell model for the current step
        cache  - a mutable dict persisted across calls by the caller (the Hook)

    Returns:
        solver - a ready-to-integrate functor (callable as solver(t, y))
    """

    # Capture the numeric optimization-space values BEFORE any prepareFunctor()
    # mutation (prepareFunctor rewrites every opt param's .val to the string
    # "self.params[k]"). On reuse no mutation happens, so this stays numeric.
    optList, initParamVals = model.getOptSpace()

    # Structure signature: number of metabolites (= length/order of the y vector
    # the compiled code indexes) plus the ordered identity of every opt-space
    # parameter. Values are intentionally excluded -- only topology matters here.
    nMet = len(model.getInitVals())
    sig = (nMet, tuple((o.type, str(o.indx)) for o in optList))

    if cache.get('builder') is None or cache.get('sig') != sig:
        # (Re)compile: first call, or the reaction/metabolite structure changed.
        builder = odecell.solver.ModelSolver(model)
        builder.prepareFunctor()  # mutates THIS model's opt params -> self.params[k]
        builder.buildCall(odeint=False, useJac=False, cythonBuild=True,
                          functor=True, verbose=0)
        # buildCall binds builder.functor = compiled cythonCompiledFunctions.stepClass
        cache['builder'] = builder
        cache['sig'] = sig

    # Cheap: instantiate the compiled stepClass with the current numeric params.
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

    #y0 = model.getInitVals()
    #print("shape: ",len(y0))
    #y0 = np.asarray(y0,dtype=np.double)

    #modelOptSpace,initParamVals=model.getOptSpace()
    #dydt = np.zeros(len(y0))
    #tout, results, info = integrate_adaptive(f_wrap(solv,time+delt,y0,dydt),None,y0,time,time+delt,atol,rtol,dx0=ts,nsteps=10000)
    #tout, results, info = integrate_adaptive(f_wrap(solv,time,y0,dydt),None,y0,time,time+delt,atol,rtol,dx0=ts,nsteps=10000)

    #tout = 0.0
    #info = "place holder"

    #solv = solv.ModelSolver(model)
    #solv.buildCall(odeint=False, useJac=False, verbose=2)
    integrator = integrate.ode(solv)#, solv.calcJac)

    # TODO: Set which integrator to use (Appears vode is default), lsoda adaptive
#     lsodaBool = True
#     if (lsodaBool):
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