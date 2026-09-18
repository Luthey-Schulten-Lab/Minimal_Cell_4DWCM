"""
Metabolism algorithm selection.

Selects how metabolism is advanced at each metabolism hook and binds ``run`` to
that implementation, so the hook does not need to know which one is active.

Authors
-------
Alfia Parvez -- adapted from the Modularize_4DWCM_MinCell reference design
"""

import processes.Rxns_ODE as ODE
import processes.Communicate as communicate
import utility.Integrate as integrate


class Metabolism:
    """
    Parameters
    ----------
    sim_properties : dict
        Simulation state and knobs.
    metabolism_algorithm : str
        ``'ODE'`` to integrate the metabolic network, or ``'skip'`` to hold
        metabolite concentrations fixed.
    """

    def __init__(self, sim_properties: dict, metabolism_algorithm: str = 'ODE') -> None:

        self.sim_properties = sim_properties
        self.metabolism_algorithm = metabolism_algorithm

        # Cython solver state. Owned here rather than by the hook: the build-once
        # cache is keyed on the model and only this class integrates the network.
        self._use_cython = True
        self._cython_failed = False
        self._solver_cache = {}

        # The counts/flux writer runs on its own cadence, not the metabolism
        # hook's, so it reports whichever integration ran most recently rather
        # than one of its own. Retain that here for it to read back.
        self.last_results = None
        self.last_model = None
        self.last_solver = None

        if metabolism_algorithm == 'ODE':
            self.run = self._run_ODE
            print("Metabolism: ordinary differential equations")
        else:
            self.run = self._skip_ODE
            print("Metabolism: skipped, concentrations held at their initial values")

    def _run_ODE(self, time):
        """
        Integrate the metabolic network for one hook interval.

        Returns
        -------
        odeResults : numpy.ndarray
            Trajectory; the last row is the updated state.
        model : odecell model
            The model the trajectory was produced with, needed to map results
            back onto species counts.
        solver : odecell solver
            Needed by the counts/flux writer to recover per-reaction fluxes.
        """

        print('Initializing ODE simulation')
        model = ODE.initModel(self.sim_properties)
        print('Initialized ODE simulation')

        initVals = integrate.getInitVals(model)

        # Cython path, with a one-way fallback to noCython for the rest of the
        # run so a failed build does not get retried every hook.
        if self._use_cython and not self._cython_failed:
            try:
                solver = integrate.setSolverCached(model, self._solver_cache)
            except Exception as exc:
                print(f"ODE: Cython solver build failed ({exc}); falling back to "
                      "noCython for the remainder of the run.")
                self._cython_failed = True
                self._solver_cache = {}
                solver = integrate.noCythonSetSolver(model)
        else:
            solver = integrate.noCythonSetSolver(model)

        return integrate.runODE(initVals, solver, model), model, solver

    def _skip_ODE(self, time):
        """
        Hold metabolism fixed.

        Not available in this repository. The reference implementation reports
        counts instead of concentrations, which needs the four-argument
        ``updateCountsODE(sim_properties, abundances, MetDict, flag)`` and a
        populated ``sim_properties['fluxes']`` table; this repository has the
        three-argument form and never fills ``fluxes`` in. Wiring the skip path
        up means porting both, so it fails loudly rather than silently writing
        the wrong species.
        """

        raise NotImplementedError(
            "Metabolism algorithm 'skip' is not available in this repository: it "
            "requires the count-based updateCountsODE and flux bookkeeping from "
            "Modularize_4DWCM_MinCell. Use -MB ODE.")

    def update_metabolism(self, time) -> None:
        """Advance metabolism and write the new counts into ``sim_properties``."""

        odeResults, model, solver = self.run(time)

        communicate.updateCountsODE(self.sim_properties, odeResults, model)

        self.last_results = odeResults
        self.last_model = model
        self.last_solver = solver

        return None
