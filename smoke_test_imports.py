#!/usr/bin/env python3
"""
Import smoke test for the package layout.

Verifies that every first-party module imports, that process-global singletons
resolve to one object, and that the helper scripts MC_CME spawns by filename
still exist. Run from the repo root inside the simulation container:

    python smoke_test_imports.py

Exits non-zero on the first failure.
"""

import importlib
import os
import sys

MODULES = [
    # leaves first, then upward through the dependency graph
    "utility.LatticeFunctions", "utility.GIP_rates", "utility.Integrate",
    "utility.FileSaving",
    "processes.Diffusion", "processes.FreeDTS_functions",
    "processes.InitRdmeDna", "processes.Rxns_CME",
    "processes.Rxns_RDME", "processes.Rxns_ODE",
    "processes.ImportInitialConditions", "processes.RegionsAndComplexes",
    "processes.MC_CME", "processes.Communicate",
    "processes.RibosomesRDME", "processes.Growth", "processes.Division",
    "processes.SpatialDnaDynamics", "processes.MC_RDME_initialization",
    "modules.SIM_State", "modules.DNA_Dynamics", "modules.Metabolism",
    "processes.Hook",
    "restart.Restart_MC_RDME_initialization", "restart.Restart_Hook",
]

# Scripts MC_CME launches by path at runtime. An import test cannot catch a
# stale path here: _ensure_worker falls back to os.system when the worker
# script is missing, so the run stays correct but loses the persistent worker.
SPAWNED_SCRIPTS = ["processes/Run_CME.py", "processes/Run_CME_Worker.py"]

failures = []

print("importing modules")
for name in MODULES:
    try:
        importlib.import_module(name)
        print(f"  ok    {name}")
    except Exception as exc:
        print(f"  FAIL  {name}: {type(exc).__name__}: {exc}")
        failures.append(name)

print("\nchecking spawned scripts exist")
for rel in SPAWNED_SCRIPTS:
    if os.path.isfile(rel):
        print(f"  ok    {rel}")
    else:
        print(f"  FAIL  {rel} not found")
        failures.append(rel)

print("\nchecking process-global singletons")
if not failures:
    import processes.MC_CME as a
    import processes.MC_CME as b
    import processes.RibosomesRDME as r1
    import processes.RibosomesRDME as r2
    for label, x, y in (("MC_CME", a, b), ("RibosomesRDME", r1, r2)):
        if x is y:
            print(f"  ok    {label} is a single module object")
        else:
            print(f"  FAIL  {label} imported twice under different identities")
            failures.append(label)

print()
if failures:
    print(f"FAILED: {len(failures)} problem(s): {', '.join(failures)}")
    sys.exit(1)
print(f"PASSED: {len(MODULES)} modules imported, layout is consistent")
