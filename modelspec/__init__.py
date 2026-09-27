"""ModelSpec: one read of every 4DWCM input, shared by the model browser, the perturbation checker and (next) the simulator.

build_spec()      parse input_data/ into a JSON-serialisable dict (genes, species, reactions, processes, code references)
perturbation      the YAML perturbation schema (knockouts, knockdowns, initial-condition and parameter edits): load, validate, resolve
impact            static first-order knockout impact over the spec (which reactions and processes lose an enzyme or subunit)

Command line: python -m modelspec {build,check,serve} -h
"""

SCHEMA_VERSION = 1
