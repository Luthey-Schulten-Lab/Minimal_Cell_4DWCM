"""
Container for the time-dependent simulation state.

Placeholder matching the Modularize_4DWCM_MinCell reference. The 4DWCM
currently carries this state in the ``sim_properties`` dictionary; this class
is where it is intended to move.

Exported state:
    counts / concentrations
    fluxes

Intermediate state:
    DNA coordinates
    replication states
"""


class SIM_State:
    def __init__(self) -> None:
        pass
