"""
Chromosome dynamics algorithm selection.

Selects how the chromosome is advanced at each DNA hook and binds ``run`` to
that implementation, so the hook does not need to know which one is active.

Authors
-------
Alfia Parvez -- adapted from the Modularize_4DWCM_MinCell reference design
"""

from processes.SpatialDnaDynamics import updateChromosome, updateChromosomeDivision


class DNA_Dynamics:
    """
    Parameters
    ----------
    sim_properties : dict
        Simulation state and knobs.
    DNA_algorithm : str
        ``'BD'`` for Brownian dynamics in LAMMPS (requires a second GPU), or
        ``'lattice'`` for the Python surrogate.
    """

    def __init__(self, sim_properties: dict, DNA_algorithm: str = 'BD') -> None:

        self.sim_properties = sim_properties
        self.DNA_algorithm = DNA_algorithm

        if DNA_algorithm == 'BD':
            self.run = self._run_BD
            print("Chromosome: Brownian dynamics and SMC looping in LAMMPS")
        else:
            self.run = self._run_lattice
            print("Chromosome: surrogate algorithm manipulating lattice and gene "
                  "particles in Python")

    def _run_BD(self, time, lattice, sim_properties, region_dict, ribo_site_dict,
                updateRegions):
        """Advance the chromosome with Brownian dynamics via btree_chromo/LAMMPS."""

        if sim_properties['division_started'] and updateRegions:
            return updateChromosomeDivision(
                time, lattice, sim_properties, region_dict, ribo_site_dict)

        return updateChromosome(
            time, lattice, sim_properties, region_dict, ribo_site_dict, updateRegions)

    def _run_lattice(self, time, lattice, sim_properties, region_dict, ribo_site_dict,
                     updateRegions):
        """
        Advance the chromosome with the Python lattice surrogate.

        Imported lazily: the surrogate is not part of this repository yet, and a
        module-level import would make selecting ``'BD'`` fail too.
        """

        try:
            from processes.DNA_Lattice import DNAMoveLattice
        except ImportError as exc:
            raise NotImplementedError(
                "The lattice surrogate (processes/DNA_Lattice.DNAMoveLattice) is not "
                "available in this repository. Use -DNA BD, or port the surrogate "
                "from Modularize_4DWCM_MinCell."
            ) from exc

        return DNAMoveLattice(time, lattice, region_dict, ribo_site_dict, updateRegions)
