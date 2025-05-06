# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================


import numpy as np

import os
import importlib

import pickle
import glob
import re

class WCM_traj:

    def __init__(self,traj_file):

        self.traj_file = traj_file

        with open(self.traj_file,'r') as f:

            # determine the number of columns
            ncols = len(f.readline().split(','))

            # read the set of species and fluxes
            f.seek(0)
            
            species_and_rxns = np.loadtxt(f,
                                          usecols=0,
                                          skiprows=1,
                                          dtype=np.str_,
                                          delimiter=',')

            rxn_start = species_and_rxns.shape[0]
            for i in range(species_and_rxns.shape[0]):
                if species_and_rxns[i][:2] == 'F_':
                    rxn_start = i
                    break

            self.species = species_and_rxns[:rxn_start]
            self.N_species = self.species.shape[0]

            self.rxns = species_and_rxns[rxn_start:]
            self.N_rxns = self.rxns.shape[0]

            for i_rxn in range(self.N_rxns):
                self.rxns[i_rxn] = self.rxns[i_rxn][2:]

            # read the set of timesteps
            f.seek(0)
            self.t = np.loadtxt(f,
                                max_rows=1,
                                usecols=range(1,ncols),
                                dtype=np.float32,
                                delimiter=',')

            self.Nt = self.t.shape[0]

            # read the particle counts
            f.seek(0)
            self.x = np.loadtxt(f,
                                usecols=range(1,ncols),
                                skiprows=1,
                                max_rows=self.N_species,
                                dtype=np.float32,
                                delimiter=',')

            print(self.traj_file)
            
            print('Species Array:')
            print(self.x.shape)

            # read the reaction fluxes
            f.seek(0)
            self.fx = np.loadtxt(f,
                                 usecols=range(1,ncols),
                                 skiprows=1+self.N_species,
                                 max_rows=self.N_rxns,
                                 dtype=np.float32,
                                 delimiter=',')

            self.fx[:,0] = self.fx[:,1]

            print('Flux Array:')
            print(self.fx.shape)

        # set the species map
        self.species_map = dict()
        for i_specie in range(self.N_species):
            self.species_map[self.species[i_specie]] = i_specie

        # set the rxns map
        self.rxns_map = dict()
        for i_rxn in range(self.N_rxns):
            self.rxns_map[self.rxns[i_rxn]] = i_rxn

        return

    def get_species(self):

        return self.species

    def get_species_map(self):

        return self.species_map

    def get_rxns(self):

        return self.rxns

    def get_rxns_map(self):

        return self.rxns_map

    def get_t(self):

        return self.t

    def get_x(self):

        return self.x

    def get_fx(self):

        return self.fx
