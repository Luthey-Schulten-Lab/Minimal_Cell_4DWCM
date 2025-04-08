import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.patches as patches
from matplotlib.lines import Line2D

from mpl_toolkits.axes_grid1 import make_axes_locatable

import os
import importlib

import pickle
import glob
import re

import locus_tag_filter as ltf
importlib.reload(ltf)

import WCM_traj as WCMt
importlib.reload(WCMt)


# use LaTeX fonts in the plot
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage[helvet]{sfmath}'

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

mm = 1/25.4
l_w_ratio = 1/1.618

class WCM_ensemble:

    def __init__(self):

        # input information for trajectories
        self.in_dir = None
        self.in_label = None
        self.reps = None

        # trajectory replicates
        self.N_reps = None        
        self.traj_files = None
        self.trajs = None

        # system information
        self.t = None
        self.Nt = None
        self.species = None
        self.N_species = None
        self.species_map = None
        self.rxns = None
        self.N_rxns = None
        self.rxns_map = None

        # merged trajectories
        self.x = None
        self.fx = None
        self.volumes = None
        self.conc_factors = None

        # locus tag info
        self.lti = ltf.locus_tag_info()

        return

    def write_merged_ensemble(self, out_dir, out_label):

        out_ensemble = dict()

        out_ensemble['N_reps'] = self.N_reps
        out_ensemble['t'] = self.t
        out_ensemble['Nt'] = self.Nt
        out_ensemble['species'] = self.species
        out_ensemble['N_species'] = self.N_species
        out_ensemble['species_map'] = self.species_map
        out_ensemble['rxns'] = self.rxns
        out_ensemble['N_rxns'] = self.N_rxns
        out_ensemble['rxns_map'] = self.rxns_map
        out_ensemble['x'] = self.x
        out_ensemble['fx'] = self.fx
        out_ensemble['volumes'] = self.volumes
        out_ensemble['conc_factors'] = self.conc_factors

        out_file = out_dir + out_label + '_WCMensemble.pkl'

        with open(out_file, 'wb') as f:
            pickle.dump(out_ensemble,f)

        return

    def read_merged_ensemble(self, in_dir, in_label):

        in_file = in_dir + in_label + '_WCMensemble.pkl'

        with open(in_file, 'rb') as f:
            in_ensemble = pickle.load(f)

        self.N_reps = in_ensemble['N_reps']
        self.t = in_ensemble['t']
        self.Nt = in_ensemble['Nt']
        self.species = in_ensemble['species']
        self.N_species = in_ensemble['N_species']
        self.species_map = in_ensemble['species_map']
        self.rxns = in_ensemble['rxns']
        self.N_rxns = in_ensemble['N_rxns']
        self.rxns_map = in_ensemble['rxns_map']
        self.x = in_ensemble['x']
        self.fx = in_ensemble['fx']
        self.volumes = in_ensemble['volumes']
        self.conc_factors = in_ensemble['conc_factors']
        
        return

    def set_traj_files(self, in_dir, in_label, reps):

        self.in_dir = in_dir
        self.in_label = in_label
        self.reps = reps
        self.N_reps = reps.shape[0]

        self.traj_files = []

        for i_rep in range(self.N_reps):

            traj_file = self.in_dir + self.in_label + \
                '.{:d}'.format(self.reps[i_rep]) + '.csv'

            self.traj_files.append(traj_file)

        return

    def get_traj_files(self):

        return self.traj_files

    def load_trajs(self):

        self.trajs = []

        for traj_file in self.traj_files:

            self.trajs.append(WCMt.WCM_traj(traj_file))

        return

    def merge_trajs(self):

        # determine the common set of timesteps and species
        for i_rep in range(self.N_reps):

            traj = self.trajs[i_rep]

            if i_rep == 0:
                
                self.t = traj.get_t()
                self.species = traj.get_species()
                self.N_species = self.species.shape[0]
                self.species_map = traj.get_species_map()
                self.rxns = traj.get_rxns()
                self.N_rxns = self.rxns.shape[0]
                self.rxns_map = traj.get_rxns_map()
                
            else:

                self.t = np.intersect1d(self.t,traj.get_t())
                

        # size the new array
        self.Nt = self.t.shape[0]
        self.x = np.zeros((self.N_species,self.Nt,self.N_reps),
                          dtype=np.int32)
        self.fx = np.zeros((self.N_rxns,self.Nt,self.N_reps),
                          dtype=np.float32)

        # set the elements of the new array
        for i_rep in range(self.N_reps):

            t_rep = self.trajs[i_rep].get_t()
            x_rep = self.trajs[i_rep].get_x()
            fx_rep = self.trajs[i_rep].get_fx()

            t_tot_rep, t_tot_idx, t_rep_idx = np.intersect1d(self.t,
                                                             t_rep,
                                                             return_indices=True)
            
            self.x[:,:,i_rep] = x_rep[:,t_rep_idx]
            self.fx[:,:,i_rep] = fx_rep[:,t_rep_idx]

        print(self.x.shape)
        print(self.fx.shape)

        self.calc_volumes_and_conc_factors()

        return

    def calc_volumes_and_conc_factors(self):

        self.volumes = self.get_specie_trace('Volume')
        self.volumes = self.volumes.astype(np.float32)

        # convert to liters
        self.volumes = (1.0E-21)*self.volumes

        # particles/mol
        N_avogadro = 6.0221409E+23

        # conversion factor to molar
        self.conc_factors = np.reciprocal(N_avogadro*self.volumes)

        return

    def get_t(self):

        return self.t

    def get_specie_trace(self,specie):
        
        return self.x[self.species_map[specie],:,:]

    def get_species_traces(self,species):

        x_subset = np.zeros((len(species),self.Nt,self.N_reps),
                            dtype=np.int32)

        for i_specie in range(len(species)):

            x_subset[i_specie] = self.get_specie_trace(species[i_specie])

        return x_subset

    def get_avg_species_traces(self,species):

        x_subset = self.get_species_traces(species)

        return np.mean(x_subset,axis=2)

    def get_conc_trace(self,specie):

        return self.conc_factors*self.get_specie_trace(specie)

    def get_conc_traces(self,species):

        c_subset = self.conc_factors*self.get_species_traces(species)

        return c_subset

    def get_avg_conc_traces(self,species):

        c_subset = self.get_conc_traces(species)

        return np.mean(c_subset,axis=2)

    def get_rxn_trace(self,rxn):

        return self.fx[self.rxns_map[rxn],:,:]

    def get_rxn_traces(self,rxns):

        fx_subset = np.zeros((len(rxns),self.Nt,self.N_reps),
                             dtype=np.float32)

        for i_rxn in range(len(rxns)):

            fx_subset[i_rxn] = self.get_rxn_trace(rxns[i_rxn])

        return fx_subset

    def get_avg_rxn_traces(self,rxns):

        fx_subset = self.get_rxns_traces(rxns)

        return np.mean(fx_subset,axis=2)

    def get_weightedsum_counts(self,species,
                               weights):
        
        y = self.get_species_traces(species)

        y = np.tensordot(weights,y,axes=1)
        
        return y

    def get_weightedsum_rxns(self,rxns,
                               weights):
        
        y = self.get_rxn_traces(rxns)

        y = np.tensordot(weights,y,axes=1)
        
        return y

    def get_group_traces(self,group_labels,
                         grouped_species):

        N_groups = len(group_labels)
        
        y = np.zeros((N_groups,self.Nt,self.N_reps),
                      dtype=np.float32)

        for i_group in range(N_groups):

            temp_label = group_labels[i_group]

            y[i_group,:,:] = self.get_weightedsum_counts(grouped_species[temp_label]['species'],
                                                         grouped_species[temp_label]['weights'])

            ref_specie = grouped_species[temp_label]['ref_specie']
            if ref_specie != None:

                yref = self.get_specie_trace(ref_specie)[0,:]

                y[i_group,:,:] += yref

        return y

    def get_group_rxn_traces(self,group_labels,
                             grouped_rxns):

        N_groups = len(group_labels)
        
        y = np.zeros((N_groups,self.Nt,self.N_reps),
                      dtype=np.float32)

        for i_group in range(N_groups):

            temp_label = group_labels[i_group]

            y[i_group,:,:] = self.get_weightedsum_rxns(grouped_rxns[temp_label]['rxns'],
                                                       grouped_rxns[temp_label]['weights'])

        return y

    def get_conc_group_count_group_rxn_traces(self,
                                              conc_species,
                                              group_labels_count,
                                              grouped_species,
                                              group_labels_rxn,
                                              grouped_rxns):


        merged_labels = conc_species + group_labels_count + group_labels_rxn

        N_concs = len(conc_species)
        N_group_counts = len(group_labels_count)
        N_group_rxns = len(group_labels_rxn)
        
        N_merged = len(merged_labels)

        y = np.zeros((N_merged,self.Nt,self.N_reps),
                     dtype=np.float32)

        # add concentrations to merged array
        y_temp = self.get_conc_traces(conc_species)
        y[0:N_concs,:,:] = y_temp

        # add grouped counts to merged array
        y_temp = self.get_group_traces(group_labels_count,
                                       grouped_species)
        y[N_concs:N_concs+N_group_counts,:,:] = y_temp.astype(np.float32)

        # add grouped reaction fluxes to merged arrayn
        y_temp = self.get_group_rxn_traces(group_labels_rxn,
                                           grouped_rxns)
        y[N_concs+N_group_counts:,:,:] = y_temp.astype(np.float32)

        return y, merged_labels

    def get_transcript_traces(self,lt):

        group_labels = []
        grouped_species = dict()

        for i in range(len(lt)):

            temp_label = 'transcript_' + lt[i]

            group_labels.append(temp_label)

            grouped_species[temp_label] = dict()

            grouped_species[temp_label]['ref_specie'] = 'R_' + lt[i]
            grouped_species[temp_label]['species'] = []
            grouped_species[temp_label]['weights'] = []

            grouped_species[temp_label]['species'].append('RPM_'+lt[i])
            grouped_species[temp_label]['weights'].append(1)

            grouped_species[temp_label]['species'].append('DM_'+lt[i])
            grouped_species[temp_label]['weights'].append(-1)

            
        y = self.get_group_traces(group_labels,
                       grouped_species)

        return group_labels, grouped_species, y

    def get_transcript_traces_nondeg(self,lt):

        group_labels = []
        grouped_species = dict()

        for i in range(len(lt)):

            temp_label = 'transcript_' + lt[i]

            group_labels.append(temp_label)

            grouped_species[temp_label] = dict()

            grouped_species[temp_label]['ref_specie'] = 'R_' + lt[i]
            grouped_species[temp_label]['species'] = []
            grouped_species[temp_label]['weights'] = []

            grouped_species[temp_label]['species'].append('RPM_'+lt[i])
            grouped_species[temp_label]['weights'].append(1)

            
        y = self.get_group_traces(group_labels,
                       grouped_species)

        return group_labels, grouped_species, y

    def get_protein_traces(self,lt):

        group_labels = []
        grouped_species = dict()

        for i in range(len(lt)):

            temp_label = 'protein_' + lt[i]

            group_labels.append(temp_label)

            grouped_species[temp_label] = dict()

            grouped_species[temp_label]['ref_specie'] = 'P_' + lt[i]
            grouped_species[temp_label]['species'] = []
            grouped_species[temp_label]['weights'] = []

            grouped_species[temp_label]['species'].append('PM_'+lt[i])
            grouped_species[temp_label]['weights'].append(1)

            
        y = self.get_group_traces(group_labels,
                       grouped_species)

        return group_labels, grouped_species, y

    def get_species_list(self):

        return self.species.tolist()

    def get_metabolite_species(self):

        metabolite_species = []

        for specie in self.get_species_list():

            if (specie.startswith('M_')) and ('trna' not in specie):

                metabolite_species.append(specie)

        return metabolite_species

    def get_nucleobase_species(self):

        nucleobase_species = []

        nb_info = dict()

        nucleobases = ['a','g','c','t','u']
        p_count = ['m','d','t']

        for nb in nucleobases:

            nb_info[nb] = dict()

        nb_info['a']['import'] = ['adn','dad_2']
        nb_info['g']['import'] = ['gsn','dgsn']
        nb_info['c']['import'] = ['cytd','dcyt']
        nb_info['t']['import'] = ['thymd']
        nb_info['u']['import'] = ['uri','duri']

        nb_info['a']['p_d_truth'] = np.array([[1,1,1],
                                              [1,1,1]],
                                             dtype=np.int32)
        nb_info['g']['p_d_truth'] = np.array([[1,1,1],
                                              [1,1,1]],
                                             dtype=np.int32)
        nb_info['c']['p_d_truth'] = np.array([[1,1,1],
                                              [1,1,1]],
                                             dtype=np.int32)
        nb_info['t']['p_d_truth'] = np.array([[0,0,0],
                                              [1,1,1]],
                                             dtype=np.int32)
        nb_info['u']['p_d_truth'] = np.array([[1,1,1],
                                              [1,1,1]],
                                             dtype=np.int32)

        specie_prefix = 'M_'
        specie_suffix = '_c'
        
        for nb in nucleobases:

            for s in nb_info[nb]['import']:

                temp_specie = specie_prefix + s + specie_suffix

                nucleobase_species.append(temp_specie)

            for i in range(2):

                ds = ''
                if i > 0:
                    ds = 'd'
                    
                
                for j in range(3):

                    if nb_info[nb]['p_d_truth'][i,j] == 1:

                        s = ds + nb + p_count[j] + 'p'
                        temp_specie = specie_prefix + s + specie_suffix
                        nucleobase_species.append(temp_specie)

                        
        return nucleobase_species

    def get_AA_species(self):

        AA_species = []

        for specie in self.get_species_list():

            if ('__L' in specie):

                AA_species.append(specie)

        return AA_species

    def get_tRNA_species(self):

        tRNA_species = []

        for specie in self.get_species_list():

            if ('trna' in specie):

                tRNA_species.append(specie)

        return tRNA_species

    def get_single_gene_species(self,locus_tag):

        single_gene_species = []

        for specie in self.get_species_list():

            if (('_' + str(locus_tag)) in specie):

                single_gene_species.append(specie)

        return single_gene_species

    def get_rxns_list(self):

        return self.rxns.tolist()

    def get_all_locus_tags(self):

        gene_species = self.get_regex_species_list(r"(\bG_\w+)")

        locus_tags = []

        for gene_specie in gene_species:

            locus_tags.append(gene_specie[2:6])

        locus_tags = np.array(locus_tags,dtype=np.str_)
        locus_tags = np.unique(locus_tags)
        locus_tags = locus_tags.tolist()

        return locus_tags

    def get_filtered_locus_tags(self,info_file):

        self.lti.load_info_file(self.get_all_locus_tags(),
                                info_file)

        return

    def get_loci_lengths(self,loci):

        return self.lti.get_lengths(loci)

    def get_tagged_loci(self,tag):

        return self.lti.get_tagged_loci(tag)

    def get_regex_species_list(self,species_regex):

        species_str = ' '.join(self.get_species_list())

        return re.findall(str(species_regex),species_str)

    def choose_colormap(self,N_items):

        if N_items <= 10:
            cmap = colormaps.get_cmap('tab10')
            c_space = np.arange(0,N_items,dtype=np.int32)
        else:
            cmap = colormaps.get_cmap('plasma')
            c_space_lower_lim = 0.0
            c_space_upper_lim = 0.9
            c_space = np.linspace(c_space_lower_lim,
                                  c_space_upper_lim,
                                  N_items,
                                  dtype.np.float32)
        return cmap, c_space

    def prepare_individual_plot_fig(self,y,ym,specie,rep_coloring):

        fig_size = [87,87*l_w_ratio]

        fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

        ax = plt.gca()

        ax.set_xlabel(r'$t$ - cell cycle time [s]',
                      fontsize=7,
                      labelpad=1.5)

        ax.set_title(r''+specie,
                     fontsize=8,
                     pad=4)

        ax.set_xlim(0,self.t[-1])

        ax.set_ylim(np.amin(y),
                    np.amax(y))

        tick_length = 4.0
        tick_width = 1.5
        ax.tick_params(labelsize=5,
                       length=tick_length,
                       width=tick_width,
                       direction='in',
                       left=True,
                       right=True,
                       bottom=True,
                       top=True,
                       which='major')

        ax.tick_params(labelsize=5,
                       length=tick_length/1.5,
                       width=tick_width/1.5,
                       direction='in',
                       left=True,
                       right=False,
                       bottom=True,
                       top=False,
                       which='minor')


        # ax.spines['right'].set_visible(False)
        # ax.spines['top'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)


        if rep_coloring == True:

            cmap, c_space = self.choose_colormap(self.N_reps)

            temp_handles = []

            for i_rep in range(self.N_reps):

                p, = ax.plot(self.t,y[:,i_rep],
                             alpha=0.5,
                             linewidth=1.0,
                             zorder=-3,
                             color=cmap(c_space[i_rep]),
                             label=r'{:d}'.format(i_rep+1))

                temp_handles.append(p)

            ax.plot(self.t,ym,
                    alpha=1.0,
                    linewidth=2.0,
                    zorder=-2,
                    color='white')

            ax.plot(self.t,ym,
                    alpha=1.0,
                    linewidth=1.4,
                    zorder=-1,
                    color='black')

            ax.legend(handles=temp_handles,
                      fontsize=6)

        else:
        
            c = 'blueviolet'

            for i_rep in range(self.N_reps):

                ax.plot(self.t,y[:,i_rep],
                            alpha=0.3,
                            linewidth=1.0,
                            zorder=-3,
                            color=c)

            ax.plot(self.t,ym,
                    alpha=1.0,
                    linewidth=2.0,
                    zorder=-2,
                    color='white')

            ax.plot(self.t,ym,
                    alpha=1.0,
                    linewidth=1.4,
                    zorder=-1,
                    color=c)

        return fig, ax


    def prepare_quantity_vs_length_plot_fig(self,y,l,label):

        l_min = np.amin(l)
        l_max = np.amax(l)
        dl = l_max - l_min

        al = np.argsort(l)
        l = l[al]
        y = y[al,:]

        # cmap = colormaps.get_cmap('plasma')
        # c_space_lower_lim = 0.0
        # c_space_upper_lim = 0.9
        # c_space = np.linspace(c_space_lower_lim,
        #                       c_space_upper_lim,
        #                       len(species))

        fig_size = [150,150]

        fig = plt.figure(figsize=(1.15*fig_size[0]*mm,fig_size[1]*mm))

        ax = plt.gca()

        length_norm = matplotlib.colors.Normalize(vmin=l_max,
                                                  vmax=l_min)
        im = colormaps.ScalarMappable(norm=length_norm,
                                      cmap='plasma')
        
        cmap = im.get_cmap()

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)

        cbar = fig.colorbar(im, cax=cax, ticks=[l_min,l_max])
        cbar.set_label(label=r'Length', fontsize=7, labelpad=0)
        cbar.ax.tick_params(labelsize=6)

        ax.set_xlabel(r'$t$ - cell cycle time [s]',
                      fontsize=7,
                      labelpad=1.5)

        ax.set_title(r''+label,
                     fontsize=8,
                     pad=4)

        ax.set_xlim(0,self.t[-1])

        ax.set_ylim(np.amin(y),
                    np.amax(y))

        tick_length = 4.0
        tick_width = 1.5
        ax.tick_params(labelsize=5,
                       length=tick_length,
                       width=tick_width,
                       direction='in',
                       left=True,
                       right=True,
                       bottom=True,
                       top=True,
                       which='major')

        ax.tick_params(labelsize=5,
                       length=tick_length/1.5,
                       width=tick_width/1.5,
                       direction='in',
                       left=True,
                       right=False,
                       bottom=True,
                       top=False,
                       which='minor')

        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)

        for i in range(y.shape[0]):

            c = (l[i] - l_min)/float(dl)

            ax.plot(self.t,y[i,:],
                    alpha=1.0,
                    linewidth=1.0,
                    zorder=-2,
                    color=cmap(c))

        return fig, ax
    

    def plot_individual_concentrations(self,fig_dir,fig_label,
                                       extension,
                                       species,
                                       rep_coloring):

        y = self.get_conc_traces(species)*1.0E+3

        ym = np.mean(y,axis=2)

        for i_specie in range(len(species)):

            specie = species[i_specie]

            fig_file = fig_dir + fig_label + '_conc.'\
                + specie + extension

            specie = specie.replace('_','\_')

            fig, ax = self.prepare_individual_plot_fig(y[i_specie,:,:],
                                                       ym[i_specie,:],
                                                       specie,
                                                       rep_coloring)

            ax.set_ylabel(r'concentration [mM]',
                          fontsize=7,
                          labelpad=1.5)

            fig.savefig(fig_file,
                        dpi=300)

            plt.close()
        
        return

    def plot_individual_counts(self,fig_dir,fig_label,
                               extension,
                               species,
                               rep_coloring):

        y = self.get_species_traces(species)

        ym = np.mean(y,axis=2)

        for i_specie in range(len(species)):

            specie = species[i_specie]

            fig_file = fig_dir + fig_label + '_count.'\
                + specie + extension

            specie = specie.replace('_','\_')

            fig, ax = self.prepare_individual_plot_fig(y[i_specie,:,:],
                                                       ym[i_specie,:],
                                                       specie,
                                                       rep_coloring)

            ax.set_ylabel(r'count [\#]',
                          fontsize=7,
                          labelpad=1.5)

            fig.savefig(fig_file,
                        dpi=300)

            plt.close()
        
        return

    def plot_individual_rxns(self,fig_dir,fig_label,
                             extension,
                             rxns,
                             rep_coloring):

        y = self.get_rxn_traces(rxns)

        ym = np.mean(y,axis=2)

        for i_rxn in range(len(rxns)):

            rxn = rxns[i_rxn]

            fig_file = fig_dir + fig_label + '_flux.'\
                + rxn + extension

            rxn = rxn.replace('_','\_')

            fig, ax = self.prepare_individual_plot_fig(y[i_rxn,:,:],
                                                       ym[i_rxn,:],
                                                       rxn,
                                                       rep_coloring)

            ax.set_ylabel(r'flux [mM/s]',
                          fontsize=7,
                          labelpad=1.5)

            fig.savefig(fig_file,
                        dpi=300)

            plt.close()
        
        return

    def plot_summed_counts(self,fig_dir,fig_label,
                           extension,
                           species,
                           sum_label,
                           rep_coloring):

        
        y = self.get_species_traces(species)
        y = np.sum(y,axis=0)

        ym = np.mean(y,axis=1)

        fig_file = fig_dir + fig_label + '_summedcount.'\
            + sum_label + extension
        sum_file = fig_dir + fig_label + 'sum_details.'\
            + sum_label + '.txt'

        sum_label = sum_label.replace('_','\_')

        fig, ax = self.prepare_individual_plot_fig(y[i_specie,:,:],
                                                   ym[i_specie,:],
                                                   sum_label,
                                                   rep_coloring)

        ax.set_ylabel(r'count [\#]',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        with open(sum_file, 'w') as f:

            for specie in species:

                f.write(specie + '\n')
        
        return

    def plot_multiple_conc_per_rep(self,fig_dir,fig_label,
                                   extension,
                                   species,
                                   multiple_label):

        cmap, c_space = self.choose_colormap(len(species))

        y = self.get_conc_traces(species)*1.0E+3

        ym = np.mean(y,axis=2)

        for i_rep in range(self.N_reps):

            fig_size = [120,120*l_w_ratio]

            fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

            ax = plt.gca()

            ax.set_xlabel(r'$t$ - cell cycle time [s]',
                          fontsize=7,
                          labelpad=1.5)

            ax.set_title(r''+multiple_label.replace('_','\_'),
                         fontsize=8,
                         pad=4)

            ax.set_xlim(0,self.t[-1])

            ax.set_ylim(np.amin(y),
                        np.amax(y))

            tick_length = 4.0
            tick_width = 1.5
            ax.tick_params(labelsize=5,
                           length=tick_length,
                           width=tick_width,
                           direction='in',
                           left=True,
                           right=True,
                           bottom=True,
                           top=True,
                           which='major')

            ax.tick_params(labelsize=5,
                           length=tick_length/1.5,
                           width=tick_width/1.5,
                           direction='in',
                           left=True,
                           right=False,
                           bottom=True,
                           top=False,
                           which='minor')

            ax.spines['left'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['top'].set_linewidth(1.5)

            fig_file = fig_dir + fig_label + '_' + multiple_label + '_multiconc.'\
                + '{:05}'.format(i_rep+1) + extension


            temp_handles = []
            
            for i_specie in range(len(species)):

                specie = species[i_specie]
                specie = specie.replace('_','\_')

                temp_color = cmap(c_space[i_specie])

                p, = ax.plot(self.t,y[i_specie,:,i_rep],
                             color=temp_color,
                             zorder=1,
                             linewidth=1.5,
                             linestyle='-',
                             alpha=1.0,
                             label=r''+specie)

                ax.plot(self.t,y[i_specie,:,i_rep],
                        color='white',
                        zorder=0,
                        linewidth=2.0,
                        linestyle='-',
                        alpha=1.0,
                        label=r''+specie)

                temp_handles.append(p)

                ax.plot(self.t,ym[i_specie,:],
                        markeredgecolor=temp_color,
                        markerfacecolor='white',
                        markeredgewidth=0.3,
                        zorder=-1,
                        marker='.',
                        markersize=6.0,
                        alpha=0.8,
                        linestyle='None',
                        markevery=25)

            ax.set_ylabel(r'concentration [mM]',
                          fontsize=7,
                          labelpad=1.5)

            ax.legend(handles=temp_handles,
                      fontsize=8)

            fig.savefig(fig_file,
                        dpi=600)

            plt.close()

        return

    def plot_multiple_rxns_per_rep(self,fig_dir,fig_label,
                                   extension,
                                   rxns,
                                   multiple_label):

        cmap, c_space = self.choose_colormap(len(rxns))

        y = self.get_rxn_traces(rxns)

        ym = np.mean(y,axis=2)

        for i_rep in range(self.N_reps):

            fig_size = [120,120*l_w_ratio]

            fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

            ax = plt.gca()

            ax.set_xlabel(r'$t$ - cell cycle time [s]',
                          fontsize=7,
                          labelpad=1.5)

            ax.set_title(r''+multiple_label.replace('_','\_'),
                         fontsize=8,
                         pad=4)

            ax.set_xlim(0,self.t[-1])

            ax.set_ylim(np.amin(y),
                        np.amax(y))

            tick_length = 4.0
            tick_width = 1.5
            ax.tick_params(labelsize=5,
                           length=tick_length,
                           width=tick_width,
                           direction='in',
                           left=True,
                           right=True,
                           bottom=True,
                           top=True,
                           which='major')

            ax.tick_params(labelsize=5,
                           length=tick_length/1.5,
                           width=tick_width/1.5,
                           direction='in',
                           left=True,
                           right=False,
                           bottom=True,
                           top=False,
                           which='minor')

            ax.spines['left'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['top'].set_linewidth(1.5)

            fig_file = fig_dir + fig_label + '_' + multiple_label + '_multirxn.'\
                + '{:05}'.format(i_rep+1) + extension


            temp_handles = []
            
            for i_rxn in range(len(rxns)):

                rxn = rxns[i_rxn]
                rxn = rxn.replace('_','\_')

                temp_color = cmap(c_space[i_rxn])

                p, = ax.plot(self.t,y[i_rxn,:,i_rep],
                             color=temp_color,
                             zorder=1,
                             linewidth=1.5,
                             linestyle='-',
                             alpha=1.0,
                             label=r''+rxn)

                ax.plot(self.t,y[i_rxn,:,i_rep],
                        color='white',
                        zorder=0,
                        linewidth=2.0,
                        linestyle='-',
                        alpha=1.0,
                        label=r''+rxn)

                temp_handles.append(p)

                ax.plot(self.t,ym[i_rxn,:],
                        markeredgecolor=temp_color,
                        markerfacecolor='white',
                        markeredgewidth=0.3,
                        zorder=-1,
                        marker='.',
                        markersize=6.0,
                        alpha=0.8,
                        linestyle='None',
                        markevery=25)

            ax.set_ylabel(r'flux [mM/s]',
                          fontsize=7,
                          labelpad=1.5)

            ax.legend(handles=temp_handles,
                      fontsize=8)

            fig.savefig(fig_file,
                        dpi=600)

            plt.close()

        return

    def plot_groups(self,fig_dir,fig_label,
                    extension,
                    grouped_species,
                    group_labels,
                    rep_coloring):

        y = self.get_group_traces(group_labels,
                                  grouped_species)
        ym = np.mean(y,axis=2)

        N_groups = len(group_labels)
        
        for i_group in range(N_groups):

            temp_label = group_labels[i_group]

            fig_file = fig_dir + fig_label + '_groupcount.'\
                + temp_label + extension

            temp_species = grouped_species[temp_label]['species']
            temp_weights = grouped_species[temp_label]['weights']

            temp_ref_specie = grouped_species[temp_label]['ref_specie']

            group_equation = '$'
            if temp_ref_specie != None:
                 group_equation += '\\normalfont{' + temp_ref_specie.replace('_','\_') + \
                    '}(t_0)'
            
            for i_specie in range(len(temp_species)):

                w = temp_weights[i_specie]
                
                if w < 0.0:
                    group_equation += '{:.1f}'.format(w)
                else:
                    if (i_specie == 0) and (temp_ref_specie == None):
                        group_equation += '{:.1f}'.format(w)
                    else:
                        group_equation += '+{:.1f}'.format(w)

                group_equation += '\\times\\normalfont{' + \
                    temp_species[i_specie].replace('_','\_') + \
                    '}'

            group_equation = group_equation + '$'

            # print(group_equation)

            fig, ax = self.prepare_individual_plot_fig(y[i_group,:,:],
                                                       ym[i_group,:],
                                                       group_equation,
                                                       rep_coloring)

            ax.set_ylabel(r'count [\#]',
                          fontsize=7,
                          labelpad=1.5)

            fig.savefig(fig_file,
                        dpi=300)

            plt.close()

        return

    def plot_rxn_groups(self,fig_dir,fig_label,
                        extension,
                        grouped_rxns,
                        group_labels,
                        rep_coloring):

        y = self.get_group_rxn_traces(group_labels,
                                      grouped_rxns)
        ym = np.mean(y,axis=2)

        N_groups = len(group_labels)
        
        for i_group in range(N_groups):

            temp_label = group_labels[i_group]

            fig_file = fig_dir + fig_label + '_groupflux.'\
                + temp_label + extension

            temp_rxns = grouped_rxns[temp_label]['rxns']
            temp_weights = grouped_rxns[temp_label]['weights']

            group_equation = '$'
            
            for i_rxn in range(len(temp_rxns)):

                w = temp_weights[i_rxn]
                
                if w < 0.0:
                    group_equation += '{:.1f}'.format(w)
                else:
                    if (i_rxn == 0):
                        group_equation += '{:.1f}'.format(w)
                    else:
                        group_equation += '+{:.1f}'.format(w)

                group_equation += '\\times\\normalfont{' + \
                    temp_rxns[i_rxn].replace('_','\_') + \
                    '}'

            group_equation = group_equation + '$'

            # print(group_equation)

            fig, ax = self.prepare_individual_plot_fig(y[i_group,:,:],
                                                       ym[i_group,:],
                                                       group_equation,
                                                       rep_coloring)

            ax.set_ylabel(r'flux [mM/s]',
                          fontsize=7,
                          labelpad=1.5)

            fig.savefig(fig_file,
                        dpi=300)

            plt.close()

        return

    def plot_group_relative_length(self,fig_dir,fig_label,
                                   extension,
                                   grouped_species,
                                   group_labels,
                                   lengths):

        cmap = colormaps.get_cmap('plasma')
        c_space_lower_lim = 0.0
        c_space_upper_lim = 0.9
        c_space = np.linspace(c_space_lower_lim,
                              c_space_upper_lim,
                              self.Nt)

        y = self.get_group_traces(group_labels,
                                  grouped_species)

        ym = np.mean(y,axis=2)

        fig_size = [120,120*l_w_ratio]

        fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

        ax = plt.gca()

        ax.set_xlabel(r'$L$ - Gene Length [bp]',
                      fontsize=7,
                      labelpad=1.5)

        ax.set_title(r'Group vs. Length',
                     fontsize=8,
                     pad=4)

        ax.set_xlim(np.amin(lengths),np.amax(lengths))

        ax.set_ylim(np.amin(ym),
                    np.amax(ym))

        tick_length = 4.0
        tick_width = 1.5
        ax.tick_params(labelsize=5,
                       length=tick_length,
                       width=tick_width,
                       direction='in',
                       left=True,
                       right=True,
                       bottom=True,
                       top=True,
                       which='major')

        ax.tick_params(labelsize=5,
                       length=tick_length/1.5,
                       width=tick_width/1.5,
                       direction='in',
                       left=True,
                       right=False,
                       bottom=True,
                       top=False,
                       which='minor')

        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)


        fig_file = fig_dir + fig_label + '_groupsvlengths'\
            + extension

        for i_t in range(self.Nt):

            temp_color = cmap(c_space[i_t])

            ax.plot(lengths,ym[:,i_t],
                    markeredgecolor=temp_color,
                    markerfacecolor=temp_color,
                    markeredgewidth=0.3,
                    zorder=-1,
                    marker='.',
                    markersize=2.0,
                    alpha=1.0,
                    linestyle='None')
                           
                           
        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        return

    def covariance_matrix(self,y,z):
        
        cov = np.zeros((y.shape[0],z.shape[0]),dtype=np.float32)

        for i_rep in range(y.shape[1]):

            cov += np.outer(y[:,i_rep],z[:,i_rep])

        cov = cov/(y.shape[1] - 1)

        return cov

    def correlation_matrix(self,y,z):

        corr = self.covariance_matrix(y,z)

        y_invstd = np.reciprocal(np.std(y,axis=1,ddof=1))
        z_invstd = np.reciprocal(np.std(z,axis=1,ddof=1))

        corr = corr*np.outer(y_invstd,z_invstd)

        return corr

    def time_correlations_vs_final(self,y):

        # Nt by N_specie by N_specie matrix
        Nt_corr = y.shape[1]
        corr_t = np.zeros((y.shape[0],y.shape[0],Nt_corr))

        for i_t in range(y.shape[1]):

            corr_t[:,:,i_t] = self.correlation_matrix(y[:,Nt_corr-1,:],
                                                      y[:,Nt_corr-1-i_t,:])

        return corr_t

    def plot_time_correlations_vs_final(self,fig_dir,fig_label,
                                        extension,
                                        corr_label,
                                        conc_species,
                                        group_labels_count,
                                        grouped_species,
                                        group_labels_rxn,
                                        grouped_rxns):

        y, merged_labels = self.get_conc_group_count_group_rxn_traces(conc_species,
                                                                      group_labels_count,
                                                                      grouped_species,
                                                                      group_labels_rxn,
                                                                      grouped_rxns)

        N_conc = len(conc_species)
        N_count = len(group_labels_count)
        N_rxn_flux = len(group_labels_rxn)

        conc_indices = range(0,N_conc)
        count_indices = range(N_conc,N_conc+N_count)
        rxn_indices = range(N_conc+N_count,N_conc+N_count+N_rxn_flux)

        conc_color = 'mistyrose'
        count_color = 'paleturquoise'
        rxn_color = 'lemonchiffon'

        y = y - np.mean(y,axis=2)[:,:,np.newaxis]

        corr_t = self.time_correlations_vs_final(y)

        N_var = len(merged_labels)

        per_plot_size = 25.4
        fig_size = [N_var*per_plot_size,N_var*per_plot_size]

        fig_file = fig_dir + fig_label + '_corrt.'\
                + corr_label + extension

        fig, axs = plt.subplots(figsize=(fig_size[0]*mm,fig_size[1]*mm),
                                ncols=N_var,
                                nrows=N_var)

        for i in range(N_var):
            for j in range(N_var):

                axs[i,j].set_xlim(0,self.t[-1])

                axs[i,j].set_ylim(-1.0,1.0)

                tick_length = 2.0
                tick_width = 0.75
                axs[i,j].tick_params(labelsize=5,
                               length=tick_length,
                               width=tick_width,
                               direction='out',
                               left=True,
                               right=False,
                               bottom=False,
                               top=False,
                               which='major')

                axs[i,j].tick_params(labelsize=5,
                               length=tick_length/1.5,
                               width=tick_width/1.5,
                               direction='in',
                               left=False,
                               right=False,
                               bottom=False,
                               top=False,
                               which='minor')

                axs[i,j].set_xticks([])
                axs[i,j].set_xticklabels(labels=[])
                axs[i,j].set_yticks([-1,1])
                axs[i,j].set_yticklabels(labels=[r'-1',r'1'])


                axs[i,j].spines['left'].set_linewidth(1.0)
                axs[i,j].spines['right'].set_linewidth(1.0)
                axs[i,j].spines['bottom'].set_linewidth(1.0)
                axs[i,j].spines['top'].set_linewidth(1.0)



                axs[i,j].spines['right'].set_visible(False)
                

                if i == 0:
                    axs[i,j].spines['bottom'].set_position('center')
                    axs[i,j].spines['top'].set_visible(False)
                    
                else:
                    axs[i,j].spines['top'].set_position('center')
                    axs[i,j].spines['bottom'].set_visible(False)
                    


                if (i == 0) or (i == N_var - 1):
                    
                    temp_label = r'' + merged_labels[j].replace('_','\_')
                    
                    if i == 0:
                        axs[i,j].xaxis.set_label_position('top')
                        axs[i,j].set_xlabel(temp_label,
                                            horizontalalignment='left',
                                            rotation=45)
                    if i == N_var - 1:
                        axs[i,j].xaxis.set_label_position('bottom')
                        axs[i,j].set_xlabel(temp_label,
                                            horizontalalignment='right',
                                            rotation=45)
                    
                if (j == 0) or (j == N_var - 1):
                    temp_label = r'' + merged_labels[i].replace('_','\_')

                    if j == 0:
                        axs[i,j].yaxis.set_label_position('left')
                        axs[i,j].set_ylabel(temp_label,
                                            horizontalalignment='right',
                                            rotation=45)
                    if j == N_var - 1:
                        axs[i,j].yaxis.set_label_position('right')
                        axs[i,j].set_ylabel(temp_label,
                                            horizontalalignment='left',
                                            rotation=45)
                        

                if (i in conc_indices) and (j in conc_indices):
                    axs[i,j].set_facecolor(conc_color)
                if (i in count_indices) and (j in count_indices):
                    axs[i,j].set_facecolor(count_color)
                if (i in rxn_indices) and (j in rxn_indices):
                    axs[i,j].set_facecolor(rxn_color)

                if i == j:
                    temp_color = 'limegreen'
                else:
                    temp_color = 'magenta'

                axs[i,j].plot(self.t,corr_t[i,j,:],
                              linewidth=1.0,
                              color=temp_color,
                              ls='-',
                              zorder=-3)

        fig.suptitle(r'Lagged time-correlation of $Y$ with respect to $X$ at final time ($T_f$) --- i.e., $\rho_{XY}(T_f,\tau)=\frac{\mathrm{Cov}[X(T_f),Y(T_f-\tau)]}{\sqrt{\mathrm{Var}[X(T_f)]\times\mathrm{Var}[Y(T_f-\tau)]}}, \forall\tau\in[0,T_f]$',
                     fontsize=32)
        fig.supylabel(r'$X$ in multivariate random process',
                      fontsize=28)
        fig.supxlabel(r'$Y$ in multivariate random process',
                      fontsize=28)

        fig.savefig(fig_file,
                    dpi=300)
        
        plt.close()
        
        return

    # def windowed_covariance_matrix(self,target_time_window,species,trace_type):
        
    #     if (target_time_window[0] > target_time_window[1]):

    #         raise ValueError('Impossible time window')

    #     i_t_min = 0
    #     i_t_max = self.Nt - 1

    #     while ((self.t[i_t_min] < target_time_window[0]) and
    #            (i_t_min < i_t_max)):

    #         i_t_min += 1

    #     while ((self.t[i_t_max] > target_time_window[1]) and
    #            (i_t_max > i_t_min)):

    #         i_t_max -= 1

    #     actual_time_window = np.array([self.t[i_t_min],self.t[i_t_max]],
    #                                   dtype=np.float32)

    #     if trace_type == 'count':
    #         y = self.get_species_traces(species)
    #     elif trace_type == 'conc':
    #         y = self.get_conc_traces(species)
    #     else:
    #         raise ValueError('Choose [count] or [conc]')
        
    #     y = y[:,i_t_min:i_t_max+1,:].astype(np.float32)
    #     y = y - np.mean(y,axis=2)[:,:,np.newaxis]

    #     cov = np.zeros((y.shape[1],y.shape[0],y.shape[0]),
    #                    dtype=np.float32)

    #     for i_t in range(y.shape[1]):

    #         cov[i_t,:,:] += self.covariance_matrix(y[:,i_t,:])

    #     cov = np.mean(cov,axis=0)

    #     return actual_time_window, cov

    # def windowed_correlation_matrix(self,target_time_window,species,trace_type):
        
    #     if (target_time_window[0] > target_time_window[1]):

    #         raise ValueError('Impossible time window')

    #     i_t_min = 0
    #     i_t_max = self.Nt - 1

    #     while ((self.t[i_t_min] < target_time_window[0]) and
    #            (i_t_min < i_t_max)):

    #         i_t_min += 1

    #     while ((self.t[i_t_max] > target_time_window[1]) and
    #            (i_t_max > i_t_min)):

    #         i_t_max -= 1

    #     actual_time_window = np.array([self.t[i_t_min],self.t[i_t_max]],
    #                                   dtype=np.float32)

    #     if trace_type == 'count':
    #         y = self.get_species_traces(species)
    #     elif trace_type == 'conc':
    #         y = self.get_conc_traces(species)
    #     else:
    #         raise ValueError('Choose [count] or [conc]')
        
    #     y = y[:,i_t_min:i_t_max+1,:].astype(np.float32)
    #     y = y - np.mean(y,axis=2)[:,:,np.newaxis]

    #     corr = np.zeros((y.shape[1],y.shape[0],y.shape[0]),
    #                     dtype=np.float32)

    #     for i_t in range(y.shape[1]):

    #         corr[i_t,:,:] += self.correlation_matrix(y[:,i_t,:])

    #     corr = np.mean(corr,axis=0)

    #     return actual_time_window, corr

    # def plot_correlations(self,fig_dir,fig_label,
    #                       extension,
    #                       target_time_window,
    #                       species,
    #                       trace_type,
    #                       corr_label):

    #     fig_size = [150,150]

    #     fig_file = fig_dir + fig_label + '_corr.'\
    #             + corr_label + extension

    #     fig = plt.figure(figsize=(1.15*fig_size[0]*mm,fig_size[1]*mm))

    #     ax = plt.gca()

    #     ax.spines['right'].set_visible(False)
    #     ax.spines['top'].set_visible(False)
    #     ax.spines['left'].set_visible(False)
    #     ax.spines['bottom'].set_visible(False)

    #     actual_time_window, corr = self.windowed_correlation_matrix(target_time_window,
    #                                                                 species,
    #                                                                 trace_type)

    #     mat_norm = matplotlib.colors.TwoSlopeNorm(vmin=-1.0,
    #                                               vcenter=0.0,
    #                                               vmax=1.0)

    #     im = ax.imshow(corr,cmap='PiYG',norm=mat_norm)

    #     divider = make_axes_locatable(ax)
    #     cax = divider.append_axes("right", size="5%", pad=0.1)

    #     cbar = fig.colorbar(im, cax=cax, ticks=[-1.0,0,1.0])
    #     cbar.set_label(label=r'Correlation', fontsize=7, labelpad=0)
    #     cbar.ax.tick_params(labelsize=6)

    #     s_expand = 0.05*float(corr.shape[0])
    #     d_offset = 0.5

    #     ax.set_xlim(xmin=0.0-d_offset-s_expand,
    #                 xmax=float(corr.shape[0])+d_offset-1.0)
    #     ax.set_ylim(ymin=float(corr.shape[0])+d_offset-1.0+s_expand,
    #                 ymax=0.0-d_offset)

    #     ticks = np.arange(0,corr.shape[0],dtype=np.int32)

    #     tick_labels = []
    #     for specie in species:
    #         temp_specie = specie
    #         #temp_specie = specie[2:-2]
    #         tick_labels.append(r''+temp_specie.replace('_','\_'))
            
    #     ax.set_xticks(ticks)
    #     ax.set_xticklabels(tick_labels,
    #                        fontsize=3,
    #                        ha='center',
    #                        rotation=90)

    #     ax.set_yticks(ticks)
    #     ax.set_yticklabels(tick_labels,
    #                        fontsize=3,
    #                        va='center',
    #                        ha='right',
    #                        rotation=0)

    #     ax.set_title(r'Time-Averaged~({:d}s - {:d}s) Correlations'.format(actual_time_window[0].astype(np.int),
    #                                                                     actual_time_window[1].astype(np.int)),
    #                  fontsize=12)

    #     fig.savefig(fig_file,
    #                 dpi=300)
        
    #     plt.close()

    #     return

    def plot_volume(self,fig_dir,fig_label,
                    extension,
                    rep_coloring):

        # absolute volume
        
        y = self.volumes*1.0E+15

        ym = np.mean(y,axis=1)

        fig_file = fig_dir + fig_label + '_absvol' + extension


        fig, ax = self.prepare_individual_plot_fig(y,
                                                   ym,
                                                   'Absolute Volume Growth',
                                                   rep_coloring)

        ax.set_ylabel(r'Volume [fL]',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        # relative volume

        y = self.volumes
        y = y*np.reciprocal(y[0,:])

        ym = np.mean(y,axis=1)

        # y = np.log2(y)
        # ym = np.log2(ym)

        fig_file = fig_dir + fig_label + '_relvol' + extension


        fig, ax = self.prepare_individual_plot_fig(y,
                                                   ym,
                                                   'Relative Volume Growth',
                                                   rep_coloring)

        # ax.set_ylabel(r'$\log_2[V(t)/V(t_0)]$',
        #               fontsize=7,
        #               labelpad=1.5)

        ax.set_ylabel(r'$V(t)/V(t_0)$',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()
        
        return

    def plot_surfacearea(self,fig_dir,fig_label,
                         extension,
                         rep_coloring):

        species = ['SA_total','SA_ptn','SA_lipid']

        z = self.get_species_traces(species)

        # absolute SA

        y = z[0,:,:].astype(np.float32)
        
        ym = np.mean(y,axis=1)

        fig_file = fig_dir + fig_label + '_absSA' + extension


        fig, ax = self.prepare_individual_plot_fig(y,
                                                   ym,
                                                   'Absolute Surface Area Growth',
                                                   rep_coloring)

        ax.set_ylabel(r'Surface Area [nm$^2$]',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        # relative SA

        y = y*np.reciprocal(y[0,:])

        ym = np.mean(y,axis=1)

        # y = np.log2(y)
        # ym = np.log2(ym)

        fig_file = fig_dir + fig_label + '_relSA' + extension

        fig, ax = self.prepare_individual_plot_fig(y,
                                                   ym,
                                                   'Relative Surface Area Growth',
                                                   rep_coloring)

        # ax.set_ylabel(r'$\log_2[SA(t)/SA(t_0)]$',
        #               fontsize=7,
        #               labelpad=1.5)

        ax.set_ylabel(r'$SA(t)/SA(t_0)$',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        # membrane protein to lipid ratio

        y = z[1:,:,:].astype(np.float32)
        y = y[0,:,:]/y[1,:,:]

        ym = np.mean(y,axis=1)

        fig_file = fig_dir + fig_label + '_ratioSA' + extension


        fig, ax = self.prepare_individual_plot_fig(y,
                                                   ym,
                                                   'Protein:Lipid Surface Area Ratio',
                                                   rep_coloring)
        ax.set_ylabel(r'Protein:Lipid',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        # relative volume

        y = z[0,:,:].astype(np.float32)

        y = (1.0/(6*np.sqrt(np.pi)))*np.power(y,3.0/2.0)
        w = self.volumes*1.0E-3 # cubic meters
        w = w*np.power(1.0E+9,3.0) # cubic nanometers
        y = y/w

        ym = np.mean(y,axis=1)

        fig_file = fig_dir + fig_label + '_redvol' + extension


        fig, ax = self.prepare_individual_plot_fig(y,
                                                   ym,
                                                   'Reduced Volume',
                                                   rep_coloring)
        ax.set_ylabel(r'$V(t)/V_{\star}[SA(t)]$',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()
        
        return

    def plot_DNAcontent(self,fig_dir,fig_label,
                        extension,
                        rep_coloring):

        # absolute chromosome
        
        y = self.get_specie_trace('chromosome').astype(np.float32)
        y = 10.0*y

        ym = np.mean(y,axis=1)

        fig_file = fig_dir + fig_label + '_absDNA' + extension


        fig, ax = self.prepare_individual_plot_fig(y,
                                                   ym,
                                                   'Absolute DNA Content Growth',
                                                   rep_coloring)

        ax.set_ylabel(r'$N_{DNA}$ - Chromosomal DNA [bp]',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        # relative chromosome

        y = y*np.reciprocal(y[0,:])

        ym = np.mean(y,axis=1)

        # y = np.log2(y)
        # ym = np.log2(ym)

        fig_file = fig_dir + fig_label + '_relDNA' + extension


        fig, ax = self.prepare_individual_plot_fig(y,
                                                   ym,
                                                   'Relative DNA Content Growth',
                                                   rep_coloring)

        # ax.set_ylabel(r'$\log_2[V(t)/V(t_0)]$',
        #               fontsize=7,
        #               labelpad=1.5)

        ax.set_ylabel(r'$N_{DNA}(t)/N_{DNA}(t_0)$',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()
        
        return

    def plot_protein_metrics(self,fig_dir,fig_label,
                             extension,
                             lt,l,
                             rep_coloring):

        l = l/3.0


        # plot the absolute number of amino acids incorporated into proteins
        
        fig_file = fig_dir + fig_label + '_absProteinAAs' + extension

        group_labels, grouped_species, y = self.get_protein_traces(lt)
        z = np.tensordot(l,y,axes=1)

        print(z.shape)

        zm = np.mean(z,axis=1)

        fig, ax = self.prepare_individual_plot_fig(z,
                                                   zm,
                                                   'Absolute Protein AA Content Growth',
                                                   rep_coloring)

        ax.set_ylabel(r'$N_{protein-AA}$ - \# Amino Acids',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        # plot the relative number of amino acids incorporated into proteins

        fig_file = fig_dir + fig_label + '_relProteinAAs' + extension

        z = z*np.reciprocal(z[0,:].astype(np.float32))

        zm = np.mean(z,axis=1)

        fig, ax = self.prepare_individual_plot_fig(z,
                                                   zm,
                                                   'Relative Protein AA Content Growth',
                                                   rep_coloring)

        ax.set_ylabel(r'$N_{protein-AA}(t)/N_{protein-AA}(t_0)$',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        # plot the protein creation relative to the protein length

        fig_file = fig_dir + fig_label + '_relProteinVsLength' + extension

        y0_inv = np.reciprocal(y[:,0,:].astype(np.float32))
        z = y*np.expand_dims(y0_inv,axis=1)

        zm = np.mean(z,axis=2)

        fig, ax = self.prepare_quantity_vs_length_plot_fig(zm,
                                                           l,
                                                           'Protein Production vs. Length')

        ax.set_ylabel(r'$N(t)/N(t_0)$',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()
        

        return

    def plot_protein_transcript_metrics(self,fig_dir,fig_label,
                                        extension,
                                        lt,l,
                                        rep_coloring):


        # plot the absolute number of nucleotides incorporated into transcripts for proteins
        
        fig_file = fig_dir + fig_label + '_absProteinTranscriptNTs' + extension

        group_labels, grouped_species, y = self.get_transcript_traces(lt)
        z = np.tensordot(l,y,axes=1)

        print(z.shape)

        zm = np.mean(z,axis=1)

        fig, ax = self.prepare_individual_plot_fig(z,
                                                   zm,
                                                   'Transcript NT Content Growth',
                                                   rep_coloring)

        ax.set_ylabel(r'$N_{transcript-NT}$ - \# Nucleotides',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()

        # plot the relative number of amino acids incorporated into proteins

        # fig_file = fig_dir + fig_label + '_relProteinTranscriptNTs' + extension

        # z = z*np.reciprocal(z[0,:].astype(np.float32))

        # zm = np.mean(z,axis=1)

        # fig, ax = self.prepare_individual_plot_fig(z,
        #                                            zm,
        #                                            'Relative Transcript NT Content Growth')

        # ax.set_ylabel(r'$N_{transcript-NT}(t)/N_{transcript-NP}(t_0)$',
        #               fontsize=7,
        #               labelpad=1.5)

        # fig.savefig(fig_file,
        #             dpi=300)

        # plt.close()

        # plot the protein creation relative to the protein length

        fig_file = fig_dir + fig_label + '_relTranscriptsVsLength' + extension

        zm = np.mean(y,axis=2)

        fig, ax = self.prepare_quantity_vs_length_plot_fig(zm,
                                                           l,
                                                           'Transcript Count vs. Length')

        ax.set_ylabel(r'$N(t)/N(t_0)$',
                      fontsize=7,
                      labelpad=1.5)

        fig.savefig(fig_file,
                    dpi=300)

        plt.close()
        

        return
 


    

