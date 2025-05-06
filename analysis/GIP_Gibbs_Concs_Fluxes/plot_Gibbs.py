"""
Functions to plot the Gibbs Free Energy Comparison
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.ticker as tck
import matplotlib.patches as patches
import numpy as np
import pandas as pd

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage[helvet]{sfmath}'

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

small_ticklabelsize=6
large_ticklabelsize=7
titlesize=20
fig_titlesize=26
axes_labelsize=8
axes_labelsize_small=8
subplot_labelsize=8
small_legendsize=6
colorbar_labelsize=8

mm = 1/25.4

def return_df(rxns_list, dG0_equilibrator, dG0_kinetics, dG_equilibrator, dG_kinetics, RTlogQ):
    """ 
    
    Save the calculated standard/ Gibbs Free Energy change to Excel  
    """
    unit = 'kJ/mol'

    df = pd.DataFrame({
        'reactions': rxns_list,
        f'dG0_thermo ({unit})': dG0_equilibrator,
        f'dG0_kinetic ({unit})': dG0_kinetics,
        f'dG_thermo ({unit})': dG_equilibrator,
        f'dG_kinetics ({unit})': dG_kinetics,
        f'RTlogQ ({unit})': RTlogQ
    })
    
    # df.to_excel(file_path, index=False)

    return df

def plot_deltaG0(dG0_equilibrator, dG0_kinetics, rxns_list, rxn_rev, balance_test, fig_path):

    N_t = 1
    
    l_w_ratio = 0.39

    fig_size = [174,174*l_w_ratio]

    heights = [1.0,1.0,1.0]

    kinetics_color = 'darkcyan'
    equilibrator_color = 'darkorange'

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    # set y-axes labels
    ax.set_ylabel(r'$\Delta G$~[kJ/mol]',fontsize=14)

    # set y-axes scales and gridlines
    ax.set_axisbelow(True)
    ax.tick_params(axis='y',labelsize=10)
    # ax.set_yscale('log')
    ax.yaxis.grid(which='major')

    # set x-axes tick labels
    rxn_labels = []; N_rxn = len(rxns_list)

    for i in range(N_rxn):
        temp_label = rxns_list[i]
        temp_label = r''+temp_label.replace('_','\_')
        rxn_labels.append(temp_label)
        
    # x coordinates
    x_coords = np.arange(0,N_rxn)

    # determine if reactions are irreversible
    rev_rxns = []
    irrev_rxns = []
    for i in range(N_rxn):
        if rxn_rev[i] == 1:
            rev_rxns.append(i)
        else:
            irrev_rxns.append(i)

    # x-coordinates for rev vs. irrev reactions in kinetics
    rev_rxns = rxn_rev == 1
    rev_x_coords = x_coords[rev_rxns]

    # x-coordinates for bal vs. unbal reactions in FBA
    bal_rxns = balance_test == 0
    bal_x_coords = x_coords[bal_rxns]

    legend_elements = []


    for i_t in range(N_t):
        # Plot thermodynamic dG
        temp_y = dG0_equilibrator[bal_rxns]
        if i_t == 0:
            current_min = np.nanmin(temp_y)
            current_max = np.nanmax(temp_y)
        else:
            current_min = min(current_min,np.nanmin(temp_y))
            current_max = max(current_max,np.nanmax(temp_y))
    #     a = ax.scatter(bal_x_coords,temp_y,color=equilibrator_color,marker=time_markers[i_t],label=r'equilibrator - '+time_labels[i_t],alpha=0.7)
        a = ax.scatter(bal_x_coords,temp_y,color=equilibrator_color,marker='+',label=r'equilibrator~($\Delta G_0$)',alpha=0.7,s=100)
        legend_elements.append(a)

    for i_t in range(N_t):
        temp_y = dG0_kinetics[rev_rxns]
        current_min = min(current_min,np.nanmin(temp_y))
        current_max = max(current_max,np.nanmax(temp_y))
    #     a = ax.scatter(rev_x_coords,temp_y,color=kinetics_color,marker=time_markers[i_t],label=r'kinetics - '+time_labels[i_t],alpha=0.7)
        a = ax.scatter(rev_x_coords,temp_y,color=kinetics_color,marker='x',label=r'kinetics~($\Delta G_0$)',alpha=0.7,s=60)
        legend_elements.append(a)

    # adjust the y limits
    intervals = 5
    current_min = intervals*np.floor_divide(current_min,intervals)
    current_max = intervals*(np.floor_divide(current_max,intervals) + 1)
        
    # ax.autoscale(enable=True,axis='y',tight=True)
    ax.set_autoscaley_on(False)
    ax.set_ylim(ymin=current_min,ymax=current_max)
    
    # set yticks
    up = current_max//(2*intervals); down = current_min//(2*intervals)+1
    yticks = 2*intervals*np.arange(down, up+1)
    ax.set_yticks(yticks)

    # draw vertical lines to separate reactions
    vline_x_coords = np.arange(1,N_rxn) - 0.5
    # ax_ylims = ax.get_ylim()
    ax.vlines(vline_x_coords,ymin=current_min,ymax=current_max,linestyles='dotted',linewidth=0.75,alpha=0.25)
    

    # adjust all x-ticks
    ax.set_xlim(xmin=-0.75,xmax=N_rxn-0.25)
    ax.set_xticks(x_coords)
    ax.set_xticklabels(rxn_labels, rotation=90, ha="center", fontsize=8)

    ax.legend(handles=legend_elements,fontsize=8,loc='lower right',title=r'$\Delta G=\Delta G_0$',title_fontsize=9)

    plt.tight_layout()

    fig.savefig(fig_path, dpi=600)

    return None




def plot_deltaG(dG0_equilibrator, dG0_kinetics, RTlogQ, rxns_list, rxn_rev, balance_test, fig_path, errorbar=False):

    N_t = 1
    
    l_w_ratio = 0.39

    fig_size = [174,174*l_w_ratio]

    heights = [1.0,1.0,1.0]

    kinetics_color = 'darkcyan'
    equilibrator_color = 'darkorange'

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    # set y-axes labels
    ax.set_ylabel(r'$\Delta G$~[kJ/mol]',fontsize=14)

    # set y-axes scales and gridlines
    ax.set_axisbelow(True)
    ax.tick_params(axis='y',labelsize=10)
    # ax.set_yscale('log')
    ax.yaxis.grid(which='major')

    # set x-axes tick labels
    rxn_labels = []; N_rxn = len(rxns_list)

    for i in range(N_rxn):
        temp_label = rxns_list[i]
        temp_label = r''+temp_label.replace('_','\_')
        rxn_labels.append(temp_label)
        
    # x coordinates
    x_coords = np.arange(0,N_rxn)


    # x-coordinates for rev vs. irrev reactions in kinetics
    rev_rxns = rxn_rev == 1
    rev_x_coords = x_coords[rev_rxns]

    # x-coordinates for bal vs. unbal reactions in FBA
    bal_rxns = balance_test==0
    bal_x_coords = x_coords[bal_rxns]

    legend_elements = []


    # Plot thermodynamic dG
    temp_y = dG0_equilibrator[bal_rxns] + RTlogQ[bal_rxns,0] # 0 means ensemble averaged conc

    current_min = np.nanmin(temp_y)
    current_max = np.nanmax(temp_y)
    a = ax.scatter(bal_x_coords,temp_y,color=equilibrator_color,marker='+',label=r'equilibrator~($\Delta G_0$)',alpha=0.7,s=100)

    if errorbar:
        yerr_upper = np.nanmax(RTlogQ[bal_rxns,:], axis=1) - RTlogQ[bal_rxns,0]
        yerr_lower = - np.nanmin(RTlogQ[bal_rxns,:], axis=1) + RTlogQ[bal_rxns,0]
        ax.errorbar(bal_x_coords,temp_y, yerr=[yerr_lower, yerr_upper], color=equilibrator_color, linestyle='None', capsize=5 )

    legend_elements.append(a)

    # repeat for kinetic deltaG
    temp_y = dG0_kinetics[rev_rxns] + RTlogQ[rev_rxns, 0]

    current_min = min(current_min,np.nanmin(temp_y))
    current_max = max(current_max,np.nanmax(temp_y))

    a = ax.scatter(rev_x_coords,temp_y,color=kinetics_color,marker='x',label=r'kinetics~($\Delta G_0$)',alpha=0.7,s=60)

    if errorbar:
        yerr_upper = np.nanmax(RTlogQ[rev_rxns,:], axis=1) - RTlogQ[rev_rxns, 0]
        yerr_lower = - np.nanmin(RTlogQ[rev_rxns,:], axis=1) + RTlogQ[rev_rxns, 0]
        ax.errorbar(rev_x_coords,temp_y, yerr=[yerr_lower, yerr_upper], color=kinetics_color, linestyle='None', capsize=5)

    legend_elements.append(a)


    intervals = 5
    current_min = intervals*np.floor_divide(current_min,intervals)
    current_max = intervals*(np.floor_divide(current_max,intervals) + 1)
        
    # ax.autoscale(enable=True,axis='y',tight=True)
    ax.set_autoscaley_on(False)
    ax.set_ylim(ymin=current_min,ymax=current_max)
    
    up = current_max//(2*intervals); down = current_min//(2*intervals)+1
    yticks = 2*intervals*np.arange(down, up+1)
    ax.set_yticks(yticks)

    # draw vertical lines to separate reactions
    vline_x_coords = np.arange(1,N_rxn) - 0.5
    # ax_ylims = ax.get_ylim()
    ax.vlines(vline_x_coords,ymin=current_min,ymax=current_max,linestyles='dotted',linewidth=0.75,alpha=0.25)

    # adjust all x-ticks
    ax.set_xlim(xmin=-0.75,xmax=N_rxn-0.25)
    ax.set_xticks(x_coords)
    ax.set_xticklabels(rxn_labels, rotation=90, ha="center", fontsize=8)

    # ax.legend(handles=legend_elements,fontsize=8,ncol=2)
    ax.legend(handles=legend_elements,fontsize=8,loc='lower right',title=r'$\Delta G=\Delta G_0+RT\log Q$',title_fontsize=9)

    plt.tight_layout()

    fig.savefig(fig_path, dpi=600)

    return None



def plot_accumu_deltaG(dG0_order_kinetics, RTlogQ_order, rxn_order_list, fig_path, errorbar=False):
    # l_w_ratio = 1/1.618
    N_t = 1

    l_w_ratio = 0.39

    fig_size = [174,174*l_w_ratio]

    heights = [1.0,1.0,1.0]

    kinetics_color = 'darkcyan'
    equilibrator_color = 'darkorange'

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    # set y-axes labels
    ax.set_ylabel(r'Cumulative $\Delta G$~[kJ/mol]',fontsize=14)

    # set y-axes scales and gridlines
    ax.set_axisbelow(True)
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.tick_params(axis='y',labelsize=10)
    # ax.set_yscale('log')
    ax.yaxis.grid(which='major',linestyle='-',linewidth=1.5)
    ax.yaxis.grid(which='minor',linestyle='--')

    # set x-axes tick labels
    rxn_labels = []; N_order_rxn = len(rxn_order_list)
    for i in range(N_order_rxn):
        temp_label = rxn_order_list[i]
        temp_label = r''+temp_label.replace('_','\_')
        rxn_labels.append(temp_label)
        
    # x coordinates
    x_coords = np.arange(0,N_order_rxn)

    legend_elements = []

    dG0_order_kinetics = dG0_order_kinetics.reshape(np.shape(dG0_order_kinetics)[0],1)
    temp_y = dG0_order_kinetics + RTlogQ_order
    temp_y = np.cumsum(temp_y, axis=0) # temp_y is 2D with dimensions rxns and replicates

    current_min = np.nanmin(temp_y)
    current_max = np.nanmax(temp_y)

#     a = ax.scatter(x_coords,temp_y,color=kinetics_color,marker=time_markers[i_t],label=r'kinetics - '+time_labels[i_t],alpha=0.7)
    a = ax.scatter(x_coords,temp_y[:,0],color=kinetics_color,marker='x',label=r'kinetics',alpha=0.7,s=80)
    if errorbar:
        yerr_upper = np.nanmax(temp_y, axis=1) - temp_y[:,0]
        yerr_lower = - np.nanmin(temp_y, axis=1) + temp_y[:,0]
        ax.errorbar(x_coords,temp_y[:,0], yerr=[yerr_lower, yerr_upper], color=kinetics_color, linestyle='None', capsize=5)

    legend_elements.append(a)        

    intervals = 5
    current_min = intervals*np.floor_divide(current_min,intervals)
    current_max = intervals*(np.floor_divide(current_max,intervals) + 1)
        
    # ax.autoscale(enable=True,axis='y',tight=True)
    ax.set_autoscaley_on(False)
    ax.set_ylim(ymin=current_min,ymax=current_max)
    # ax.set_ylim(ymin=-50,ymax=50)
    # draw vertical lines to separate reactions
    vline_x_coords = np.arange(1,N_order_rxn) - 0.5
    # ax_ylims = ax.get_ylim()
    ax.vlines(vline_x_coords,ymin=current_min,ymax=current_max,linestyles='dotted',linewidth=0.75,alpha=0.25)

    # adjust all x-ticks
    ax.set_xlim(xmin=-0.75,xmax=N_order_rxn-0.25)
    ax.set_xticks(x_coords)
    ax.set_xticklabels(rxn_labels, rotation=90, ha="center", fontsize=8)

    # set x-axis label
    # ax.set_xlabel(r'Reaction',fontsize=14)

    # fig.suptitle(r'$\Delta G=\Delta G_0^{\prime}+RT\log Q$',fontsize=16)

    # ax.legend(handles=legend_elements,fontsize=8,ncol=2)
    # ax.legend(handles=legend_elements,fontsize=8)

    plt.tight_layout()

    fig.savefig(fig_path, dpi=600)

    plt.show()
    return None