
from AC_consensus import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

from copy import deepcopy

def cutoff_fM_TH_30(axins, Qp, xs, ys):
    x_ticks = []
    Q, p = Qp
    x_cutoff_l, _, bit_sec = find_cutoff_fM([10**9,0.3,0.0], 128, 80, 0.01, input=(Q, p))
    x, y = xs, ys
    x.append(np.float64(x_cutoff_l))
    y.append(bit_sec)
    axins.scatter(x,y)
    x_ticks.append(float(x_cutoff_l))

    axins.set_xlim(0.066, 0.067)
    axins.axhline(y=80, color='lightgray', linestyle='--')
    axins.set_xticks(x_ticks)
    axins.set_ylim(79, 81)
    axins.set_yticks([80])

def cutoff_fM_TH_multiple(configs, axins, Qps, xs, ys):
    x_ticks = []
    for i, config in enumerate(configs):
        x_cutoff_l, _, bit_sec = find_cutoff_fM(config, 128, 80, 0.01, input=Qps[i])
        x, y = xs[i], ys[i]
        x.append(np.float64(x_cutoff_l))
        y.append(bit_sec)
        axins.scatter(x,y)
        x_ticks.append(float(x_cutoff_l))
    axins.set_xlim(0.065, 0.088)
    axins.axhline(y=80, color='lightgray', linestyle='--')
    axins.set_xticks(x_ticks)
    axins.set_ylim(79.75, 80.25)
    axins.set_yticks([80])

def cutoff_fI_TH(config, ax, Qp, x, y):
    Q, p = Qp
    axins = inset_axes(ax, width="20%", height="20%", loc='center left', bbox_to_anchor=(0.07, 0.15, 1, 1), bbox_transform=ax.transAxes)

    _, x_cutoff_left, bit_sec_left = find_cutoff_fI_left(config, 128, 80, -0.1, input=(Q,p))
    x_cutoff_right, _, bit_sec_right = find_cutoff_fI_right(config, 128, 80, 0.02, input=(Q,p))

    x.append(x_cutoff_left)
    y.append(bit_sec_left)
    axins.scatter(x, y)#, s=size, color='r')
    axins.set_xlim(-0.075, -0.055)
    axins.axhline(y=80, color='lightgray', linestyle='--')
    axins.set_xticks([np.float64(x_cutoff_left)])
    axins.set_ylim(79, 81)
    axins.set_yticks([80])

    
    axins = inset_axes(ax, width="20%", height="20%", loc='center right', bbox_to_anchor=(0, 0.15, 1, 1), bbox_transform=ax.transAxes)

    x[-1] = x_cutoff_right
    y[-1] = bit_sec_right
    axins.scatter(x, y) #, s=size)

    axins.set_xlim(0.032, 0.052)
    axins.axhline(y=80, color='lightgray', linestyle='--')
    axins.set_xticks([np.float64(x_cutoff_right)])
    axins.set_ylim(79, 81)
    axins.set_yticks([80])

def cutoff_fM_COD_1(ax, Qp, xs, ys):
    axins = inset_axes(ax, width="20%", height="20%", loc='center right', bbox_to_anchor=(0, 0.19, 1, 1), bbox_transform=ax.transAxes)
    x_ticks = []
    Q, p = Qp
    mes = lambda x, y: x if x > y else 1
    x_cutoff_l, _, bit_sec = find_cutoff_fM([10**9, 0.01, 0.2], 256, 128, 0.01, cod=True, measure=mes, input=(Q, p))
    x, y = xs, ys
    x.append(np.float64(x_cutoff_l))
    y.append(bit_sec)
    axins.scatter(x,y)
    x_ticks.append(float(x_cutoff_l))

    axins.set_xlim(0.0684, 0.0690)
    axins.axhline(y=128, color='lightgray', linestyle='--')
    axins.set_xticks(x_ticks)
    axins.set_ylim(127, 129)
    axins.set_yticks([128])


def cutoff_fM_COD_1_split(ax, Qp, xs, ys):
    axins = inset_axes(ax, width="30%", height="20%", loc='center right', bbox_to_anchor=(0, 0.19, 1, 1), bbox_transform=ax.transAxes)
    x_ticks = []
    Q, p = Qp
    mes = lambda x, y: x if x > y else 1
    x_cutoff_live, _, bit_sec_live = find_cutoff_fM([10**9, 0.01, 0.2], 256, 128, 0.01, cod=True, measure=lambda x, y: x, input=(Q, p))
    x, y = xs, ys[0]
    x.append(np.float64(x_cutoff_live))
    y.append(bit_sec_live)
    axins.scatter(x,y)
    x_ticks.append(float(x_cutoff_live))

    x_cutoff_safe, _, bit_sec_safe = find_cutoff_fM([10**9, 0.01, 0.2], 256, 128, 0.01, cod=True, measure=lambda x, y: y, input=(Q, p))
    x, y = xs, ys[1]
    x[-1] = np.float64(x_cutoff_safe)
    y.append(bit_sec_safe)
    axins.scatter(x,y)
    x_ticks.append(float(x_cutoff_safe))


    axins.set_xlim(0.064, 0.1)
    axins.axhline(y=128, color='lightgray', linestyle='--')
    axins.set_xticks(x_ticks)
    axins.set_ylim(127, 129)
    axins.set_yticks([128])

def cutoff_fM_cod_multiple(configs, axins, Qps, xs, ys):
    x_ticks = []
    for i, config in enumerate(configs):
        x_cutoff_l, _, bit_sec = find_cutoff_fM(config, 256, 128, 0.01, cod=True, input=Qps[i])
        x, y = xs[i], ys[i]
        x.append(np.float64(x_cutoff_l))
        y.append(bit_sec)
        axins.scatter(x,y)
        x_ticks.append(float(x_cutoff_l))
    axins.set_xlim(0.01, 0.075)
    axins.axhline(y=128, color='lightgray', linestyle='--')
    axins.set_xticks(x_ticks)
    axins.set_ylim(127, 129)
    axins.set_yticks([128])

def cutoff_fI_cod(config, ax, Qp, x, ys):
    Q, p = Qp
    axins = inset_axes(ax, width="20%", height="20%", loc='center left', bbox_to_anchor=(0.07, 0.15, 1, 1), bbox_transform=ax.transAxes)

    _, x_cutoff_left, bit_sec_left = find_cutoff_fI_left(config, 256, 128, -0.1, measure=lambda x, y: y, cod=True, input=(Q,p))
    x_cutoff_right, _, bit_sec_right = find_cutoff_fI_right(config, 256, 128, 0.02, measure=lambda x, y: x, cod=True, input=(Q,p))

    x.append(x_cutoff_left)
    ys[0].append(bit_sec_left)
    axins.scatter(x, ys[0], color='r')#, s=size, color='r')
    axins.set_xlim(-0.075, -0.055)
    axins.axhline(y=128, color='lightgray', linestyle='--')
    axins.set_xticks([np.float64(x_cutoff_left)])
    axins.set_ylim(127, 129)
    axins.set_yticks([128])

    
    axins = inset_axes(ax, width="20%", height="20%", loc='center right', bbox_to_anchor=(0, 0.15, 1, 1), bbox_transform=ax.transAxes)

    x[-1] = x_cutoff_right
    ys[1].append(bit_sec_right)
    axins.scatter(x, ys[1]) #, s=size)

    axins.set_xlim(0.032, 0.052)
    axins.axhline(y=128, color='lightgray', linestyle='--')
    axins.set_xticks([np.float64(x_cutoff_right)])
    axins.set_ylim(127, 129)
    axins.set_yticks([128])



def plot_cutoff_fM(configs, ax, Qps, xs, ys, cod, split=False):
    if cod:
        if len(configs) == 1:
            fM, fI = configs[0][1], configs[0][2]
            if (fM, fI) == (0.01, 0.2):
                cutoff_fM_COD_1_split(ax, Qps[0], xs[0], ys[0]) if split else cutoff_fM_COD_1(ax, Qps[0], xs[0], ys[0]) 
        else:
            axins = inset_axes(ax, width="35%", height="20%", loc='center right', bbox_to_anchor=(0, 0.19, 1, 1), bbox_transform=ax.transAxes)
            cutoff_fM_cod_multiple(configs, axins, Qps, xs, ys)
            
        
    else:
        if len(configs) == 1:
            axins = inset_axes(ax, width="20%", height="20%", loc='center right', bbox_to_anchor=(0, 0.19, 1, 1), bbox_transform=ax.transAxes)
            fM, fI = configs[0][1], configs[0][2]
            if (fM, fI) == (0.3, 0):
                cutoff_fM_TH_30(axins, Qps[0], xs[0], ys[0])
        else:
            axins = inset_axes(ax, width="35%", height="20%", loc='center right', bbox_to_anchor=(0, 0.19, 1, 1), bbox_transform=ax.transAxes)
            cutoff_fM_TH_multiple(configs, axins, Qps, xs, ys)

def plot_cutoff_fI(config, ax, Qp, xs, ys, cod, split):

    if cod:
        if split: 
            cutoff_fI_cod(config, ax, Qp, xs, ys)
        else:
            return -1
    
    else:
        cutoff_fI_TH(config, ax, Qp, xs, ys)





def plot_fM(configs, b, cod=False, n_steps=35, split=False,  measure=lambda x,y: max(x,y), plot_cutoff=False, save_fig=False, save_name=""):
    fig, ax = plt.subplots()
    fig.set_dpi(150)
    xs, ys = [], []
    qsps = []
    
    for config in configs:
        N, fM, fI = config
        x,y, Q, p = func_fM(N, fM, fI, b, n_steps, cod=cod, measure=measure)
        print(f'Number of plotted points = {len(x)}')
        xs.append(x)
        qsps.append((Q,p))
        
        if split:
            y1_l = list(map(lambda z: -gmpy2.log2(z[0]), y))
            y1_r = list(map(lambda z: -gmpy2.log2(z[1]), y))
            ys.append((y1_l, y1_r))

            ax.scatter(x, y1_l, label=f'Liveness', s=2)
            ax.scatter(x, y1_r, label=f'Safety', s=2)
        else:
            y1 = list(map(lambda z: -gmpy2.log2(measure(z[0],z[1])), y))

            ys.append(y1)
            ax.scatter(x, y1, label=fr'$f_M = {fM}$, $f_I = {fI}$, $Q = {Q}$', s=2 )
        

        
    ax.legend(loc=1)
    legend = ax.get_legend()
    for handle in legend.legend_handles:
        handle.set_sizes([24.0])
    
    ax.set_xlabel(r"$Δf_M = f_M' - f_M$")
    ax.set_ylabel(u"Bit Security")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.axhline(y=60, color='lightgray', linestyle='--')
    ax.axhline(y=80, color='lightgray', linestyle='--')
    ax.axhline(y=128, color='lightgray', linestyle='--')
    ax.axhline(y=256, color='lightgray', linestyle='--')
    ax.axhline(y=512, color='lightgray', linestyle='--')
    
    if b >= 256:
        ax.set_yticks([32, 60, 80, 128, 180, 256])
        ax.set_ylim((0, 300))
    else:
        ax.set_yticks([32, 60, 80, 128])
        ax.set_ylim((50, 140))

    ax.set_xticks(np.linspace(0, .34, 35))
    ax.tick_params(axis='x', labelsize=8)
    ax.set_xlim((0, 0.13))


    if plot_cutoff:
        plot_cutoff_fM(configs, ax, qsps, xs, ys, cod, split=split)

        
        
        

    

    if save_fig:
        plt.savefig("images/" + save_name, dpi=600)
    else:
        plt.show()


def plot_fI(configs, b, n_steps=31, cod=False, split=False, measure=lambda x,y: max(x,y), plot_cutoff=False, save_fig=False, save_name=""):
    fig, ax = plt.subplots()
    #fig.set_size_inches(11.0, 8.0)
    fig.set_dpi(150)
    for config in configs:
        N, fM, fI = config 
        
        x, y, Q, p = func_fI(N, fM, fI, b, n_steps, cod,  measure)
        
        print(f'Number of plotted points = {len(x)}')
        


        if split:
            y_r = list(map(lambda z: -gmpy2.log2(z[0]), y))
            y_l = list(map(lambda z: -gmpy2.log2(z[1]), y))
            ax.scatter(x, y_r, label=f'Liveness', s=2)
            ax.scatter(x, y_l, label=f'Safety', s=2)
            if plot_cutoff:
                plot_cutoff_fI(config, ax, (Q,p), x, (y_r, y_l), cod, split=split)

            
        else:
            y1 = list(map(lambda z: -gmpy2.log2(measure(z[0],z[1])), y))
            ax.scatter(x, y1, s=2, label=f'$f_I = {fI}, f_M = {fM}, Q = {Q}$')
            
            if plot_cutoff:
                plot_cutoff_fI(config, ax, (Q,p), x, y1, cod, split)


        


    if b >= 256:
        ax.set_yticks([32, 60, 80, 128, 256])
        ax.set_ylim((0, 300))
    else:
        ax.set_yticks([0, 32, 60, 80, 100, 128, 150])
        ax.set_ylim((20, 160))


    ax.legend(loc=1)
    legend = ax.get_legend()
    for handle in legend.legend_handles:
        handle.set_sizes([24.0])
    ax.set_xlabel(u"$Δf_I = f_I' - f_I$")
    ax.set_ylabel(u"Bit Security")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.axhline(y=60, color='lightgray', linestyle='--')
    ax.axhline(y=80, color='lightgray', linestyle='--')
    ax.axhline(y=128, color='lightgray', linestyle='--')
    ax.axhline(y=256, color='lightgray', linestyle='--')
    ax.axhline(y=512, color='lightgray', linestyle='--')

    x = np.linspace(-0.15, 0.15, 31)
    ax.set_xticks(x[0::3])
    ax.tick_params(axis='x', labelsize=9)
    ax.set_xlim((-0.15, 0.15))

    

    if save_fig:
        plt.savefig("images/" + save_name, dpi=600)
    else:
        plt.show()