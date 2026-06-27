import AC_consensus
import AC_consensus_plots

from quorums.quorum_sizes import *
import matplotlib.pyplot as plt
from gmpy2 import *


import matplotlib as mpl
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#1f77b4', "r", "#09b000", "#f30bd0", 'b'])
 

ctx = get_context()
ctx.precision = 1000


N = (10**9) 
safe_over_live = lambda x, y: x if x > y else 1
max_tail = lambda x,y: max(x,y)


configs = [(N, 0.3, 0.0)]
AC_consensus_plots.plot_fM(configs, 128, cod=False, n_steps=210,  measure=max_tail, plot_cutoff=True, save_fig=True, save_name="M30_AC_MAX")

configs = [(N, 0.3, 0.0), (N, 0.15, 0.0), (N, 0.01, 0.0)]
#AC_consensus_plots.plot_fM(configs, 128, cod=False, n_steps=210, measure=max_tail, split=False, plot_cutoff=True, save_fig=True, save_name="M1_15_30_AC_MAX")

configs = [(N, 0.1, 0.2)]
#AC_consensus_plots.plot_fI(configs, 128, n_steps=210, cod=False, measure=max_tail, split=False, plot_cutoff=True, save_fig=True, save_name="I20_M10_AC_MAX")
#AC_consensus_plots.plot_fI(configs, 128, n_steps=210, cod=False, measure=max_tail, split=True, plot_cutoff=False, save_fig=True, save_name="I20_M10_AC_SPLIT")

configs = [(N, 0.01, 0.2)]
#AC_consensus_plots.plot_fM(configs, 256, cod=True, n_steps=210, measure=safe_over_live, split=False, plot_cutoff=True, save_fig=True, save_name="M1_COD")
#AC_consensus_plots.plot_fM(configs, 256, cod=True, n_steps=210, measure=safe_over_live, split=True, plot_cutoff=True, save_fig=True, save_name="M1_COD_SPLIT")

configs = [(N, 0.01, 0.2), (N, 0.1, 0.2), (N, 0.2, 0.2)]
#AC_consensus_plots.plot_fM(configs, 256, cod=True, n_steps=210, measure=safe_over_live, plot_cutoff=True, save_fig=True, save_name="M1_M10_M20_COD")

configs = [(N, 0.1, 0.2)]
#AC_consensus_plots.plot_fI(configs, 256, cod=True, n_steps=210, measure=safe_over_live, save_fig=True, save_name="I20_M10_COD")
#AC_consensus_plots.plot_fI(configs, 256, cod=True, n_steps=210, measure=safe_over_live, split=True, plot_cutoff=True, save_fig=True, save_name="I20_M10_COD_SPLIT")













