import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import math

# set general text options for plotting
plt.rc('text'           , usetex=True)
plt.rc('font'           , size=12)
plt.rc('legend'         , fontsize=12)
plt.rc('text.latex'     , preamble=r'\usepackage{cmbright}')

#

almost_black            = '#262626'

f2x= 3.7
gamma = 1.0

epsilon  = np.arange(1.0,2.0,0.01,'float')

ECS = -f2x/(-f2x/1.78+gamma*(epsilon-1))

# ---------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,4))  

axes.plot(epsilon, ECS, color='blue')

axes.set_xlabel(r'Ocean heat uptake efficacy, $\epsilon$')
axes.set_ylabel('ECS')

plt.ylim((0,4))

axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
axes.spines['bottom'].set_position('zero')
axes.xaxis.set_label_position("bottom")

spines_to_remove        = ['top', 'right'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom', 'left']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)
    
plt.tight_layout()
plt.savefig('../plots/ECS_time-dependent_feedback.pdf', dpi=600)
plt.close()

