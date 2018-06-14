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

So     = 1360.0
albedo = 0.29
sigma  = 5.67e-8

emissivity  = np.arange(0.3,0.71,0.01,'float')
temperature = (So*(1-albedo)/(4*emissivity*sigma))**0.25
LW  = 5.67e-8*temperature**4.
SW  = 1360./4.*(1-albedo)
LWe = emissivity*sigma*temperature**4.

lambda_planck = -4*emissivity*sigma*temperature**3.
lambda_solar   = -4*0.61*sigma*temperature**3.
test = -4*sigma*255**3*emissivity**(1/4)

meraner_T        = np.array([2.79,6.70,12.46,22.68])+287.0
meraner_lambda_T = np.array([-4.27,-4.11,-4.30,-4.28])

#

fig, axes = plt.subplots(1,1, figsize=(6,4))  

hfixed, = axes.plot(temperature, lambda_planck, color='blue', label='Variable emissivity')
#axes.plot(temperature, test, color='green')
hvariable, = axes.plot(temperature, lambda_solar, color='black', label='Fixed emissivity')
hmeraner = axes.scatter(meraner_T, meraner_lambda_T, color='red', label=r'ECHAM6 temperature feedback: 2, 4, 8 and 16xCO$_2$')

axes.legend(handles=[hmeraner,hfixed,hvariable], fontsize=10, loc=2, frameon=False)

#axes.text(280,-1.0,r'ECHAM6 temperature feedback: 2, 4, 8 and 16xCO$_2$', color='red')
#axes.text(280,-1.4,'Solar forced, fixed emissivity', color='black')
#axes.text(280,-1.8,'Emissivity forced', color='blue')


axes.text(temperature[0], lambda_planck[0]+0.1, r'$\epsilon = $ '+str(emissivity[0]), color='blue',va='bottom', ha='center')
axes.text(temperature[-1], lambda_planck[-1]-0.1, r'$\epsilon = $ '+str(emissivity[-1]), color='blue',va='top', ha='center')


axes.set_xlabel('Surface temperature (K)')
axes.set_ylabel(r'Feedback parameter (Wm$^{-2}$K$^{-1}$)')

plt.ylim((-6,0))
plt.xlim((275,349))

axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('top')
axes.spines['top'].set_position('zero')
axes.xaxis.set_label_position("top")

spines_to_remove        = ['bottom', 'right'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'top', 'left']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

    
plt.tight_layout()
plt.savefig('../plots/Planck_state_dependence.pdf', dpi=600)
plt.close()

