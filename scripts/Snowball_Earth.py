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

temperature = np.arange(200,320,1,'float')
albedo = np.empty_like(temperature)

ai = 0.7
ao = 0.29
ti = 250.
to = 288.

for t in temperature:
    i = np.where(temperature == t)
    if t <= ti:
        albedo[i] = ai
    elif t >= to:
        albedo[i] = ao
    else:
        albedo[i] = ao + (ai-ao)*(t-to)**2./(ti-to)**2.

LW  = 5.67e-8*temperature**4.
SW  = 1360./4.*(1-albedo)
LWe = 0.61*5.67e-8*temperature**4.




#

fig, axes = plt.subplots(1,1, figsize=(6,4))  

axes.plot(temperature, SW, color='blue')
axes.plot(temperature[3:], LW[3:], color='black')
axes.plot(temperature[15:], LWe[15:], color='red')

axes.text(270,350,r'$\sigma T^4$', color='black')
axes.text(293,300,r'$\epsilon \sigma T^4$', color='red')
axes.text(300,210,r'$\frac{S_o}{4}(1-\alpha)$', color='blue')

#axes.text(270,150,str(ai)+' , T$_s < $ '+str(ti) \n str(a0) ', T$_s < $ '+str(ti) )



#axes.fill_between(sinlats, sw_all, -lw_all, where=-lw_all >= sw_all, facecolor='blue', interpolate=True)
#axes.fill_between(sinlats, sw_all, -lw_all, where=-lw_all < sw_all, facecolor='red', interpolate=True)
#axes.fill_between(sinlats, y1, y2, where=y2 <= y1, facecolor='red', interpolate=True)
#axes.plot(sinlats,  sw_all, color='black')
#axes.plot(sinlats, -lw_all, color='black')
#axes.text(0,15,'CERES EBAF Edition 2.8, 2001-2014', color=almost_black, ha='center')
#axes.text(0,330,'Shortwave', color=almost_black, ha='center')
#axes.text(0,220,'-Longwave', color=almost_black, ha='center')

axes.set_xlabel('Surface temperature (K)')
axes.set_ylabel(r'Top-of-atmosphere net fluxes (Wm$^{-2}$)')

plt.ylim((50,400))
#plt.xticks(np.arange(-1.0,1.5,0.5))
#axes.set_xticklabels(('90S','30S','Equator','30N','90N'))



axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
axes.spines['bottom'].set_position('zero')


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
plt.savefig('../plots/Snowball_Earth_budget.pdf', dpi=600)
plt.close()

