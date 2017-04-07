import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
import math

# set general text options for plotting
plt.rc('text'           , usetex=True)
plt.rc('font'           , size=12)
plt.rc('legend'         , fontsize=12)
plt.rc('text.latex'     , preamble=r'\usepackage{cmbright}')

#

almost_black            = '#262626'

f = netcdf.netcdf_file('../Data/CERES_EBAF-TOA_Ed2.8_Subset_2001-2014_timmean_zonmean.nc')
g = netcdf.netcdf_file('../Data/CERES_EBAF-Surface_Ed2.8_Subset_2001-2014_timmean_zonmean.nc')

lats    = f.variables['lat'][:]
sinlats = np.sin(lats*math.pi/180)

sw_all =  f.variables['solar_mon'][0,:,0] - f.variables['toa_sw_all_mon'][0,:,0]
albedo =  f.variables['toa_sw_all_mon'][0,:,0]/f.variables['solar_mon'][0,:,0]
lw_all = -f.variables['toa_lw_all_mon'][0,:,0]

#

fig, axes = plt.subplots(1,1, figsize=(6,4))  

axes.fill_between(sinlats, sw_all, -lw_all, where=-lw_all >= sw_all, facecolor='blue', interpolate=True)
axes.fill_between(sinlats, sw_all, -lw_all, where=-lw_all < sw_all, facecolor='red', interpolate=True)
#axes.fill_between(sinlats, y1, y2, where=y2 <= y1, facecolor='red', interpolate=True)
axes.plot(sinlats,  sw_all, color='black')
axes.plot(sinlats, -lw_all, color='black')
axes.text(0,15,'CERES EBAF Edition 2.8, 2001-2014', color=almost_black, ha='center')
axes.text(0,330,'Shortwave', color=almost_black, ha='center')
axes.text(0,220,'-Longwave', color=almost_black, ha='center')

#axes.set_xlabel('Latitude')
axes.set_ylabel(r'Top-of-atmosphere net flux (Wm$^{-2}$)')

#plt.ylim((0,350))
plt.xticks(np.arange(-1.0,1.5,0.5))
axes.set_xticklabels(('90S','30S','Equator','30N','90N'))



#axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
axes.spines['bottom'].set_position('zero')

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['top'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom', 'left', 'right']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

plt.tight_layout()
plt.savefig('../plots/CERES_Ebaf_zonalmean.pdf', dpi=600)
plt.close()

#-------------------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,4))  

axes.plot(sinlats,  albedo, color='black')
#axes.plot(sinlats,  g.variables['sfc_sw_up_all_mon'][0,:,0]/g.variables['sfc_sw_down_all_mon'][0,:,0], color='blue')

axes.set_ylabel('Planetary albedo')

plt.ylim((0,0.8))
plt.xticks(np.arange(-1.0,1.5,0.5))
axes.set_xticklabels(('90S','30S','Equator','30N','90N'))

#axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
#axes.spines['bottom'].set_position('zero')

axes.set_yticks(np.arange(0.1,0.9,0.1))

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['top'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom', 'left', 'right']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

plt.tight_layout()
plt.savefig('../plots/CERES_Ebaf_albedo.pdf', dpi=600)
plt.close()


#-------------------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,4))  

axes.plot(sinlats,  f.variables['solar_mon'][0,:,0], color='black')
#axes.plot(sinlats,  sw_all, color='blue')

surface_down = g.variables['sfc_sw_down_all_mon'][0,:,0]
surface_up   = g.variables['sfc_sw_up_all_mon'][0,:,0]

axes.fill_between(sinlats, surface_down, surface_up, where= surface_up < surface_down, facecolor='lightgray', interpolate=True)
axes.plot(sinlats,  surface_down, color='red')
axes.plot(sinlats,  surface_up,   color='blue')

axes.set_ylabel('Shortwave fluxes')

#plt.ylim((0,0.8))
plt.xticks(np.arange(-1.0,1.5,0.5))
axes.set_xticklabels(('90S','30S','Equator','30N','90N'))

#axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
#axes.spines['bottom'].set_position('zero')

#axes.set_yticks(np.arange(0.1,0.9,0.1))

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['top'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom', 'left', 'right']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

plt.tight_layout()
plt.savefig('../plots/CERES_Ebaf_shortwave.pdf', dpi=600)
plt.close()







f.close()
g.close()
