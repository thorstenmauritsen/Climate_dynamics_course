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

amip = netcdf.netcdf_file('../Data/ATM_echam-6.3.02p4_amip_1979-2008_mean.nc')
amip4xCO2 = netcdf.netcdf_file('../Data/ATM_echam-6.3.02p4_amip4xCO2_1979-2008_mean.nc')

amip_diff = netcdf.netcdf_file('../Data/BOT_echam-6.3.02p4_amip4xCO2_1979-2008_timmean_fldmean_anomaly.nc')

#hf_amip = netcdf.netcdf_file('../Data/echam-6.3.02p4_amip_highfrequency_197901-2.01_echamm_fldmean.nc')
hf_amip4xCO2 = netcdf.netcdf_file('../Data/echam-6.3.02p4_amip4xCO2_highfrequency_197901-2.01_echamm_fldmean_daymean_anomaly.nc')

# ---------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(4,6))  


axes.plot(amip.variables['t'][0,:,0],amip.variables['geopoth'][0,:,0]/1000, color='black')
axes.plot(amip4xCO2.variables['t'][0,:,0],amip4xCO2.variables['geopoth'][0,:,0]/1000, color='red')

axes.set_xlabel('Global mean temperature (K)')
axes.set_ylabel(r'Height (km)')

axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
axes.spines['bottom'].set_position('zero')

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['top', 'right'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom', 'left']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

plt.tight_layout()
plt.savefig('../plots/Stratospheric_adjustment_profiles.pdf', dpi=300)
plt.close()

#-------------------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,4))  

time = hf_amip4xCO2.variables['time'][:]-hf_amip4xCO2.variables['time'][0]
#toa  = hf_amip4xCO2.variables['trad0'][:,0,0]-hf_amip.variables['trad0'][:,0,0] + hf_amip4xCO2.variables['srad0'][:,0,0]-hf_amip.variables['srad0'][:,0,0]
#toaf = hf_amip4xCO2.variables['traf0'][:,0,0]-hf_amip.variables['traf0'][:,0,0] + hf_amip4xCO2.variables['sraf0'][:,0,0]-hf_amip.variables['sraf0'][:,0,0]
#sfc  = hf_amip4xCO2.variables['trads'][:,0,0]-hf_amip.variables['trads'][:,0,0] + hf_amip4xCO2.variables['srads'][:,0,0]-hf_amip.variables['srads'][:,0,0]
#trop = hf_amip4xCO2.variables['tradl'][:,0,0]-hf_amip.variables['tradl'][:,0,0] + hf_amip4xCO2.variables['sradl'][:,0,0]-hf_amip.variables['sradl'][:,0,0]

toa  = hf_amip4xCO2.variables['trad0'][:,0,0] + hf_amip4xCO2.variables['srad0'][:,0,0]
trop = hf_amip4xCO2.variables['tradl'][:,0,0] + hf_amip4xCO2.variables['sradl'][:,0,0]
sfc  = hf_amip4xCO2.variables['trads'][:,0,0] + hf_amip4xCO2.variables['srads'][:,0,0] + hf_amip4xCO2.variables['ahfs'][:,0,0] + hf_amip4xCO2.variables['ahfl'][:,0,0] 

toaf  = hf_amip4xCO2.variables['traf0'][:,0,0] + hf_amip4xCO2.variables['sraf0'][:,0,0]
tropf = hf_amip4xCO2.variables['trafl'][:,0,0] + hf_amip4xCO2.variables['srafl'][:,0,0]
sfcf  = hf_amip4xCO2.variables['trafs'][:,0,0] + hf_amip4xCO2.variables['srafs'][:,0,0] + hf_amip4xCO2.variables['ahfs'][:,0,0] + hf_amip4xCO2.variables['ahfl'][:,0,0] 

meanf = amip_diff.variables['traf0'][:,0,0] + amip_diff.variables['sraf0'][:,0,0]
axes.plot((0,60),(meanf,meanf),color='lightgray', linestyle='--')

axes.plot(time,toaf, color='black')
axes.text(60,5,'Top-of-atmosphere',color='black',ha='right')
axes.plot(time,tropf, color='blue')
axes.text(60,4.3,'200 hPa level',color='blue',ha='right')
#axes.plot(time,sfcf, color='red')
#axes.plot(amip4xCO2.variables['t'][0,:,0],amip4xCO2.variables['geopoth'][0,:,0]/1000, color='red')

plt.ylim(0,10)

axes.set_xlabel('Time (days)')
axes.set_ylabel(r'Clear-sky downward radiation flux change (Wm$^{-2}$)')

axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
axes.spines['bottom'].set_position('zero')

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['top', 'right'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom', 'left']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

plt.tight_layout()
plt.savefig('../plots/Stratospheric_adjustment_timeseries.pdf', dpi=300)
plt.close()







#amip.close()
#amip4xCO2.close()
