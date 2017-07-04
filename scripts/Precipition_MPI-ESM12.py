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


abrupt2xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.00_abrupt2xCO2_anomaly.nc')
abrupt4xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.00_abrupt4xCO2_anomaly.nc')
abrupt8xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.00_abrupt8xCO2_anomaly.nc')
abrupt16xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.00_abrupt16xCO2_anomaly.nc')



#abrupt2xCO2_prp = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt2xCO2_echam6_prp_mean_anomalies.nc')
#abrupt2xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt2xCO2_echam6_BOT_mean_anomalies.nc')
#abrupt4xCO2_prp = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt4xCO2_echam6_prp_mean_anomalies.nc')
#abrupt4xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt4xCO2_echam6_BOT_mean_anomalies.nc')
#abrupt8xCO2_prp = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt8xCO2_echam6_prp_mean_anomalies.nc')
#abrupt8xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt8xCO2_echam6_BOT_mean_anomalies.nc')
#abrupt16xCO2_prp = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt16xCO2_echam6_prp_mean_anomalies.nc')
#abrupt16xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt16xCO2_echam6_BOT_mean_anomalies.nc')

#abrupt16xCO2 = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt16xCO2_echam6_Chris_anomalies.nc')

#amip4xCO2 = netcdf.netcdf_file('../Data/BOT_echam-6.3.02p4_amip4xCO2_1979-2008_timmean_fldmean_anomaly.nc')

#abrupt2xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.00_abrupt2xCO2_anomaly.nc')

#-------------------------------------------------------------

id=np.arange(0,999)
#idfit=np.arange(99,999)

control_t = 287.63

T_2xCO2 = abrupt2xCO2_bot.variables['tsurf'][id,0,0]
T_4xCO2 = abrupt4xCO2_bot.variables['tsurf'][id,0,0]
T_8xCO2 = abrupt8xCO2_bot.variables['tsurf'][id,0,0]
T_16xCO2 = abrupt16xCO2_bot.variables['tsurf'][id,0,0]

mmday = 24*3600
control_p = 3.3007e-05*mmday

P_2xCO2 = mmday*abrupt2xCO2_bot.variables['precip'][id,0,0]
P_4xCO2 = mmday*abrupt4xCO2_bot.variables['precip'][id,0,0]
P_8xCO2 = mmday*abrupt8xCO2_bot.variables['precip'][id,0,0]
P_16xCO2 = mmday*abrupt16xCO2_bot.variables['precip'][id,0,0]

control_qvi = 22.894

qvi_2xCO2 = abrupt2xCO2_bot.variables['qvi'][id,0,0]+control_qvi
qvi_4xCO2 = abrupt4xCO2_bot.variables['qvi'][id,0,0]+control_qvi
qvi_8xCO2 = abrupt8xCO2_bot.variables['qvi'][id,0,0]+control_qvi
qvi_16xCO2 = abrupt16xCO2_bot.variables['qvi'][id,0,0]+control_qvi

# ---------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,6))  

color2                  = 'royalblue'
color4                  = "deepskyblue"
color8                  = 'darkkhaki'
color16                 = 'salmon'

q2xCO2 = axes.scatter(T_2xCO2, qvi_2xCO2,color=color2,label=r'2xCO$_2$')
q4xCO2 = axes.scatter(T_4xCO2, qvi_4xCO2,color=color4,label=r'4xCO$_2$')
q8xCO2 = axes.scatter(T_8xCO2, qvi_8xCO2,color=color8,label=r'8xCO$_2$')
q16xCO2 = axes.scatter(T_16xCO2, qvi_16xCO2,color=color16,label=r'16xCO$_2$')


dt=np.arange(0,18,0.1)

alpha=2.5e6/461./273.15**2

line1, = axes.plot(dt,alpha*control_qvi*dt+control_qvi,color='black',linestyle='--', label=r'Linear, $\alpha = $'+str(round(alpha,3)))
line2, = axes.plot(dt,control_qvi*np.exp(alpha*dt),color='black', label=r'Exponential, $\alpha = $'+str(round(alpha,3)))
#line3, = axes.plot(dt,alpha*control_qvi*dt+control_qvi,color='black',linestyle=':', label=r'Linear, $\alpha = $'+str(round(alpha,3)))

control = axes.scatter(0,control_qvi,color='black', label='Control')

plt.legend(handles=[control, q2xCO2, q4xCO2, q8xCO2, q16xCO2,line1, line2], loc=4, frameon=False)

#plt.legend(handles=[line1, line2]) #, loc=4, frameon=False)

plt.xlim(-2,18)
plt.ylim(0,100)
#plt.ylim(-0.4,0.8)

axes.set_ylabel(r'Vertically integrated water vapor (gm$^{-2}$)')
axes.set_xlabel(r'Global mean temperature change (K)')

axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
axes.spines['bottom'].set_position('zero')
axes.spines['left'].set_position('zero')

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
plt.savefig('../plots/MPI-ESM12_water_vapor.pdf', dpi=300)
plt.close()


# ---------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,10))  

color2                  = 'royalblue'
color4                  = "deepskyblue"
color8                  = 'darkkhaki'
color16                 = 'salmon'

axes.scatter(T_2xCO2, P_2xCO2+control_p,color=color2)
axes.scatter(T_4xCO2, P_4xCO2+control_p,color=color4)
axes.scatter(T_8xCO2, P_8xCO2+control_p,color=color8)
axes.scatter(T_16xCO2, P_16xCO2+control_p,color=color16)

axes.plot(dt, control_p*alpha*dt+control_p ,color='black',linestyle='--')
axes.plot(dt, control_p*np.exp(alpha*dt) ,color='black')
axes.scatter(0,control_p,color='black')

plt.xlim(-2,18)
#plt.ylim(-5e-6,1e-5)
plt.ylim(-0.0,5.0)

plt.legend(handles=[control, q2xCO2, q4xCO2, q8xCO2, q16xCO2,line1,line2], loc=4, frameon=False)

axes.set_xlabel(r'Global mean temperature change, $T_s$ (K)')
axes.set_ylabel(r'Global mean precipitation, $P$ (mm day$^{-1}$)')

axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')
axes.spines['bottom'].set_position('zero')
axes.spines['left'].set_position('zero')

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
plt.savefig('../plots/MPI-ESM12_precipitation.pdf', dpi=300)
plt.close()

# ---------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,6))  

color2                  = 'royalblue'
color4                  = "deepskyblue"
color8                  = 'darkkhaki'
color16                 = 'salmon'

wm2tommperday = 24*3600/2.5e6

axes.scatter(T_2xCO2, -wm2tommperday*(abrupt2xCO2_bot.variables['trad0'][id,0,0]-abrupt2xCO2_bot.variables['trads'][id,0,0]),color=color2)
axes.scatter(T_4xCO2, -wm2tommperday*(abrupt4xCO2_bot.variables['trad0'][id,0,0]-abrupt4xCO2_bot.variables['trads'][id,0,0]),color=color4)
axes.scatter(T_8xCO2, -wm2tommperday*(abrupt8xCO2_bot.variables['trad0'][id,0,0]-abrupt8xCO2_bot.variables['trads'][id,0,0]),color=color8)
thermal = axes.scatter(T_16xCO2, -wm2tommperday*(abrupt16xCO2_bot.variables['trad0'][id,0,0]-abrupt16xCO2_bot.variables['trads'][id,0,0]),color=color16, label='Thermal emission')

axes.scatter(T_2xCO2, -wm2tommperday*(abrupt2xCO2_bot.variables['srad0'][id,0,0]-abrupt2xCO2_bot.variables['srads'][id,0,0]),color=color2,marker='^')
axes.scatter(T_4xCO2, -wm2tommperday*(abrupt4xCO2_bot.variables['srad0'][id,0,0]-abrupt4xCO2_bot.variables['srads'][id,0,0]),color=color4,marker='^')
axes.scatter(T_8xCO2, -wm2tommperday*(abrupt8xCO2_bot.variables['srad0'][id,0,0]-abrupt8xCO2_bot.variables['srads'][id,0,0]),color=color8,marker='^')
shortwave = axes.scatter(T_16xCO2, -wm2tommperday*(abrupt16xCO2_bot.variables['srad0'][id,0,0]-abrupt16xCO2_bot.variables['srads'][id,0,0]),color=color16,marker='^', label='Solar absorption')

axes.scatter(T_2xCO2, wm2tommperday*abrupt2xCO2_bot.variables['ahfs'][id,0,0],color=color2,marker='x')
axes.scatter(T_4xCO2, wm2tommperday*abrupt4xCO2_bot.variables['ahfs'][id,0,0],color=color4,marker='x')
axes.scatter(T_8xCO2, wm2tommperday*abrupt8xCO2_bot.variables['ahfs'][id,0,0],color=color8,marker='x')
sensible = axes.scatter(T_16xCO2, wm2tommperday*abrupt16xCO2_bot.variables['ahfs'][id,0,0],color=color16,marker='x', label='Sensible heat flux')

plt.legend(handles=[thermal, shortwave, sensible], loc=2, frameon=False)

plt.xlim(0,18)
#plt.ylim(-5e-6,1e-5)
#plt.ylim(-0.4,0.8)

axes.set_xlabel(r'Global mean temperature change, $T_s$ (K)')
axes.set_ylabel(r'Equivalent precipitation change, (mm day$^{-1}$)')

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
plt.savefig('../plots/MPI-ESM12_atmosphere_budget.pdf', dpi=300)
plt.close()
