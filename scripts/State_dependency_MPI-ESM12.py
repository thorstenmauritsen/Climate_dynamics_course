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


abrupt2xCO2_prp = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt2xCO2_echam6_prp_mean_anomalies.nc')
abrupt2xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt2xCO2_echam6_BOT_mean_anomalies.nc')
abrupt4xCO2_prp = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt4xCO2_echam6_prp_mean_anomalies.nc')
abrupt4xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt4xCO2_echam6_BOT_mean_anomalies.nc')
abrupt8xCO2_prp = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt8xCO2_echam6_prp_mean_anomalies.nc')
abrupt8xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt8xCO2_echam6_BOT_mean_anomalies.nc')
abrupt16xCO2_prp = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt16xCO2_echam6_prp_mean_anomalies.nc')
abrupt16xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt16xCO2_echam6_BOT_mean_anomalies.nc')

amip4xCO2 = netcdf.netcdf_file('../Data/BOT_echam-6.3.02p4_amip4xCO2_1979-2008_timmean_fldmean_anomaly.nc')

# ---------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(5,7))  

id=np.arange(0,299)
idfit=np.arange(0,299)

#N = bot.variables['trad0'][id,0,0]+bot.variables['srad0'][id,0,0]
#T = bot.variables['tsurf'][id,0,0]

T_2xCO2 = abrupt2xCO2_bot.variables['tsurf'][id,0,0]
T_4xCO2 = abrupt4xCO2_bot.variables['tsurf'][id,0,0]
T_8xCO2 = abrupt8xCO2_bot.variables['tsurf'][id,0,0]
T_16xCO2 = abrupt16xCO2_bot.variables['tsurf'][id,0,0]

N_tmp_2xCO2 = abrupt2xCO2_prp.variables['dR_tmp_trad0'][id,0,0]  + abrupt2xCO2_prp.variables['dR_tmp_srad0'][id,0,0]
N_tmp_4xCO2 = abrupt4xCO2_prp.variables['dR_tmp_trad0'][id,0,0] + abrupt4xCO2_prp.variables['dR_tmp_srad0'][id,0,0]
N_tmp_8xCO2 = abrupt8xCO2_prp.variables['dR_tmp_trad0'][id,0,0] + abrupt8xCO2_prp.variables['dR_tmp_srad0'][id,0,0]
N_tmp_16xCO2 = abrupt16xCO2_prp.variables['dR_tmp_trad0'][id,0,0] + abrupt16xCO2_prp.variables['dR_tmp_srad0'][id,0,0]

N_vap_2xCO2 = abrupt2xCO2_prp.variables['dR_vap_trad0'][id,0,0]  + abrupt2xCO2_prp.variables['dR_vap_srad0'][id,0,0]
N_vap_4xCO2 = abrupt4xCO2_prp.variables['dR_vap_trad0'][id,0,0] + abrupt4xCO2_prp.variables['dR_vap_srad0'][id,0,0]
N_vap_8xCO2 = abrupt8xCO2_prp.variables['dR_vap_trad0'][id,0,0] + abrupt8xCO2_prp.variables['dR_vap_srad0'][id,0,0]
N_vap_16xCO2 = abrupt16xCO2_prp.variables['dR_vap_trad0'][id,0,0] + abrupt16xCO2_prp.variables['dR_vap_srad0'][id,0,0]

N_cld_2xCO2 = abrupt2xCO2_prp.variables['dR_cld_trad0'][id,0,0]  + abrupt2xCO2_prp.variables['dR_cld_srad0'][id,0,0]
N_cld_4xCO2 = abrupt4xCO2_prp.variables['dR_cld_trad0'][id,0,0] + abrupt4xCO2_prp.variables['dR_cld_srad0'][id,0,0]
N_cld_8xCO2 = abrupt8xCO2_prp.variables['dR_cld_trad0'][id,0,0] + abrupt8xCO2_prp.variables['dR_cld_srad0'][id,0,0]
N_cld_16xCO2 = abrupt16xCO2_prp.variables['dR_cld_trad0'][id,0,0] + abrupt16xCO2_prp.variables['dR_cld_srad0'][id,0,0]

N_alb_2xCO2 = abrupt2xCO2_prp.variables['dR_alb_trad0'][id,0,0]  + abrupt2xCO2_prp.variables['dR_alb_srad0'][id,0,0]
N_alb_4xCO2 = abrupt4xCO2_prp.variables['dR_alb_trad0'][id,0,0] + abrupt4xCO2_prp.variables['dR_alb_srad0'][id,0,0]
N_alb_8xCO2 = abrupt8xCO2_prp.variables['dR_alb_trad0'][id,0,0] + abrupt8xCO2_prp.variables['dR_alb_srad0'][id,0,0]
N_alb_16xCO2 = abrupt16xCO2_prp.variables['dR_alb_trad0'][id,0,0] + abrupt16xCO2_prp.variables['dR_alb_srad0'][id,0,0]


#axes.scatter(T_2xCO2, N_tmp_2xCO2,color='green')
#axes.scatter(T_4xCO2, N_tmp_4xCO2,color='green')
#axes.scatter(T_8xCO2, N_tmp_8xCO2,color='green')
axes.scatter(T_16xCO2, N_tmp_16xCO2,color='green')

#axes.scatter(T_2xCO2, N_vap_2xCO2,color='chocolate')
#axes.scatter(T_4xCO2, N_vap_4xCO2,color='chocolate')
#axes.scatter(T_8xCO2, N_vap_8xCO2,color='chocolate')
axes.scatter(T_16xCO2, N_vap_16xCO2,color='chocolate')

#axes.scatter(T_2xCO2, N_cld_2xCO2,color='blue')
#axes.scatter(T_4xCO2, N_cld_4xCO2,color='blue')
#axes.scatter(T_8xCO2, N_cld_8xCO2,color='blue')
axes.scatter(T_16xCO2, N_cld_16xCO2,color='blue')

#axes.scatter(T_2xCO2, N_alb_2xCO2,color='orange')
#axes.scatter(T_4xCO2, N_alb_4xCO2,color='orange')
#axes.scatter(T_8xCO2, N_alb_8xCO2,color='orange')
axes.scatter(T_16xCO2, N_alb_16xCO2,color='orange')



plt.xlim(0,15)
plt.ylim(-60,60)

fit_T_2xCO2   = np.polyfit(T_2xCO2[idfit],N_tmp_2xCO2[idfit],1)
fit_T_4xCO2   = np.polyfit(T_4xCO2[idfit],N_tmp_4xCO2[idfit],1)
fit_T_8xCO2   = np.polyfit(T_8xCO2[idfit],N_tmp_8xCO2[idfit],1)
fit_T_16xCO2   = np.polyfit(T_16xCO2[idfit],N_tmp_16xCO2[idfit],1)

axes.set_xlabel(r'Global mean temperature change, $T_s$ (K)')
axes.set_ylabel(r'Top-of-atmosphere imbalance, $N$ (Wm$^{-2}$)')

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
plt.savefig('../plots/MPI-ESM12_abrupt_forcing.pdf', dpi=300)
plt.close()

#-------------------------------------------------------------

#N_T   = prp.variables['dR_tmp_trad0'][id,0,0] + prp.variables['dR_tmp_srad0'][id,0,0]
#N_W   = prp.variables['dR_vap_trad0'][id,0,0] + prp.variables['dR_vap_srad0'][id,0,0]
#N_C   = prp.variables['dR_cld_trad0'][id,0,0] + prp.variables['dR_cld_srad0'][id,0,0]
#N_A   = prp.variables['dR_alb_trad0'][id,0,0] + prp.variables['dR_alb_srad0'][id,0,0]
#N_CO2 = prp.variables['dR_co2_trad0'][id,0,0] + prp.variables['dR_co2_srad0'][id,0,0]
#N_PLK = prp.variables['dR_plk_trad0'][id,0,0] + prp.variables['dR_plk_srad0'][id,0,0]
#N_LR  = prp.variables['dR_lr_trad0'][id,0,0]  + prp.variables['dR_lr_srad0'][id,0,0]
#N_STR = prp.variables['dR_str_trad0'][id,0,0] + prp.variables['dR_str_srad0'][id,0,0]

#fit_T   = np.polyfit(T,N_T,1)
#fit_W   = np.polyfit(T,N_W,1)
#fit_C   = np.polyfit(T,N_C,1)
#fit_A   = np.polyfit(T,N_A,1)
#fit_CO2 = np.polyfit(T,N_CO2,1)
#fit10_C = np.polyfit(T[0:9],N_C[0:9],1)
#fit_PLK = np.polyfit(T,N_PLK,1)
#fit_LR  = np.polyfit(T,N_LR,1)
#fit_STR = np.polyfit(T,N_STR,1)
#fit10_LR= np.polyfit(T[0:9],N_LR[0:9],1)

#-------------------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,8))  

xmax=6

axes.scatter(T,N, color='black')
axes.scatter(T,N_T, color='green')
axes.scatter(T,N_W, color='chocolate')
axes.scatter(T,N_C, color='blue')
axes.scatter(T,N_A, color='orange')
#axes.scatter(T,N_CO2, color='purple')

axes.plot((0,-fit[1]/fit[0]),(fit[1],0),color='gray')
axes.plot((0,xmax),(fit_T[1],fit_T[1]+fit_T[0]*xmax),color='green')
axes.plot((0,xmax),(fit_W[1],fit_W[1]+fit_W[0]*xmax),color='chocolate')
axes.plot((0,xmax),(fit_C[1],fit_C[1]+fit_C[0]*xmax),color='blue')
axes.plot((0,xmax),(fit_A[1],fit_A[1]+fit_A[0]*xmax),color='orange')
axes.plot((0,xmax),(fit_CO2[1],fit_CO2[1]+fit_CO2[0]*xmax),color='purple')

axes.plot((0,T10),(fit10[1],fit10[1]+T10*fit10[0]),color='gray',linestyle='--')
axes.plot((0,T10),(fit10_C[1],fit10_C[1]+T10*fit10_C[0]),color='blue',linestyle='--')

axes.scatter(T_hansen,N_hansen, color='red')
axes.plot((0,T_hansen),(N_hansen-fit[0]*T_hansen,N_hansen),color='red')


plt.xlim(0,xmax)
plt.ylim(-20,15)

axes.text(0.2,-10,r'Total imbalance, $N$',color='black')
axes.text(0.2,-11.5,r'Temperature, $N_T$',color='green')
axes.text(0.2,-13,r'Water vapor, $N_W$',color='chocolate')
axes.text(0.2,-14.5,r'Clouds, $N_C$',color='blue')
axes.text(0.2,-16,r'Surface albedo, $N_A$',color='orange')
axes.text(0.2,-17.5,r'Direct CO$_2$ effect, $F^*_{4xCO_2}$',color='purple')

axes.set_xlabel(r'Global mean temperature change, $T_s$ (K)')
axes.set_ylabel(r'Top-of-atmosphere imbalance components, $N_i$ (Wm$^{-2}$)')

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
#plt.savefig('../plots/MPI-ESM12_abrupt4xCO2_PRP.pdf', dpi=300)
plt.close()

#-------------------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,8))  

xmax=6

#axes.scatter(T,N, color='black')
#axes.scatter(T,N_T, color='green')
axes.scatter(T,N_PLK, color='gray')
axes.scatter(T,N_LR, color='steelblue')
axes.scatter(T,N_STR, color='lightgray')

#axes.scatter(T,N_W, color='chocolate')
#axes.scatter(T,N_C, color='blue')
#axes.scatter(T,N_A, color='orange')
#axes.scatter(T,N_CO2, color='purple')

axes.plot((0,xmax),(fit_PLK[1],fit_PLK[1]+fit_PLK[0]*xmax),color='black')
axes.plot((0,xmax),(fit_LR[1],fit_LR[1]+fit_LR[0]*xmax),color='black')
axes.plot((0,T10),(fit10_LR[1],fit10_LR[1]+T10*fit10_LR[0]),color='black',linestyle='--')
axes.plot((0,xmax),(fit_STR[1],fit_STR[1]+fit_STR[0]*xmax),color='black')

plt.xlim(0,xmax)
plt.ylim(-20,15)

axes.text(xmax,4,r'Stratospheric temperature',color='black',ha='right')
axes.text(xmax-0.5,-6,r'Lapse-rate',color='black',ha='right')
axes.text(4.5,-16,r'Planck',color='black',ha='right')

axes.set_xlabel(r'Global mean temperature change, $T_s$ (K)')
axes.set_ylabel(r'Top-of-atmosphere imbalance components, $N_i$ (Wm$^{-2}$)')

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
#plt.savefig('../plots/MPI-ESM12_abrupt4xCO2_PRP_temperatures.pdf', dpi=300)
plt.close()

#-------------------------------------------------------------
