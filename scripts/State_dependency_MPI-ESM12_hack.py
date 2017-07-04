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
#abrupt16xCO2_prp = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt16xCO2_echam6_prp_mean_anomalies.nc')
#abrupt16xCO2_bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt16xCO2_echam6_BOT_mean_anomalies.nc')

abrupt16xCO2 = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt16xCO2_echam6_Chris_anomalies.nc')

#amip4xCO2 = netcdf.netcdf_file('../Data/BOT_echam-6.3.02p4_amip4xCO2_1979-2008_timmean_fldmean_anomaly.nc')

#-------------------------------------------------------------

id=np.arange(0,299)
idfit=np.arange(99,299)

T_2xCO2 = abrupt2xCO2_bot.variables['tsurf'][id,0,0]
T_4xCO2 = abrupt4xCO2_bot.variables['tsurf'][id,0,0]
T_8xCO2 = abrupt8xCO2_bot.variables['tsurf'][id,0,0]
T_16xCO2 = abrupt16xCO2.variables['tsurf'][id,0,0]

N_tmp_2xCO2 = abrupt2xCO2_prp.variables['dR_tmp_trad0'][id,0,0]  + abrupt2xCO2_prp.variables['dR_tmp_srad0'][id,0,0]
N_tmp_4xCO2 = abrupt4xCO2_prp.variables['dR_tmp_trad0'][id,0,0] + abrupt4xCO2_prp.variables['dR_tmp_srad0'][id,0,0]
N_tmp_8xCO2 = abrupt8xCO2_prp.variables['dR_tmp_trad0'][id,0,0] + abrupt8xCO2_prp.variables['dR_tmp_srad0'][id,0,0]
N_tmp_16xCO2 = abrupt16xCO2.variables['dR_tmp_trad0'][id,0,0] + abrupt16xCO2.variables['dR_tmp_srad0'][id,0,0]

N_vap_2xCO2 = abrupt2xCO2_prp.variables['dR_vap_trad0'][id,0,0]  + abrupt2xCO2_prp.variables['dR_vap_srad0'][id,0,0]
N_vap_4xCO2 = abrupt4xCO2_prp.variables['dR_vap_trad0'][id,0,0] + abrupt4xCO2_prp.variables['dR_vap_srad0'][id,0,0]
N_vap_8xCO2 = abrupt8xCO2_prp.variables['dR_vap_trad0'][id,0,0] + abrupt8xCO2_prp.variables['dR_vap_srad0'][id,0,0]
N_vap_16xCO2 = abrupt16xCO2.variables['dR_vap_trad0'][id,0,0] + abrupt16xCO2.variables['dR_vap_srad0'][id,0,0]

N_cld_2xCO2 = abrupt2xCO2_prp.variables['dR_cld_trad0'][id,0,0]  + abrupt2xCO2_prp.variables['dR_cld_srad0'][id,0,0]
N_cld_4xCO2 = abrupt4xCO2_prp.variables['dR_cld_trad0'][id,0,0] + abrupt4xCO2_prp.variables['dR_cld_srad0'][id,0,0]
N_cld_8xCO2 = abrupt8xCO2_prp.variables['dR_cld_trad0'][id,0,0] + abrupt8xCO2_prp.variables['dR_cld_srad0'][id,0,0]
N_cld_16xCO2 = abrupt16xCO2.variables['dR_cld_trad0'][id,0,0] + abrupt16xCO2.variables['dR_cld_srad0'][id,0,0]

N_alb_2xCO2 = abrupt2xCO2_prp.variables['dR_alb_trad0'][id,0,0]  + abrupt2xCO2_prp.variables['dR_alb_srad0'][id,0,0]
N_alb_4xCO2 = abrupt4xCO2_prp.variables['dR_alb_trad0'][id,0,0] + abrupt4xCO2_prp.variables['dR_alb_srad0'][id,0,0]
N_alb_8xCO2 = abrupt8xCO2_prp.variables['dR_alb_trad0'][id,0,0] + abrupt8xCO2_prp.variables['dR_alb_srad0'][id,0,0]
N_alb_16xCO2 = abrupt16xCO2.variables['dR_alb_trad0'][id,0,0] + abrupt16xCO2.variables['dR_alb_srad0'][id,0,0]

fit_tmp_2xCO2   = np.polyfit(T_2xCO2[idfit],N_tmp_2xCO2[idfit],1)
fit_tmp_4xCO2   = np.polyfit(T_4xCO2[idfit],N_tmp_4xCO2[idfit],1)
fit_tmp_8xCO2   = np.polyfit(T_8xCO2[idfit],N_tmp_8xCO2[idfit],1)
fit_tmp_16xCO2   = np.polyfit(T_16xCO2[idfit],N_tmp_16xCO2[idfit],1)

fit_vap_2xCO2   = np.polyfit(T_2xCO2[idfit],N_vap_2xCO2[idfit],1)
fit_vap_4xCO2   = np.polyfit(T_4xCO2[idfit],N_vap_4xCO2[idfit],1)
fit_vap_8xCO2   = np.polyfit(T_8xCO2[idfit],N_vap_8xCO2[idfit],1)
fit_vap_16xCO2   = np.polyfit(T_16xCO2[idfit],N_vap_16xCO2[idfit],1)

fit_cld_2xCO2   = np.polyfit(T_2xCO2[idfit],N_cld_2xCO2[idfit],1)
fit_cld_4xCO2   = np.polyfit(T_4xCO2[idfit],N_cld_4xCO2[idfit],1)
fit_cld_8xCO2   = np.polyfit(T_8xCO2[idfit],N_cld_8xCO2[idfit],1)
fit_cld_16xCO2   = np.polyfit(T_16xCO2[idfit],N_cld_16xCO2[idfit],1)

fit_alb_2xCO2   = np.polyfit(T_2xCO2[idfit],N_alb_2xCO2[idfit],1)
fit_alb_4xCO2   = np.polyfit(T_4xCO2[idfit],N_alb_4xCO2[idfit],1)
fit_alb_8xCO2   = np.polyfit(T_8xCO2[idfit],N_alb_8xCO2[idfit],1)
fit_alb_16xCO2   = np.polyfit(T_16xCO2[idfit],N_alb_16xCO2[idfit],1)


# ---------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(5,7))  



#N = bot.variables['trad0'][id,0,0]+bot.variables['srad0'][id,0,0]
#T = bot.variables['tsurf'][id,0,0]





axes.scatter(T_2xCO2, N_tmp_2xCO2,color='green')
#axes.scatter(T_4xCO2, N_tmp_4xCO2,color='green')
#axes.scatter(T_8xCO2, N_tmp_8xCO2,color='green')
axes.scatter(T_16xCO2, N_tmp_16xCO2,color='green',marker='x')

axes.scatter(T_2xCO2, N_vap_2xCO2,color='chocolate')
#axes.scatter(T_4xCO2, N_vap_4xCO2,color='chocolate')
#axes.scatter(T_8xCO2, N_vap_8xCO2,color='chocolate')
axes.scatter(T_16xCO2, N_vap_16xCO2,color='chocolate',marker='x')

axes.scatter(T_2xCO2, N_cld_2xCO2,color='blue')
#axes.scatter(T_4xCO2, N_cld_4xCO2,color='blue')
#axes.scatter(T_8xCO2, N_cld_8xCO2,color='blue')
axes.scatter(T_16xCO2, N_cld_16xCO2,color='blue',marker='x')

axes.scatter(T_2xCO2, N_alb_2xCO2,color='orange')
#axes.scatter(T_4xCO2, N_alb_4xCO2,color='orange')
#axes.scatter(T_8xCO2, N_alb_8xCO2,color='orange')
axes.scatter(T_16xCO2, N_alb_16xCO2,color='orange',marker='x')



plt.xlim(0,15)
plt.ylim(-60,40)





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

