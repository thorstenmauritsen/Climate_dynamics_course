#--------------------------------------------------------------------------
# The two-layer model as a function. The two first arguments are mandatory
# in twolayermodel, whereas the rest is optional:
#
# forcing_years: a list of the years in the input_forcing
# input_forcing: the yearly forcing value
# ECS:           the equilibrium climate sensitivity
# gamma:         deep ocean heat uptake efficiency
# T_ml0:         the mixed-layer initial temperature
# T_deep0:       the deep ocean initial temperature
# b:             state-dependency of feedback parameter (typically -0.1 to 0.1)
# efficacy:      ocean heat uptake efficacy parameter (typically 1.0 to 2.0)
#
#--------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math
from scipy.io import netcdf

# set general text options for plotting
plt.rc('text'           , usetex=True)
plt.rc('font'           , size=12)
plt.rc('legend'         , fontsize=12)
plt.rc('text.latex'     , preamble=r'\usepackage{cmbright}')

almost_black = '#262626'

# Define time parameters in seconds:
hour       = 3600
day        = 24*hour
month      = 30*day
year       = 12*month

# Forcing from a doubling of CO2:
f2x        = 3.7

def twolayermodel(forcing_years, input_forcing, ECS=3.0, gamma=1.0, T_ml0=0.0, T_deep0=0.0, b=0.0, efficacy=1.0):

    result = {}

    nyears     = np.size(forcing_years)
    delta_time = month
    nstep      = round(nyears*year/delta_time)
    time       = delta_time*np.arange(0,nstep,dtype=float) + min(forcing_years)*year
    timeyears  = time/year
    forcing    = np.interp(timeyears,forcing_years,input_forcing)

    # Parameters:
    density    = 1000.
    c_w        = 4181.
    C_ml       = 75*0.7*density*c_w
    C_deep     = 1000*0.7*density*c_w
    lambda_0   = -f2x/ECS

    # Initialize state variables:
    T_ml       = np.zeros(nstep)
    T_deep     = np.zeros(nstep)
    T_ml[0]    = T_ml0
    T_deep[0]  = T_deep0
    imbalance  = np.zeros(nstep)
    
    # Integrate:
    imbalance[0] = np.nan 
    for t in range(0, nstep-1):
        T_ml[t+1]      = T_ml[t] + (forcing[t]+(lambda_0+b*T_ml[t])*T_ml[t]-efficacy*gamma*(T_ml[t]-T_deep[t]))*delta_time/C_ml
        T_deep[t+1]    = T_deep[t] + gamma*(T_ml[t]-T_deep[t])*delta_time/C_deep
        imbalance[t+1] = (C_ml*(T_ml[t+1]-T_ml[t]) + C_deep*(T_deep[t+1]-T_deep[t]))/delta_time
    
    # Output result and settings for this run:
    result['time'] = timeyears
    result['forcing'] = forcing
    result['T_ml'] = T_ml
    result['T_deep'] = T_deep
    result['imbalance'] = imbalance

    result['ECS'] = ECS
    result['C_ml'] = C_ml
    result['C_deep'] = C_deep
    result['gamma'] = gamma
    result['lambda_0'] = lambda_0
    result['b'] = b
    result['efficacy'] = efficacy
        
    return result

#--------------------------------------------------------------------------

# Forcing scenario:
F_year, F_co2, F_ghgother, F_o3trop, F_o3strat, F_aero, F_luc, F_h2ostrat, F_bcsnow, F_contrails, F_solar, F_volc = \
  np.loadtxt('../Data/WG1_AR5_Annex_II_sheet_1-2_forcings.txt', skiprows=1, unpack=True)

F_volc = F_volc*1.00 
F_aero = F_aero*1.0
  
F_total = F_co2 + F_ghgother + F_o3trop + F_o3strat + F_aero + F_luc + F_h2ostrat + F_bcsnow + F_contrails + F_solar + F_volc #- np.mean(F_volc)
F_ghg   = F_co2 + F_ghgother + F_h2ostrat + F_o3trop + F_o3strat 
F_aero  = F_aero + F_bcsnow + F_contrails


ECS=1.78

id  = np.where(F_year >= 1850)
id1 = np.where(F_year == 1850)

exp1 = twolayermodel(F_year[id], F_total[id]-F_total[id1], ECS=ECS)
exp2 = twolayermodel(F_year[id], F_ghg[id]-F_ghg[id1], ECS=ECS)
exp3 = twolayermodel(F_year[id], F_aero[id]-F_aero[id1], ECS=ECS)
exp4 = twolayermodel(F_year[id], F_luc[id]-F_luc[id1], ECS=ECS)
exp5 = twolayermodel(F_year[id], F_volc[id], ECS=ECS)
exp6 = twolayermodel(F_year[id], F_solar[id], ECS=ECS)

#--------------------------------------------------------------------------

# Observed temperature:
f_HadCRU    = netcdf.netcdf_file('../Data/HadCRUT.4.5.0.0.median_means.nc')
HadCRU_year = f_HadCRU.variables['time'][:]
HadCRU_t    = f_HadCRU.variables['temperature_anomaly'][:,0,0]

#--------------------------------------------------------------------------

# Plot

fig, axes = plt.subplots(1,1, figsize=(6,4))

axes.plot(HadCRU_year,HadCRU_t-np.mean(HadCRU_t)+np.mean(exp1['T_ml']),color='lightgray')

axes.plot(exp2['time'],exp2['T_ml'],color='green')
axes.plot(exp3['time'],exp3['T_ml'],color='violet')
axes.plot(exp4['time'],exp4['T_ml'],color='blue')
axes.plot(exp5['time'],exp5['T_ml'],color='chocolate')
axes.plot(exp6['time'],exp6['T_ml'],color='orange')
axes.plot(exp1['time'],exp1['T_ml'],color='black',linewidth=1.5)

plt.xlim((1850,2020))
plt.ylim((-0.6,1.2))

axes.set_xlabel('Time (years)')
axes.set_ylabel('Temperature (K)')

axes.text(1860,1.1,'Observed, HadCRUT',color='gray')
axes.text(1860,1.0,'Total forcing',color='black')
axes.text(1860,0.9,'Greenhouse gases, ozone, stratospheric water',color='green')
axes.text(1860,0.8,'Aerosol, BC on snow, contrails',color='violet')
axes.text(1860,0.7,'Land-use change',color='blue')
axes.text(1860,0.6,'Volcanoes',color='chocolate')
axes.text(1860,0.5,'Solar',color='orange')

#axes.xaxis.set_ticks_position('top')
axes.spines['top'].set_position('zero')
axes.yaxis.set_ticks_position('left')

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['right'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom', 'top', 'left']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

plt.tight_layout()
plt.savefig('../plots/Temperature_two_layer_model_historical.pdf', dpi=300)
plt.close()

#--------------------------------------------------------------------------

# Plot

fig, axes = plt.subplots(1,1, figsize=(6,5))

axes.plot(exp2['time'],exp2['forcing'],color='green')
axes.plot(exp3['time'],exp3['forcing'],color='violet')
axes.plot(exp4['time'],exp4['forcing'],color='blue')
axes.plot(exp5['time'],exp5['forcing'],color='chocolate')
axes.plot(exp6['time'],exp6['forcing'],color='orange')
axes.plot(exp1['time'],exp1['forcing'],color='black',linewidth=1.5)

plt.xlim((1850,2020))
plt.ylim((-4.0,4.0))

axes.set_xlabel('Time (years)')
axes.set_ylabel(r'Forcing relative to 1850 (Wm$^{-2}$)')

axes.text(1860,2.7,'Total forcing',color='black')
axes.text(1860,2.4,'Greenhouse gases, ozone, stratospheric water',color='green')
axes.text(1860,2.1,'Aerosol, BC on snow, contrails',color='violet')
axes.text(1860,1.8,'Land-use change',color='blue')
axes.text(1860,1.5,'Volcanoes',color='chocolate')
axes.text(1860,1.2,'Solar',color='orange')

#axes.xaxis.set_ticks_position('top')
axes.spines['top'].set_position('zero')
axes.yaxis.set_ticks_position('left')

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['right'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom', 'top', 'left']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

plt.tight_layout()
plt.savefig('../plots/Forcing_two_layer_model_historical.pdf', dpi=300)
plt.close()


