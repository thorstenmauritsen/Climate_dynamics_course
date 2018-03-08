import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np                
import scipy.stats as stat
import math
from scipy.io import netcdf
#from netCDF4 import Dataset

# set general text options for plotting
plt.rc('text'           , usetex=True)
plt.rc('font'           , size=12)
plt.rc('legend'         , fontsize=12)
plt.rc('text.latex'     , preamble=r'\usepackage{cmbright}')

almost_black            = '#262626'

# Define time parameters in seconds:
hour       = 3600
day        = 24*hour
month      = 30*day
year       = 12*month

# Forcing from a doubling of CO2:
f2x        = 3.7

#--------------------------------------------------------------------------

def twolayermodel(forcing_years, input_forcing, ECS=3.0, gamma=0.8, T_ml0=0.0):

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
    C_ml       = 50*density*c_w
    C_deep     = 1200*density*c_w
    lambda_0   = -f2x/ECS
    b          = 0.0
    efficacy   = 1.

    # Initialize state variables:
    T_ml       = np.zeros(nstep)
    T_ml[0]    = T_ml0
    T_deep     = np.zeros(nstep)
    imbalance  = np.zeros(nstep)
    
    # Integrate:
    imbalance[0] = np.nan 
    for t in range(0, nstep-1):
        T_ml[t+1]      = T_ml[t] + (forcing[t]+(lambda_0+b*T_ml[t])*T_ml[t]-efficacy*gamma*(T_ml[t]-T_deep[t]))*delta_time/C_ml
        T_deep[t+1]    = T_deep[t] + gamma*(T_ml[t]-T_deep[t])*delta_time/C_deep
        imbalance[t+1] = (C_ml*(T_ml[t+1]-T_ml[t]) + C_deep*(T_deep[t+1]-T_deep[t]))/delta_time
    
    # Output result and settings:
    result['forcing'] = forcing
    result['T_ml'] = T_ml
    result['T_deep'] = T_deep
    result['imbalance'] = imbalance
    result['time'] = timeyears

    result['ECS'] = ECS
    result['C_ml'] = C_ml
    result['C_deep'] = C_deep
    result['gamma'] = gamma
    result['lambda_0'] = lambda_0
    result['b'] = b
    result['efficacy'] = efficacy
        
    return result

#--------------------------------------------------------------------------

nyears     = 100
forcing_years = np.arange(0,nyears)
input_forcing    = f2x*np.ones(nyears)

exp1 = twolayermodel(forcing_years, input_forcing, ECS=2.77 )
exp2 = twolayermodel(forcing_years, input_forcing, gamma=0, ECS=2.77)

mpiesm = netcdf.netcdf_file('../Data/mpiesm-1.2.00_abrupt2xCO2_anomaly.nc', 'r')

#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(8,4))

axes.plot([0.0,10.0], [0.0,10.0*f2x/exp1['C_ml']*year], color='gray', linestyle='--')

axes.plot(exp1['time'],exp1['T_ml'],color='blue')
axes.plot(exp2['time'],exp2['T_ml'],color='black')

axes.scatter(np.arange(0,nyears), mpiesm.variables['tsurf'][0:nyears,0,0], color='lightblue')

axes.scatter(-exp2['C_ml']/exp2['lambda_0']/year,exp2['ECS']*(1-1/math.e), color='k')


axes.set_xlabel('Time (years)')
axes.set_ylabel('Temperature (K)')

plt.xlim((0,nyears))
plt.ylim((0,4.0))

axes.xaxis.set_ticks_position('bottom')

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
plt.savefig('../plots/Two_layer_temperature_abrupt_forcing.pdf', dpi=300)
plt.close()

mpiesm.close


