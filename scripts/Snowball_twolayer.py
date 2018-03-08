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

def twolayermodel(forcing_years, input_forcing, gamma=0.8, S=1.0, T_ml0=0.0, T_deep0=0.0, epsilon=0.61):

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

    b          = 0.0
    efficacy   = 1.
    So         = 1360.*S
    sigma      = 5.67e-8
    

    ai = 0.7
    ao = 0.29
    ti = 250.
    to = 288.
    

    # Initialize state variables:
    T_ml       = np.zeros(nstep)
    T_ml[0]    = T_ml0
    T_deep     = np.zeros(nstep)
    T_deep[0]  = T_deep0
    imbalance  = np.zeros(nstep)
    
    # Integrate:
    imbalance[0] = np.nan 
    for t in range(0, nstep-1):
        
        if T_ml[t] <= ti:
            albedo = ai
        elif T_ml[t] >= to:
            albedo = ao
        else:
            albedo = ao + (ai-ao)*(T_ml[t]-to)**2./(ti-to)**2.
            
        N              = So/4.0*(1-albedo) - epsilon*sigma*T_ml[t]**4.0
        T_ml[t+1]      = T_ml[t] + (forcing[t]+N-efficacy*gamma*(T_ml[t]-T_deep[t]))*delta_time/C_ml
        T_deep[t+1]    = T_deep[t] + gamma*(T_ml[t]-T_deep[t])*delta_time/C_deep
        imbalance[t+1] = (C_ml*(T_ml[t+1]-T_ml[t]) + C_deep*(T_deep[t+1]-T_deep[t]))/delta_time
    
    # Output result and settings:
    result['forcing'] = forcing
    result['T_ml'] = T_ml
    result['T_deep'] = T_deep
    result['imbalance'] = imbalance
    result['time'] = timeyears


    result['C_ml'] = C_ml
    result['C_deep'] = C_deep
    result['gamma'] = gamma
    result['efficacy'] = efficacy
        
    return result

#--------------------------------------------------------------------------

nyears     = 3000
forcing_years = np.arange(0,nyears)
input_forcing    = 0*np.ones(nyears)
random_forcing   = np.random.normal(0,1,nyears)

exp1 = twolayermodel(forcing_years, input_forcing, T_ml0=289, T_deep0=289)
exp2 = twolayermodel(forcing_years, input_forcing, S=0.88, T_ml0=289, T_deep0=289)
exp3 = twolayermodel(forcing_years, input_forcing, S=0.87, T_ml0=289, T_deep0=289)
#exp4 = twolayermodel(forcing_years, random_forcing, S=0.875, T_ml0=289, T_deep0=289)

exp5 = twolayermodel(forcing_years, input_forcing, S=1.0, epsilon=0.46, T_ml0=233, T_deep0=233)
exp6 = twolayermodel(forcing_years, input_forcing, S=1.0, T_ml0=233, T_deep0=233)

#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(6,4))

axes.plot(exp1['time'],exp1['T_ml'],color='black')
axes.plot(exp2['time'],exp2['T_ml'],color='purple')
axes.plot(exp3['time'],exp3['T_ml'],color='blue')
#axes.plot(exp4['time'],exp4['T_ml'],color='pink')
axes.plot(exp5['time'],exp5['T_ml'],color='orange')
axes.plot(exp6['time'],exp6['T_ml'],color='red')

axes.set_xlabel('Time (years)')
axes.set_ylabel('Temperature (K)')

plt.xlim((0,nyears))

tyear = 2000

axes.text(tyear,265,r'$S=1.00\cdot S_o$', color='black')
axes.text(tyear,260,r'$S=1.00\cdot S_o$', color='red')
axes.text(tyear,250,r'$S=0.87\cdot S_o$', color='blue')
axes.text(tyear,255,r'$S=0.88\cdot S_o$', color='purple')
#axes.text(tyear,265,r'$S=0.87\cdot S_o$ + variability', color='pink')
axes.text(tyear,245,r'$S=1.00\cdot S_o$, $\epsilon$=0.46', color='orange')


axes.xaxis.set_ticks_position('bottom')
axes.yaxis.set_ticks_position('left')

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
plt.savefig('../plots/Two_layer_snowball.pdf', dpi=300)
plt.close()




