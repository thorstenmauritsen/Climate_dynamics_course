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

# Define time parameters in seconds:
hour       = 3600
day        = 24*hour
month      = 30*day
year       = 12*month

# Forcing from a doubling of CO2:
f2x        = 4.5

# set general text options for plotting
plt.rc('text'           , usetex=True)
plt.rc('font'           , size=12)
plt.rc('legend'         , fontsize=12)
plt.rc('text.latex'     , preamble=r'\usepackage{cmbright}')

almost_black            = '#262626'

def twolayermodel(forcing_years, input_forcing, ECS=3.0, gamma=0.8, T_ml0=0.0, T_deep0=0.0, b=0.0, efficacy=1.0):

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

nyears        = 500
forcing_years = np.arange(0,nyears)
input_forcing = f2x*np.ones(nyears)
random_forcing= np.random.normal(0,1,nyears)

forcing_years_long = np.arange(0,nyears*2)
input_forcing_long = f2x*np.ones(nyears*2)

b=0.045
ECS=2.6
efficacy=1.0

exp1 = twolayermodel(forcing_years, input_forcing,   ECS=ECS, b=b, efficacy=efficacy)
exp2 = twolayermodel(forcing_years, 2*input_forcing, ECS=ECS, b=b, efficacy=efficacy)
exp3 = twolayermodel(forcing_years, 3*input_forcing, ECS=ECS, b=b, efficacy=efficacy)
exp4 = twolayermodel(forcing_years, 4*input_forcing, ECS=ECS, b=b, efficacy=efficacy)
exp4b= twolayermodel(forcing_years_long, 4*input_forcing_long,   ECS=ECS, b=b, efficacy=efficacy)
#exp5 = twolayermodel(forcing_years, 5*input_forcing, ECS=ECS, b=b, efficacy=efficacy, gamma=0)

#exp5 = twolayermodel(forcing_years, 4*input_forcing+random_forcing, ECS=ECS, b=b)

#exp5 = twolayermodel(forcing_years, input_forcing+random_forcing,   ECS=ECS, b=b)
#exp6 = twolayermodel(forcing_years, 2*input_forcing+random_forcing, ECS=ECS, b=b)
#exp7 = twolayermodel(forcing_years, 3*input_forcing+random_forcing, ECS=ECS, b=b)
#exp8 = twolayermodel(forcing_years, 4*input_forcing+random_forcing, ECS=ECS, b=b)

#expo = twolayermodel(forcing_years, 1.2*input_forcing, ECS=f2x, b=0.05)

#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(5,3.5))

color2                  = 'royalblue'
color4                  = "deepskyblue"
color8                  = 'darkkhaki'
color16                 = 'salmon'
lw=3

axes.plot(exp1['T_ml'],exp1['imbalance'],lw=lw,color=color2)
axes.plot(exp2['T_ml'],exp2['imbalance'],lw=lw,color=color4)
axes.plot(exp3['T_ml'],exp3['imbalance'],lw=lw,color=color8)
axes.plot(exp4b['T_ml'],exp4b['imbalance'],lw=lw/2,color=color16,linestyle='--')
axes.plot(exp4['T_ml'],exp4['imbalance'],lw=lw,color=color16)
#axes.plot(exp5['T_ml'],exp5['imbalance'],lw=lw,color=color16)

#axes.plot(expo['T_ml'],expo['imbalance'],lw=lw,color=color16)

#axes.scatter(exp5['T_ml'],exp5['imbalance'],color=color16)
#axes.scatter(exp6['T_ml'],exp6['imbalance'])
#axes.scatter(exp7['T_ml'],exp7['imbalance'])
#axes.scatter(exp8['T_ml'],exp8['imbalance'])

#axes.set_xlabel('Time (years)')
axes.set_xlabel('Temperature (K)')
axes.set_ylabel(r'Imblance (Wm$^{-2}$)')

plt.xlim((0,21))
#plt.ylim((0,16.5))
plt.ylim((0,0.5+4*f2x))


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
plt.savefig('../plots/Two_layer_like_MPI-ESM12.pdf', dpi=300)
plt.close()




#--------------------------------------------------------------------------

#nyears        = 20000
#forcing_years = np.arange(0,nyears)
#input_forcing = f2x*np.ones(nyears) #np.linspace(0.0,1.0,nyears)
#random_forcing= np.random.normal(0,1,nyears)

#b=0.05
#ECS=f2x
#efficacy=1.0

#exp10 = twolayermodel(forcing_years, 1.11*input_forcing,   ECS=ECS, b=b, efficacy=efficacy)
#exp11 = twolayermodel(forcing_years, 1.11*input_forcing+random_forcing,   ECS=ECS, b=b, efficacy=efficacy)
#exp12 = twolayermodel(forcing_years, 1.12*input_forcing,   ECS=ECS, b=b, efficacy=efficacy)

#--------------------------------------------------------------------------
# Plot

#fig, axes = plt.subplots(1,1, figsize=(5,3.5))#

#axes.plot(exp12['time'],exp12['T_ml'],lw=1,color='blue')
#axes.plot(exp11['time'],exp11['T_ml'],lw=1,color='red')
#axes.plot(exp10['time'],exp10['T_ml'],lw=lw,color='black')

#axes.set_xlabel('Time (years)')
#axes.set_ylabel('Temperature (K)')
#axes.set_ylabel(r'Imblance (Wm$^{-2}$)')

#plt.xlim((0,21))
#plt.ylim((0,30))
#plt.ylim((0,0.5+4*f2x))


#axes.xaxis.set_ticks_position('bottom')
#axes.yaxis.set_ticks_position('left')

#for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
#    ticks.set_color(almost_black)

#spines_to_remove        = ['top', 'right'] 
#for spine in spines_to_remove:
#    axes.spines[spine].set_visible(False)

#spines_to_keep = [ 'bottom', 'left']     
#for spine in spines_to_keep:
#    axes.spines[spine].set_linewidth(0.5)
#    axes.spines[spine].set_color(almost_black)



#plt.tight_layout()
#plt.savefig('../plots/Linear_instability.pdf', dpi=300)
#plt.close()

