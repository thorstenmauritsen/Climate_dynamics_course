import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np                

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
# Run experiments

nyears     = 100
forcing_years = np.arange(0,nyears)
input_forcing    = f2x*np.ones(nyears)

exp1 = twolayermodel(forcing_years, input_forcing, ECS=1.0 )
exp2 = twolayermodel(forcing_years, input_forcing, ECS=3.0)
exp3 = twolayermodel(forcing_years, input_forcing, ECS=10.0)

#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(5,4))

axes.plot(exp1['time'],exp1['T_ml'],color='blue')
axes.plot(exp2['time'],exp2['T_ml'],color='black')
axes.plot(exp3['time'],exp3['T_ml'],color='red')

axes.set_xlabel('Time (years)')
axes.set_ylabel('Temperature (K)')

plt.xlim(xmin = 0)
plt.ylim(ymin = 0)

plt.tight_layout()
plt.savefig('Two_layer_temperature_abrupt_forcing.pdf', dpi=300)
plt.close()

