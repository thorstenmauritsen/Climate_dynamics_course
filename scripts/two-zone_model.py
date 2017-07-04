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

#----------------------------------------------------------------------

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

#----------------------------------------------------------------------

def twozonemodel(forcing_years, input_forcing, T_1_0=0.0, T_2_0=0.0, coupling=0.0):

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
    C_1        = 75*0.7*density*c_w
    C_2        = 2*75*0.7*density*c_w
    lambda_1   = -3.
    lambda_2   = -1.
   
    
    # Initialize state variables:
    T_1        = np.zeros(nstep)
    T_2        = np.zeros(nstep)
    T_1[0]     = T_1_0
    T_2[0]     = T_2_0
    imbalance  = np.zeros(nstep)
    w_1        = 0.5
    w_2        = 1.0-w_1

    # Calculate ECS:
    ECS = -w_1*f2x/lambda_1 -w_2*f2x/lambda_2
    
    # Integrate:
    imbalance[0] = np.nan 
    for t in range(0, nstep-1):
        T_1[t+1]       = T_1[t] + (forcing[t]+lambda_1*T_1[t]+coupling*(T_2[t]-T_1[t]))*delta_time/C_1
        T_2[t+1]       = T_2[t] + (forcing[t]+lambda_2*T_2[t]-coupling*(T_2[t]-T_1[t]))*delta_time/C_2
        imbalance[t+1] = (w_1*C_1*(T_1[t+1]-T_1[t]) + w_2*C_2*(T_2[t+1]-T_2[t]))/delta_time
        
    # Output result and settings for this run:
    result['time'] = timeyears
    result['forcing'] = forcing
    result['T']   = w_1*T_1+w_2*T_2
    result['T_1'] = T_1
    result['T_2'] = T_2
    result['imbalance'] = imbalance
    
    result['ECS'] = ECS
    result['C_1'] = C_1
    result['C_2'] = C_2

    return result

#--------------------------------------------------------------------------

nyears        = 500
forcing_years = np.arange(0,nyears)
input_forcing = f2x*np.ones(nyears)
random_forcing= np.random.normal(0,1,nyears)

forcing_years_long = np.arange(0,nyears*2)
input_forcing_long = f2x*np.ones(nyears*2)

coupling = 0.0

exp1 = twozonemodel(forcing_years, 1*input_forcing, coupling=coupling)
exp2 = twozonemodel(forcing_years, 2*input_forcing, coupling=coupling)
exp3 = twozonemodel(forcing_years, 3*input_forcing, coupling=coupling)
exp4 = twozonemodel(forcing_years, 4*input_forcing, coupling=coupling)



exp5 = twolayermodel(forcing_years, 2*input_forcing, ECS=3.0, efficacy=1.5)
exp6 = twolayermodel(forcing_years, 2*input_forcing, ECS=2.6, b=0.045)
exp7 = twozonemodel(forcing_years, 2*input_forcing, coupling=3.9)
exp0 = twolayermodel(forcing_years, 2*input_forcing, ECS=3.0)

#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(5,3.5))

color2                  = 'royalblue'
color4                  = "deepskyblue"
color8                  = 'darkkhaki'
color16                 = 'salmon'
lw=3

axes.plot(exp1['T'],exp1['imbalance'],lw=lw,color=color2)
axes.plot(exp2['T'],exp2['imbalance'],lw=lw,color=color4)
axes.plot(exp3['T'],exp3['imbalance'],lw=lw,color=color8)
axes.plot(exp4['T'],exp4['imbalance'],lw=lw,color=color16)

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
plt.savefig('../plots/Two_zone_model.pdf', dpi=300)
plt.close()



#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(5,4))

lw=2

bot = netcdf.netcdf_file('../Data/mpiesm-1.2.prp_abrupt4xCO2_echam6_BOT_mean_anomalies.nc')

id=np.arange(0,299)

N = bot.variables['trad0'][id,0,0]+bot.variables['srad0'][id,0,0]
T = bot.variables['tsurf'][id,0,0]

fit = np.polyfit(T,N,1)
#fit_fn = np.poly1d(fit) 

# Gregory method:
axes.scatter(T,N, color='lightgray')
axes.plot((0,-fit[1]/fit[0]),(fit[1],0),color='gray')

# 10-year fit:
fit10 = np.polyfit(T[0:9],N[0:9],1)
T10 = T[9]
axes.plot((0,T10),(fit10[1],fit10[1]+T10*fit10[0]),color='gray',linestyle='--')

axes.plot(exp2['T'],exp2['imbalance'],lw=lw,color=color4)
axes.plot(exp7['T'],exp7['imbalance'],lw=lw,linestyle='--', color=color4)
axes.plot(exp5['T_ml'],exp5['imbalance'],lw=lw,color=color16)
axes.plot(exp6['T_ml'],exp6['imbalance'],lw=lw,color='green')

axes.plot(exp0['T_ml'],exp0['imbalance'],lw=2,color='black')


axes.text(6,9.1, r'MPI-ESM1.2, abrupt 4xCO$_2$', color='gray', ha='right')
axes.text(6,8.5, r'Linear model, ECS='+str(exp0['ECS'])+'K', color='black', ha='right')
axes.text(6,7.9, r'Winton-Held model, $\epsilon$='+str(exp5['efficacy']), color=color16, ha='right')
axes.text(6,7.3, r'Two-zone model', color=color4, ha='right')
axes.text(6,6.7, r'Quadratic model, b='+str(exp6['b']), color='green', ha='right') #+r'Wm$^{-2}$K$^{-2}$'

plt.xlim(0,6)
plt.ylim(0,10)

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
plt.savefig('../plots/Time-dependent_models_abrupt4xCO2.pdf', dpi=300)
plt.close()




#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(5,4))

axes.scatter((bot.variables['time'][id]-bot.variables['time'][0])/10000,T, color='lightgray')

axes.plot(exp2['time'],exp2['T'],lw=lw,color=color4)
axes.plot(exp7['time'],exp7['T'],lw=lw,linestyle='--',color=color4)
axes.plot(exp5['time'],exp5['T_ml'],lw=lw,color=color16)
axes.plot(exp6['time'],exp6['T_ml'],lw=lw,color='green')
axes.plot(exp0['time'],exp0['T_ml'],lw=lw,color='black')


axes.text(300,2.5, r'MPI-ESM1.2, abrupt 4xCO$_2$', color='gray', ha='right')
axes.text(300,2.0, r'Linear model, ECS='+str(exp0['ECS'])+'K', color='black', ha='right')
axes.text(300,1.5, r'Winton-Held model, $\epsilon$='+str(exp5['efficacy']), color=color16, ha='right')
axes.text(300,1.0, r'Two-zone model', color=color4, ha='right')
axes.text(300,0.5, r'Quadratic model, b='+str(exp6['b']), color='green', ha='right') #+r'Wm$^{-2}$K$^{-2}$'


plt.xlim(0,300)
plt.ylim(0,7)

axes.set_xlabel(r'Time (y)')
axes.set_ylabel(r'Global mean temperature change, $T_s$ (K)')

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
plt.savefig('../plots/Time-dependent_models_abrupt4xCO2_timeseries.pdf', dpi=300)
plt.close()
