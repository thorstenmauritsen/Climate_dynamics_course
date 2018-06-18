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
import scipy.signal as signal
import math
from scipy.io import netcdf

# Define time parameters in seconds:
hour       = 3600
day        = 24*hour
month      = 30*day
year       = 12*month

# Forcing from a doubling of CO2:
f2x        = 3.7

# set general text options for plotting
plt.rc('text'           , usetex=True)
plt.rc('font'           , size=12)
plt.rc('legend'         , fontsize=12)
plt.rc('text.latex'     , preamble=r'\usepackage{cmbright}')

almost_black            = '#262626'


def yearmean(inyears,x):

    y1         = min(np.floor(inyears))
    y2         = max(np.floor(inyears))
    outyears   = np.arange(y1,y2+1,1)
    nyears     = np.size(outyears)
    y          = np.zeros(nyears)
    i=0
    for year in outyears:
        y[i]   = np.average(x[np.where(np.floor(inyears)==year)])
        i      = i+1

    return y


def twolayermodel(forcing_years, input_forcing, ECS=3.0, gamma=0.8, T_ml0=0.0, T_deep0=0.0, b=0.0, efficacy=1.0, sigma=10):

    result = {}

    nyears     = np.size(forcing_years)
    delta_time = month
    nstep      = round(nyears*year/delta_time)
    time       = delta_time*np.arange(0,nstep,dtype=float) + min(forcing_years)*year
    timeyears  = time/year
    forcing    = np.interp(timeyears,forcing_years,input_forcing)
    noise      = np.random.normal(0.0, sigma/(year/delta_time)**0.5, nstep)

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
        T_ml[t+1]      = T_ml[t] + (forcing[t] + noise[t] +(lambda_0+b*T_ml[t])*T_ml[t]-efficacy*gamma*(T_ml[t]-T_deep[t]))*delta_time/C_ml
        T_deep[t+1]    = T_deep[t] + gamma*(T_ml[t]-T_deep[t])*delta_time/C_deep
        imbalance[t+1] = (C_ml*(T_ml[t+1]-T_ml[t]) + C_deep*(T_deep[t+1]-T_deep[t]))/delta_time

    # Output result and settings for this run:
    result['time'] = timeyears
    result['forcing'] = forcing
    result['T_ml'] = T_ml
    result['T_deep'] = T_deep
    result['imbalance'] = imbalance

#    result['time'] = forcing_years
#    result['forcing'] = input_forcing
#    result['T_ml'] = yearmean(timeyears,T_ml)
#    result['T_deep'] = yearmean(timeyears,T_deep)
#    result['imbalance'] = yearmean(timeyears,imbalance)

    result['ECS'] = ECS
    result['C_ml'] = C_ml
    result['C_deep'] = C_deep
    result['gamma'] = gamma
    result['lambda_0'] = lambda_0
    result['b'] = b
    result['efficacy'] = efficacy
    result['sigma'] = sigma

    result['ar1t'] = np.min(np.corrcoef(T_ml[:-1],T_ml[1:]))**12

    return result

#--------------------------------------------------------------------------

nyears        = 100000
forcing_years = np.arange(0,nyears)
input_forcing = 0*np.linspace(0.0,1.0,nyears)

exp1 = twolayermodel(forcing_years, input_forcing)
exp2 = twolayermodel(forcing_years, input_forcing+0.3)
exp3 = twolayermodel(forcing_years, input_forcing+1)

#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(7,3.5))

color1                  = 'black'
color2                  = "royalblue"
color3                  = 'darkkhaki'
color4                  = 'salmon'
color5                  = 'orange'

lw=1.0

axes.plot(exp1['time'],exp1['T_ml'],lw=lw,color=color1)
axes.plot(exp2['time'],exp2['T_ml'],lw=lw,color=color2)
axes.plot(exp3['time'],exp3['T_ml'],lw=lw,color=color3)


#axes.set_xlabel('Time (years)')
axes.set_xlabel('Time (y)')
axes.set_ylabel(r'Temperature (K)')

xmax = max(axes.get_xlim())
plt.xlim((0,xmax))
ymax = max(axes.get_ylim())


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
plt.savefig('../plots/variability_shift_timeseries.pdf', dpi=300)
plt.close()


#--------------------------------------------------------------------------
# Plot

fig, ax = plt.subplots(1,2, figsize=(7,3.5))

axes = ax[0]

xbins = np.linspace(-0.7,2,100)

ns = 10000

axes.hist(exp1['T_ml'][ns:], bins=xbins, histtype='stepfilled', color=color1)
axes.hist(exp2['T_ml'][ns:], bins=xbins, histtype='stepfilled', color=color2, alpha=0.8)
axes.hist(exp3['T_ml'][ns:], bins=xbins, histtype='stepfilled', color=color3, alpha=0.8)



ymax = max(axes.get_ylim())

axes.plot((0.5, 0.5),(0, ymax/3), linestyle='--',color='black')

axes.set_xlabel('Temperature (K)')


axes.xaxis.set_ticks_position('bottom')
axes.yaxis.set_ticks_position('none')
#axes.yaxes.set_yticklabels()
axes.get_yaxis().set_visible(False)

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['top', 'right', 'left']
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom']
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)


axes = ax[1]


x=np.linspace(0,1.2,1000)
std=np.std(exp1['T_ml'][ns:])

axes.plot(x,100*stats.norm.sf(0.5, loc=x, scale=std),color='black')

axes.scatter(np.mean(exp1['T_ml'][ns:]), 100 - stats.percentileofscore(exp1['T_ml'][ns:],0.5), zorder=10, color=color1)
axes.scatter(np.mean(exp2['T_ml'][ns:]), 100 - stats.percentileofscore(exp2['T_ml'][ns:],0.5), zorder=10, color=color2)
axes.scatter(np.mean(exp3['T_ml'][ns:]), 100 - stats.percentileofscore(exp3['T_ml'][ns:],0.5), zorder=10, color=color3)


axes.set_ylabel('Probability of exceeding 0.5 K')
axes.set_xlabel('Mean temperature (K)')

ymax = max(axes.get_ylim())
plt.ylim((0,ymax))

axes.xaxis.set_ticks_position('bottom')
axes.yaxis.set_ticks_position('right')
axes.yaxis.set_label_position("right")

for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['top', 'left']
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

spines_to_keep = [ 'bottom', 'right']
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)


plt.tight_layout()
plt.savefig('../plots/variability_shift_distribution.pdf', dpi=300)
plt.close()
