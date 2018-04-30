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


def power_spectrum(exp):
    frequency  = np.linspace(0,12./2.,int(nyears/2))
    nbins= 40
    bins = np.exp(np.linspace(-7,2,nbins))
    P_ml = np.fft.fft(signal.detrend(exp['T_ml']))
    P_ml = P_ml[0:int(nyears/2)]
    data = abs(P_ml)**2
    XX = np.array([np.mean((bins[binInd], bins[binInd+1])) for binInd in range(nbins-1)])
    YY = np.array([np.mean(data[(frequency > bins[binInd]) & (frequency <= bins[binInd+1])]) for binInd in range(nbins-1)])
    return XX, YY


def twolayermodel(forcing_years, input_forcing, ECS=3.0, gamma=0.8, T_ml0=0.0, T_deep0=0.0, b=0.0, efficacy=1.0, sigma=0.1):

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
        T_ml[t+1]      = T_ml[t] + (forcing[t]+(lambda_0+b*T_ml[t])*T_ml[t]-efficacy*gamma*(T_ml[t]-T_deep[t]))*delta_time/C_ml + noise[t]
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

exp1 = twolayermodel(forcing_years, input_forcing,   ECS=1)
exp2 = twolayermodel(forcing_years, input_forcing,   ECS=3)
exp3 = twolayermodel(forcing_years, input_forcing,   ECS=9)
exp4 = twolayermodel(forcing_years, input_forcing,   ECS=27)
exp5 = twolayermodel(forcing_years, input_forcing,   ECS=1e20)

mlo1 = twolayermodel(forcing_years, input_forcing, gamma=0,   ECS=1)
mlo2 = twolayermodel(forcing_years, input_forcing, gamma=0,   ECS=3)
mlo3 = twolayermodel(forcing_years, input_forcing, gamma=0,   ECS=9)
mlo4 = twolayermodel(forcing_years, input_forcing, gamma=0,   ECS=27)
mlo5 = twolayermodel(forcing_years, input_forcing, gamma=0,   ECS=1e20)

#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(7,3.5))

color1                  = 'black'
color2                  = "royalblue"
color3                  = 'darkkhaki'
color4                  = 'salmon'
color5                  = 'orange'

lw=1.0

#axes.plot(exp4['time'],exp4['T_ml'],lw=lw,color=color4)

axes.plot(exp1['time'],exp1['T_ml'],lw=lw,color=color1)
axes.plot(exp2['time'],exp2['T_ml'],lw=lw,color=color2)
axes.plot(exp3['time'],exp3['T_ml'],lw=lw,color=color3)
axes.plot(exp4['time'],exp4['T_ml'],lw=lw,color=color4)
axes.plot(exp5['time'],exp5['T_ml'],lw=lw,color=color5)

axes.plot(mlo1['time'],mlo1['T_ml'],lw=lw,color='green')
axes.plot(mlo2['time'],mlo2['T_ml'],lw=lw,color=color2)
axes.plot(mlo3['time'],mlo3['T_ml'],lw=lw,color=color3)
axes.plot(mlo4['time'],mlo4['T_ml'],lw=lw,color=color4)
axes.plot(mlo5['time'],mlo5['T_ml'],lw=lw,color=color5)


#axes.text(10,3.5,'ECS = '+str(exp1['ECS'])+r' K, $\epsilon$ = '+str(exp1['efficacy']), color=color1)


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
plt.savefig('../plots/Global_variability.pdf', dpi=300)
plt.close()


#------------------------------------------------


fig, axes = plt.subplots(1,2, figsize=(10,5))

x,y = power_spectrum(exp1)
axes[1].loglog(1/x, y, color=color1)

x,y = power_spectrum(exp2)
axes[1].loglog(1/x, y, color=color2)

x,y = power_spectrum(exp3)
axes[1].loglog(1/x, y, color=color3)

x,y = power_spectrum(exp4)
axes[1].loglog(1/x, y, color=color4)

x,y = power_spectrum(exp5)
axes[1].loglog(1/x, y, color=color5)

axes[1].yaxis.set_ticks_position('right')
axes[1].set_xlabel('Time-scale (years)')
axes[1].set_title("Two-layer model")


x,y = power_spectrum(mlo1)
l1, = axes[0].loglog(1/x, y, color=color1, label='ECS = '+str(mlo1['ECS'])+ 'K')

x,y = power_spectrum(mlo2)
l2, = axes[0].loglog(1/x, y, color=color2, label='ECS = '+str(mlo2['ECS'])+ 'K')

x,y = power_spectrum(mlo3)
l3, = axes[0].loglog(1/x, y, color=color3, label='ECS = '+str(mlo3['ECS'])+ 'K')

x,y = power_spectrum(mlo4)
l4, = axes[0].loglog(1/x, y, color=color4, label=r'ECS = '+str(mlo4['ECS'])+ 'K')

x,y = power_spectrum(mlo5)
l5, = axes[0].loglog(1/x, y, color=color5, label=r'ECS $\approx$ $\infty$')

axes[0].set_xlabel('Time-scale (years)')
axes[0].set_ylabel(r'Power spectral density (K$^2$)')
axes[0].set_title("Mixed-layer model")

axes[0].legend(handles=[l1, l2, l3, l4, l5], fontsize=10)


plt.setp(axes, xlim=(5e-2,2e3), ylim=(1e4,1e12))
plt.setp(axes, xticks=(0.1, 1, 10, 100, 1000), xticklabels=('0.1', '1', '10', '100', '1000'))


plt.tight_layout()
plt.savefig('../plots/Global_variability_spectrum.pdf', dpi=300)
plt.close()

#----------------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(6,3.5))


mpiesm11_picontrol = netcdf.netcdf_file('../Data/tsurf-grandens/MPI-ESM1.1_piControl_lkm0001_fldmean.nc')
mpiesm11_ar1t = np.min(np.corrcoef(mpiesm11_picontrol.variables['tsurf'][:-1,0,0],mpiesm11_picontrol.variables['tsurf'][1:,0,0]))
mpim = axes.scatter(mpiesm11_ar1t,2.77, color='black', label='MPI-ESM1.1')

#axes.plot((0.3,0.8),(2.1,4.3),linestyle='--', color='black')
hadcrut, = axes.plot((0.42,0.42),(0.0,6.0),linestyle='-.', color='black', label='Observed')

colormlm = 'red'
colortlm = 'blue'

axes.scatter(mlo1['ar1t'], mlo1['ECS'], color=colormlm)
axes.scatter(mlo2['ar1t'], mlo2['ECS'], color=colormlm)
axes.scatter(mlo3['ar1t'], mlo3['ECS'], color=colormlm)

x=np.arange(0.1,0.95,0.01)
theoryb = -f2x/(np.log(x)*mlo1['C_ml']/year+mlo1['gamma']*mlo1['efficacy'])
theoryb[np.where(theoryb<0)] = np.nan
mlm, = axes.plot(x,theoryb, color=colormlm, label='Mixed-layer model')

axes.scatter(exp1['ar1t'], exp1['ECS'], color=colortlm)
axes.scatter(exp2['ar1t'], exp2['ECS'], color=colortlm)
axes.scatter(exp3['ar1t'], exp3['ECS'], color=colortlm)



theorya = -f2x/(np.log(x)*exp1['C_ml']/year+exp1['gamma']*exp1['efficacy'])
theorya[np.where(theorya<0)] = np.nan
tlm, = axes.plot(x,theorya, color=colortlm, label='Two-layer model')

axes.legend(handles=[mlm,tlm,mpim, hadcrut], fontsize=10)

axes.set_xlabel(r'Year-1 lag correlation ($\alpha_{1T}$)')
axes.set_ylabel('ECS (K)')

plt.setp(axes, xlim=(0,1), ylim=(0,6))

plt.tight_layout()
plt.savefig('../plots/Global_variability_lag-correlation.pdf', dpi=300)
plt.close()
