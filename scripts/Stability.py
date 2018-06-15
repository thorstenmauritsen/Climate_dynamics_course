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
    nstep = np.size(exp['T_ml'])
    frequency  = np.linspace(0,12./2.,int(nstep/2))
    nbins= 40
    bins = np.exp(np.linspace(-9,2,nbins))
    P_ml = np.fft.fft(signal.detrend(exp['T_ml']))
    P_ml = P_ml[0:int(nstep/2)]
    data = abs(P_ml)**2/nstep
    XX = np.array([np.mean((bins[binInd], bins[binInd+1])) for binInd in range(nbins-1)])
    YY = np.array([np.mean(data[(frequency > bins[binInd]) & (frequency <= bins[binInd+1])]) for binInd in range(nbins-1)])
    return XX, YY

def qubicon(forcing_years, input_forcing, ECS=3.7, gamma=0.8, T_ml0=0.0, T_deep0=0.0, b=0.0,
c=0.0, efficacy=1.0, sigma=2.42):

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
        T_ml[t+1]      = T_ml[t] + (forcing[t] + noise[t] +(lambda_0+b*T_ml[t] + c*T_ml[t]**2)*T_ml[t]-efficacy*gamma*(T_ml[t]-T_deep[t]))*delta_time/C_ml
        T_deep[t+1]    = T_deep[t] + gamma*(T_ml[t]-T_deep[t])*delta_time/C_deep
        imbalance[t+1] = (C_ml*(T_ml[t+1]-T_ml[t]) + C_deep*(T_deep[t+1]-T_deep[t]))/delta_time

        if imbalance[t+1] > 1.e4:
            T_ml[t+1:] = np.nan
            T_deep[t+1:] = np.nan
            imbalance[t+1:] = np.nan
            print('Your TOA imbalance is more than 10,000 Wm-2' + \
            'I think it\'s time to stop, don\'t you?')
            break

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
nyears        = 10000
forcing_years = np.arange(0,nyears)
input_forcing = 0*np.linspace(0.0,1.0,nyears)


qub1 = qubicon(forcing_years, input_forcing, gamma=0,ECS=-3.7/3
,c=-1.0
,sigma=30.0)


##### Periodic forcing
input_forcing = 0.8 * np.sin(np.linspace(0.0,2*np.pi,nyears))

qub2 = qubicon(forcing_years, input_forcing, gamma=0,ECS=-3.7/3
,c=-1.0
,sigma=30.0)



###### mixed layer
input_forcing = np.empty(nyears)
input_forcing[:int(nyears/2.)] = np.linspace(1.0,9./8.+0.0001,int(nyears/2.))
input_forcing[int(nyears/2.):] = 9./8.+0.0001


# No forcing
pml1 = qubicon(forcing_years, 0.*input_forcing, gamma=0,ECS=3.7/1.5
,b=1./2.
,sigma=2.42)

pml2 = qubicon(forcing_years, input_forcing, gamma=0,ECS=3.7/1.5
,b=1./2., T_ml0 = 1.0
,sigma=2.42)

# No variability
pml3 = qubicon(forcing_years, input_forcing, gamma=0,ECS=3.7/1.5
,b=1./2.,T_ml0 = 1.0
,sigma=0.)

###### Testing variability
input_forcing = np.ones(nyears)

vml1 = qubicon(forcing_years, -2.*input_forcing, gamma=0,ECS=3.7/1.5
,b=1./2., T_ml0 = -1.0
,sigma=2.42)

vml2 = qubicon(forcing_years, 0.*input_forcing, gamma=0,ECS=3.7/1.5
,b=1./2., T_ml0 = 0.0
,sigma=2.42)

vml3 = qubicon(forcing_years, 1.*input_forcing, gamma=0,ECS=3.7/1.5
,b=1./2.,T_ml0 = 1.0
,sigma=2.42)

vtl1 = qubicon(forcing_years, -2.*input_forcing, gamma=0.5,ECS=3.7/1.5
,b=1./2., T_ml0 = -1.0, T_deep0 = -1.0
,sigma=2.42)

vtl2 = qubicon(forcing_years, 0.*input_forcing, gamma=0.5,ECS=3.7/1.5
,b=1./2., T_ml0 = 0.0, T_deep0 = 0.
,sigma=2.42)

vtl3 = qubicon(forcing_years, 1.*input_forcing, gamma=0.5,ECS=3.7/1.5
,b=1./2.,T_ml0 = 1.0, T_deep0 = 1.0
,sigma=2.42)


#### Two layer
nyears        = 10000
forcing_years = np.arange(0,nyears)
input_forcing = np.empty(nyears)
input_forcing[:int(nyears/2.)] = np.linspace(1.0,1.5+0.0001,int(nyears/2.))
input_forcing[int(nyears/2.):] = 1.5+0.0001
# No forcing
ptl1 = qubicon(forcing_years, 0*input_forcing, gamma=1./2,ECS=3.7/1.5
,b=1./2.
,sigma=2.42)

ptl2 = qubicon(forcing_years, input_forcing, gamma=1./2,ECS=3.7/1.5
,b=1./2., T_ml0 = 1.0, T_deep0=0.
,sigma=2.42)

# No variability
ptl3 = qubicon(forcing_years, input_forcing, gamma=1./2,ECS=3.7/1.5
,b=1./2.,T_ml0 = 1.0, T_deep0=0.
,sigma=0.)

#--------------------------------------------------------------------------
# Plot

fig, axes = plt.subplots(1,1, figsize=(7,3.5))

color1                  = 'black'
color2                  = "royalblue"
color3                  = 'darkkhaki'
color4                  = 'salmon'
color5                  = 'orange'

lw=1.0

axes.plot(qub1['time'],qub1['T_ml'],lw=lw,color=color1,alpha=0.6)
axes.plot(qub2['time'],qub2['T_ml'],lw=lw,color=color5,alpha=0.6)


ax2 = axes.twinx()
ax2.plot(qub2['time'],qub2['forcing'],lw=lw,color=color5,alpha=1.0)
ax2.set_ylabel('Forcing')
#ax2.tick_params('y')


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
plt.savefig('../plots/Stability_qubicon.pdf', dpi=300)
plt.close()


#------------------------------------------------

# Plot

fig, axes = plt.subplots(1,1, figsize=(7,3.5))

lin1, = axes.plot(pml1['time'],pml1['T_ml'],lw=lw,color=color1,zorder=5
,label='$\sigma_{yr} = 2.42$ Wm$^{-2}$')
lin2, = axes.plot(pml3['time'],pml3['T_ml'],lw=2,color=color5,zorder=6,alpha=0.8
,label=r'F = [1, $\approx \frac{9}{8}$] Wm$^{-2}$')
lin3, = axes.plot(pml2['time'],pml2['T_ml'],lw=lw,color=color2,zorder=4,alpha=1.0
,label='F,  $\sigma_{yr}$')


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

plt.setp(axes, ylim=(-0.3,3.))
axes.legend(fontsize=10,frameon=False)

plt.tight_layout()
plt.savefig('../plots/Stability_mlm.pdf', dpi=300)
plt.close()


#------------------------------------------------

# Plot

fig, axes = plt.subplots(1,1, figsize=(7,3.5))


lin1, = axes.plot(ptl1['time'],ptl1['T_ml'],lw=lw,color=color1,zorder=5
,label='$\sigma_{yr} = 2.42$ Wm$^{-2}$')
lin2, = axes.plot(ptl3['time'],ptl3['T_ml'],lw=2,color=color5,zorder=6,alpha=0.8
,label=r'F = [1, $\approx \frac{9}{8}$] Wm$^{-2}$')
lin3, = axes.plot(ptl2['time'],ptl2['T_ml'],lw=lw,color=color2,zorder=4,alpha=1.0
,label='F,  $\sigma_{yr}$')
#lin4, = axes.plot(ptl1['time'],ptl1['T_ml'],lw=lw,color='gray',zorder=2,alpha=1.0
#,label='F,  $\sigma_{yr}$, two-layer')

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

plt.setp(axes, ylim=(-0.3,3.))
axes.legend(fontsize=10,frameon=False)

plt.tight_layout()
plt.savefig('../plots/Stability_tlm.pdf', dpi=300)
plt.close()

#------------------------------------------------

# Plot Histogram

fig, axes = plt.subplots(1,1, figsize=(7,3.5))


bins = np.arange(-1.,1.,0.02)
#print(bins)
n, obins, patches = plt.hist(vml1['T_ml'],
                            bins + np.mean(vml1['T_ml']), normed=1,
                            facecolor=color1, alpha=0.6,
                            label='T = Mixed layer')
n, obins, patches = plt.hist(vml2['T_ml'],
                            bins + np.mean(vml2['T_ml']), normed=1,
                            facecolor=color1, alpha=0.6)
n, obins, patches = plt.hist(vml3['T_ml'],
                            bins + np.mean(vml3['T_ml']), normed=1,
                            facecolor=color1, alpha=0.6)


n, obins, patches = plt.hist(vtl1['T_ml'],
                            bins + np.mean(vtl1['T_ml']), normed=1,
                            facecolor=color5, alpha=0.6,
                            label='T = Two layer')
n, obins, patches = plt.hist(vtl2['T_ml'],
                            bins + np.mean(vtl2['T_ml']), normed=1,
                            facecolor=color5, alpha=0.6)
n, obins, patches = plt.hist(vtl3['T_ml'],
                            bins + np.mean(vtl3['T_ml']), normed=1,
                            facecolor=color5, alpha=0.6)

#plt.title("Histogram with 'auto' bins")
axes.set_ylabel(r'Frequency [%]',fontsize=12)
axes.set_xlabel(r'Temperature [K]',fontsize=12)



xmax = max(axes.get_xlim())
xmin = min(axes.get_xlim())
plt.xlim((xmin,xmax))
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

axes.legend(fontsize=10,frameon=False)
plt.setp(axes, xlim=(-1.5,1.5))

plt.tight_layout()
plt.savefig('../plots/Stability_vari_1.pdf', dpi=300)
plt.close()

#----------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(6,3.5))

TT = np.linspace(-4.,4.,200)
NT = -TT**3 + 3*TT
VT = 1./4. * TT**4 - 3./2. * TT**2

NT2 = -TT**3 + 3*TT + 1
NT3 = -TT**3 + 3*TT - 1

VT2 = 1./4. * TT**4 - 3./2. * TT**2 + TT
VT3 = 1./4. * TT**4 - 3./2. * TT**2 - TT

colormlm = 'red'
colortlm = 'blue'

lin1, = ax2.plot(TT,VT, color=color1, label='F = 0',alpha=0.7)
lin3, = ax2.plot(TT,VT2, color=color2, label='F =  1',alpha=0.7)
lin5, = ax2.plot(TT,VT3, color=color5, label='F = -1',alpha=0.7)


#plt.xticks((-2,-1,1,2))
#plt.yticks((-3,-2,-1,1,2,3))

lin2, = ax1.plot(TT,NT, color=color1, label='F = 0'  ,alpha=0.7)
lin4, = ax1.plot(TT,NT2, color=color2, label='F =  1',alpha=0.7)
lin6, = ax1.plot(TT,NT3, color=color5, label='F = -1',alpha=0.7)

plt.setp(ax1, xticks=(-2,-1,1,2),yticks=(-3,-2,-1,1,2,3))
plt.setp(ax2, xticks=(-3,-2,-1,1,2,3),yticks=(-4,-2,2,4))


#plt.xticks((-3,-2,-1,1,2,3))
#plt.yticks((-4,-2,2,4))

ax2.legend(handles=[lin1,lin3,lin5], frameon=False,
fontsize=10,bbox_to_anchor=(0.92,1))

for axes in (ax1,ax2):
    spines_to_remove        = ['top', 'right']
    for spine in spines_to_remove:
        axes.spines[spine].set_visible(False)

    axes.spines['bottom'].set_position('zero')
    axes.spines['left'].set_position('zero')
    axes.yaxis.set_ticks_position('none')
    axes.xaxis.set_ticks_position('none')
    axes.xaxis.set_label_coords(1.1, -0.1)
    axes.yaxis.set_label_coords(-0.1, 0.5)

fig.text(0.5, 0.04, 'Temperature [K]', ha='center')

ax1.set_ylabel('Wm$^{-2}$')

plt.setp(ax1,title='Imbalance')
plt.setp(ax2,title='Potential')
plt.setp(ax2, xlim=(-3.,3.), ylim=(-5,5))
plt.setp(ax1, xlim=(-2.,2.), ylim=(-3.5,3.5))
plt.tight_layout()
plt.savefig('../plots/Stability_theory.pdf', dpi=300)
plt.close()
