import matplotlib.pyplot as plt
import numpy as np                
import scipy.stats as stat
from netCDF4 import Dataset

# set general text options for plotting
plt.rc('text'           , usetex=True)
plt.rc('font'           , size=12)
plt.rc('legend'         , fontsize=12)
plt.rc('text.latex'     , preamble=r'\usepackage{cmbright}')

#--------------------------------------------------------------------------
def summaryStats(samples):
        statstring   =  ' '+                                                   \
                        str('%.2f' % round(stat.scoreatpercentile(samples,50),2))+' ['  \
                       +str('%.2f' % round(stat.scoreatpercentile(samples,5),2)) +' - ' \
                       +str('%.2f' % round(stat.scoreatpercentile(samples,95),2))+'] '
        return(statstring)

#--------------------------------------------------------------------------

nsamp = int(5e6)
from95tostd = stat.norm.ppf(0.95)
almost_black            = '#262626'

#----------------------------------------------------------------------
# Load HadCRUT temperatures from file. The file as downloaded from CRU has been pre-processes with CDO:
# cdo -a settunits,years -selyear,1850/2015 -fldmean -yearmean HadCRUT.4.4.0.0.median.nc HadCRUT.4.4.0.0.median_means.nc
#
f           = Dataset('../Data/HadCRUT.4.5.0.0.median_means.nc', mode='r')
HadCRU_year = f.variables['time'][:]
HadCRU_t    = f.variables['temperature_anomaly'][:,0,0]

# Adjust the temperature relative to 1850-1899:
HadCRU_t    = HadCRU_t - np.mean(HadCRU_t[np.where((HadCRU_year >= 1850) & (HadCRU_year <= 1899))])

# Calculate temperature of the first reference period:
HadCRU_t1    = np.mean(HadCRU_t[np.where((HadCRU_year >= 1859) & (HadCRU_year <= 1882))])    

# Calculate temperature of the second reference period:
HadCRU_t2    = np.mean(HadCRU_t[np.where((HadCRU_year >= 2005) & (HadCRU_year <= 2015))])    

# Calculate temperature change to the second reference period:
HadCRU_dt   = HadCRU_t2 - HadCRU_t1

# Standard deviation on temperature change from uncertainty and variability (Lewis and Curry 2014):
HadCRU_sd   = 0.08

# Randomly sample temperatures:
delta_t     = np.random.normal(loc=HadCRU_dt, scale=HadCRU_sd, size=nsamp)

#----------------------------------------------------------------------
# Radiation imbalance for 2005-2015 from Johnson et al. (2016), doi:10.1038/nclimate3043:
Q2_mean = 0.71
Q2_sd   = 0.10/from95tostd

# Imbalance for 1859-1882 from Lewis and Curry (2014), based on GCM:
Q1_mean = 0.15
Q1_sd   = 0.075

# Ramdomly sample imbalance and imbalance difference:
imbalance = np.random.normal(loc=Q2_mean, scale=Q2_sd, size=nsamp)
delta_imbalance = imbalance - np.random.normal(loc=Q1_mean, scale=Q1_sd, size=nsamp)

#----------------------------------------------------------------------
# Load AR5 forcings from file:
F_year, F_co2, F_ghgother, F_o3trop, F_o3strat, F_aero, F_luc, F_h2ostrat, F_bcsnow, F_contrails, F_solar, F_volc = \
  np.loadtxt('../Data/WG1_AR5_Annex_II_sheet_1-2_forcings.txt', skiprows=1, unpack=True)

# Calculate total forcing, and assume volcanoes are on average zero:
F_total = F_co2 + F_ghgother + F_o3trop + F_o3strat + F_aero + F_luc + F_h2ostrat + F_bcsnow + F_contrails + F_solar + F_volc
F_total = F_total - np.mean(F_volc)

# Identify indices (i1 and i2) for the reference periods. For the later 
# reference period of 2005-2015 we choose the values in the year 2010 as 
# the average forcing representative of the period.
i1 = np.where((F_year >= 1859) & (F_year <= 1882))
i2 = np.where(F_year == 2010)

# Forcing from a doubling of CO2 that is consistent with :
F2x_mean               =   3.71

# Mean change in forcings between reference periods:
delta_F_ghg_mean       =   np.mean(F_co2[i2])-np.mean(F_co2[i1]) \
                         + np.mean(F_ghgother[i2])-np.mean(F_ghgother[i1])
delta_F_o3_mean        =   np.mean(F_o3trop[i2])-np.mean(F_o3trop[i1]) \
                         + np.mean(F_o3strat[i2])-np.mean(F_o3strat[i1])
delta_F_aero_mean      =   np.mean(F_aero[i2])-np.mean(F_aero[i1])
F_aero_mean            =   np.mean(F_aero[i2])
delta_F_luc_mean       =   np.mean(F_luc[i2])-np.mean(F_luc[i1])
delta_F_h2ostrat_mean  =   np.mean(F_h2ostrat[i2])-np.mean(F_h2ostrat[i1])
delta_F_bcsnow_mean    =   np.mean(F_bcsnow[i2])-np.mean(F_bcsnow[i1])
delta_F_contrails_mean =   np.mean(F_contrails[i2])-np.mean(F_contrails[i1])
delta_F_nat_mean       =   np.mean(F_solar[i2])-np.mean(F_solar[i1]) \
                         + np.mean(F_volc[i2])-np.mean(F_volc[i1])

delta_F_mean           =   delta_F_ghg_mean + delta_F_o3_mean + delta_F_aero_mean \
                         + delta_F_luc_mean + delta_F_h2ostrat_mean + delta_F_bcsnow_mean \
                         + delta_F_contrails_mean + delta_F_nat_mean 

#----------------------------------------------------------------------
    
    
#----------------------------------------------------------------------
# Define standard deviations on all forcings. These are taken from
# AR5 wherever possible.

# Well mixed greenhouse gases: 
# "The RF of WMGHG is 2.83 (2.54 to 3.12)"
ghg_sd = (3.12 - 2.54)/(2*from95tostd)

# Forcing from CO2-doubling alone, 20 percent as in Lewis and Curry (2014):
F2x_sd = ghg_sd*F2x_mean/delta_F_ghg_mean

# Ozone:
# " total RF estimated from modelled ozone changes is 0.35 (0.15 to 0.55) W m–2, 
#   with RF due to tropospheric ozone changes of 0.40 (0.20 to 0.60) W m–2 and 
#   due to stratospheric ozone changes of –0.05 (–0.15 to +0.05) W m–2" 
o3_sd = (0.55 - 0.14)/(2*from95tostd)
    
# Aerosols: 
# "is estimated as ERF of –0.9 (–1.9 to –0.1) W m–2."
aero_sd = (1.9 - 0.1)/(2*from95tostd)
    
# Land use change:  
# "leads to an RF of –0.15 ± 0.10 W m–2", here we assume this means 5-95 percentiles
luc_sd = 0.1/from95tostd 

# Stratospheric water vapor:
# "is 0.07 (0.02 to 0.12) W m–2" 
h2ostrat_sd = (0.12 - 0.02)/(2*from95tostd)

# Black carbon (BC) on snow and ice:
# "is 0.04 (0.02 to 0.09) W m–2. "
bcsnow_sd = (0.09 - 0.02)/(2*from95tostd)

# Contrails: no uncertainty range found.
contrails_sd = 0.0

# Natural forcings (volcanic and solar): no uncertainty range found.
nat_sd = 0.0

#----------------------------------------------------------------------
# Randomly sample forcings as Gaussian distributions:

delta_F_aero    =  np.random.normal(loc=delta_F_aero_mean, scale=aero_sd, size=nsamp)
total_F_aero    =  delta_F_aero + (F_aero_mean-delta_F_aero_mean)
delta_F_ghg     =  np.random.normal(loc=delta_F_ghg_mean, scale=ghg_sd, size=nsamp)
delta_F_nonaero =  np.random.normal(loc=delta_F_o3_mean, scale=o3_sd, size=nsamp) \
                  +np.random.normal(loc=delta_F_luc_mean, scale=luc_sd, size=nsamp) \
                  +np.random.normal(loc=delta_F_h2ostrat_mean, scale=h2ostrat_sd, size=nsamp) \
                  +np.random.normal(loc=delta_F_bcsnow_mean, scale=bcsnow_sd, size=nsamp) \
                  +delta_F_contrails_mean + delta_F_nat_mean
delta_F         =  delta_F_aero + delta_F_ghg + delta_F_nonaero

# Make errors in forcing from CO2 doubling correlate with errors in greenhouse gas forcing:
F2x             = (delta_F_ghg-delta_F_ghg_mean)*F2x_sd/ghg_sd + F2x_mean

print(' ')
print('----------------------------------------')
print('GHG forcing:         '+summaryStats(delta_F_ghg))
print('Aerosol forcing:    ' +summaryStats(delta_F_aero))
print('Non-aerosol forcing: '+summaryStats(delta_F_nonaero))
print('----------------------------------------')
print('Total forcing change:'+summaryStats(delta_F))
print(' ')

print('----------------------------------------')
print('Aerosol forcing uncertainty is '\
          + str('%.2f' % round(100.*np.std(total_F_aero)**2/(np.std(delta_F - delta_imbalance)**2),2)) \
          +'% of the variance in the denominator of ECS')


#----------------------------------------------------------------------
# Calculate ECS and TCR from samples:

TCR  = F2x*delta_t/delta_F
ECS  = F2x*delta_t/(delta_F-delta_imbalance)

# Set infinite values to not any number:
id_inf      = np.where((TCR<0.0) | (ECS<0.0))
TCR[id_inf] = float('nan')
ECS[id_inf] = float('nan')
id          = np.where(np.isfinite(ECS))

print('----------------------')
print('Climate sensitivities:')
print('TCR:'+summaryStats(TCR[id]))
print('ECS:'+summaryStats(ECS[id]))
print('----------------------')
print(' ')



#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Plots
#----------------------------------------------------------------------
#----------------------------------------------------------------------
    
colora = 'black'
colorb = 'red'
colorc = 'purple'
colord = 'deepskyblue'
colore = 'orange'


# Calculate histograms:
bins        = np.arange(0,10,0.02)
bin_centers = (bins[1:]+bins[:-1])/2.
        
pTCR, hTCR = np.histogram(TCR[id],bins)
pECS, hECS = np.histogram(ECS[id],bins)
ratio      = TCR[id]/ECS[id]
pratio, hratio = np.histogram(ratio,np.arange(0,1,0.005))
    
#----------------------------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(5,3))


axes.plot(bin_centers, pTCR, color='deepskyblue', label='TCR')
axes.plot(bin_centers, pECS, color='black', label='ECS')

axes.set_xlim([0, 6])
xmax = max(axes.get_xlim())
ymax = max(axes.get_ylim())
axes.set_ylim([0, ymax])

axes.set_xlabel('Temperature [K]')
axes.set_ylabel('Probability')

axes.text(xmax*0.8,ymax*0.7,'ECS: '+summaryStats(ECS[id]), color='black', va='top', ha='right', fontsize=10)
axes.text(xmax*0.8,ymax*0.61, 'TCR: '+summaryStats(TCR[id]), color='deepskyblue', va='top', ha='right', fontsize=10)


for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

spines_to_remove        = ['top', 'right'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

axes.set_yticks([])
axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')

spines_to_keep = [ 'bottom', 'left']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

plt.tight_layout()
plt.savefig('../plots/ECS_inferred.pdf', dpi=300)
plt.close()

#----------------------------------------------------------------------

fig, axes = plt.subplots(1,1, figsize=(5,3))

# Calculate histograms:
bins        = np.arange(-5,5,0.02)
bin_centers = (bins[1:]+bins[:-1])/2.

pF, hF = np.histogram(delta_F,bins)
pF_aero, hF_aero = np.histogram(delta_F_aero,bins)
pF_naero, hF_nonaero = np.histogram(delta_F-delta_F_aero,bins)

#axes.plot(bin_centers, pTCR, color='deepskyblue', label='TCR')
axes.plot(bin_centers, pF, color='black',lw=1.5)
axes.plot(bin_centers, pF_aero, color='red')
axes.plot(bin_centers, pF_naero, color='green')

axes.set_xlim([-3, 5])
xmin = min(axes.get_xlim())
ymax = max(axes.get_ylim())
axes.set_ylim([0, ymax])

axes.set_xlabel(r'Forcing (Wm$^{-2}$)')
#axes.set_ylabel('Probability')

#axes.text(xmax*0.8,ymax*0.7,'ECS: '+summaryStats(ECS[id]), color='black', va='top', ha='right', fontsize=10)
axes.text(xmin,ymax*0.9, 'Total forcing', color='black', va='top', ha='left', fontsize=10)
axes.text(xmin,ymax*0.8, 'Non-aerosol forcing', color='green', va='top', ha='left', fontsize=10)
axes.text(xmin,ymax*0.7, 'Aerosol forcing', color='red', va='top', ha='left', fontsize=10)


for ticks in axes.xaxis.get_ticklines() + axes.yaxis.get_ticklines():
    ticks.set_color(almost_black)

axes.spines['left'].set_position('zero')

spines_to_remove        = ['top', 'right'] 
for spine in spines_to_remove:
    axes.spines[spine].set_visible(False)

axes.set_yticks([])
axes.yaxis.set_ticks_position('left')
axes.xaxis.set_ticks_position('bottom')

spines_to_keep = [ 'bottom', 'left']     
for spine in spines_to_keep:
    axes.spines[spine].set_linewidth(0.5)
    axes.spines[spine].set_color(almost_black)

plt.tight_layout()
plt.savefig('../plots/ECS_inferred_forcing.pdf', dpi=300)
plt.close()
