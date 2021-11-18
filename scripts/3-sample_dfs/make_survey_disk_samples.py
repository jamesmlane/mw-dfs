#!/usr/bin/env python
# coding: utf-8

# ------------------------------------------------------------------------
#
# TITLE - make_solar_disk_samples
# AUTHOR - James Lane
# PROJECT - mw-dfs
#
# ------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Draw samples from the disk DFs at the positions of APOGEE stars. Creating 
the disk portion of the survey sample from Lane+ 2021
'''

__author__ = "James Lane"

### Imports

## Basic
import numpy as np, sys, os, warnings, dill as pickle
import astropy.units as apu

## galpy
from galpy import potential
from galpy import orbit
from galpy import actionAngle as aA
from galpy.util import conversion
from galpy import df
ro = 8.178 # Gravity 2019
vo = 220
zo = 0.0208 # Bennet + Bovy 2018

sys.path.append('../../src/')
from mw_dfs import df as project_df
from mw_dfs import potential as project_potential

### Keywords & Setup

# Directories
data_dir = '../../data/gaia_apogee_processed/'
out_dir = '../../data/df_samples/'
# Ensure directories exist
if not os.path.exists(out_dir+'df_stash/'):
    os.makedirs(out_dir+'df_stash/')
##fi

# Sampling
oversamp_fac = 10 # Factor by which to oversample the disk survey set
print('Oversample factor is: '+str(oversamp_fac))

# Apply Gaia uncertainties to survey disk orbits?
survey_disk_gaia_errs = True

# MW Potential
mwpot = potential.MWPotential2014

### Load APOGEE Data

print('Loading Gaia-APOGEE disk data')
gaia2_thin = np.load(data_dir+'gaia_data_thin.npy', allow_pickle=True)
allstar_thin = np.load(data_dir+'allstar_thin.npy', allow_pickle=True)
with open(data_dir+'orbits_thin.pkl','rb') as f:
    os_thin = pickle.load(f)    
##wi
gaia2_thick = np.load(data_dir+'gaia_data_thick.npy', allow_pickle=True)
allstar_thick = np.load(data_dir+'allstar_thick.npy', allow_pickle=True)
with open(data_dir+'orbits_thick.pkl','rb') as f:
    os_thick = pickle.load(f)
##wi

# Only pick needed fields to cut down on size of files being moved around
gaia2_req_fields = ['ra_error','dec_error','pmra_error','pmdec_error',
                    'ra_dec_corr','ra_pmra_corr','ra_pmdec_corr',
                    'dec_pmra_corr','dec_pmdec_corr','pmra_pmdec_corr']
allstar_req_fields = ['weighted_dist_error','VERR']
gaia2_thin = gaia2_thin[gaia2_req_fields]
gaia2_thick = gaia2_thick[gaia2_req_fields]
allstar_thin = allstar_thin[allstar_req_fields]
allstar_thick = allstar_thick[allstar_req_fields]

# Make into galpy orbits
os_thin = orbit.Orbit(vxvv=np.array([os_thin.R().value/ro, os_thin.vR().value/vo, 
    os_thin.vT().value/vo, os_thin.z().value/ro, os_thin.vz().value/vo, 
    os_thin.phi().value]).T, ro=ro, vo=vo, zo=zo)
os_thick = orbit.Orbit(vxvv=np.array([os_thick.R().value/ro, os_thick.vR().value/vo, 
    os_thick.vT().value/vo, os_thick.z().value/ro, os_thick.vz().value/vo, 
    os_thick.phi().value]).T, ro=ro, vo=vo, zo=zo)

### Sampling

print('Making survey disk samples')
potential.turn_physical_off(mwpot)
# Using actionAngle adiabatic for stability, but produces same results as 
# Staeckel
aAA = aA.actionAngleAdiabatic(pot=mwpot, c=True, ro=ro, vo=vo, zo=zo, use_physical=True)

# Make the DFs
sigma_vR, sigma_vz  = 40*apu.km/apu.s, 30*apu.km/apu.s
thin_disk_Rscale, thin_disk_vRscale, thin_disk_vzscale = 3*apu.kpc, ro*apu.kpc, ro*apu.kpc
thin_disk_df = df.quasiisothermaldf(hr=thin_disk_Rscale,sr=sigma_vR,sz=sigma_vz,
    hsr=thin_disk_vRscale,hsz=thin_disk_vzscale, pot=mwpot, aA=aAA, ro=ro, vo=vo)
thick_disk_df = df.quasiisothermaldf(hr=thin_disk_Rscale*(2/3),sr=sigma_vR*1.5,
    sz=sigma_vz*1.5,hsr=thin_disk_vRscale, hsz=thin_disk_vzscale, 
    pot=mwpot,aA=aAA, ro=ro, vo=vo)

# Grab variables here so we don't have to keep querying
os_thin_R = os_thin.R().repeat(oversamp_fac)
os_thin_z = os_thin.z().repeat(oversamp_fac)
os_thin_phi = os_thin.phi().repeat(oversamp_fac)
os_thick_R = os_thick.R().repeat(oversamp_fac)
os_thick_z = os_thick.z().repeat(oversamp_fac)
os_thick_phi = os_thick.phi().repeat(oversamp_fac)

# Draw samples
thin_disk_samples = thin_disk_df.sampleV_interpolate(R=os_thin_R.value/ro, 
    z=os_thin_z.value/ro, R_pixel=0.01, z_pixel=0.01)
thick_disk_samples = thick_disk_df.sampleV_interpolate(R=os_thick_R.value/ro, 
    z=os_thick_z.value/ro, R_pixel=0.01, z_pixel=0.01, R_min=0.1, R_max = 2, 
    z_max = 0.5 )

# Make into orbits
vxvv_thin = np.array([os_thin_R.value/ro, thin_disk_samples[:,0].value/vo, 
    thin_disk_samples[:,1].value/vo, os_thin_z.value/ro, 
    thin_disk_samples[:,2].value/vo, os_thin_phi.value])
vxvv_thick = np.array([os_thick_R.value/ro, thick_disk_samples[:,0].value/vo, 
    thick_disk_samples[:,1].value/vo, os_thick_z.value/ro, 
    thick_disk_samples[:,2].value/vo, os_thick_phi.value])
os_thin_noerr = orbit.Orbit(vxvv=vxvv_thin.T, ro=ro, vo=vo, zo=zo )
os_thick_noerr = orbit.Orbit(vxvv=vxvv_thick.T, ro=ro, vo=vo, zo=zo )
disk_orbs_noerr = [os_thin_noerr,os_thick_noerr]

# Now apply observational uncertainties to orbits
disk_orbs = []
gaia2_disk = [gaia2_thin.repeat(oversamp_fac), gaia2_thick.repeat(oversamp_fac)]
allstar_disk = [allstar_thin.repeat(oversamp_fac), allstar_thick.repeat(oversamp_fac)]
for i in range(len(disk_orbs_noerr)):
    o = project_df.sample_orbits_with_uncertainties(disk_orbs_noerr[i], 
        gaia2_disk[i], allstar_disk[i], ro=ro, vo=vo, zo=zo, only_velocities=True)
    disk_orbs.append(o)
###i

# Calculating kinematic quantities
print('Calculating kinematic quantities')
disk_actions, disk_eELzs = project_df.calculate_conserved_quantities(disk_orbs, 
    mwpot, aA_type='staeckel', ro=ro, vo=vo, zo=zo)

### Save Disk

oversample_suffix = '_os_factor'+str(oversamp_fac)+'.pkl'    
with open(out_dir+'disk_orbits_survey'+oversample_suffix,'wb') as f:
    pickle.dump(disk_orbs,f)
##wi
with open(out_dir+'disk_actions_survey'+oversample_suffix,'wb') as f:
    pickle.dump(disk_actions,f)
##wi
with open(out_dir+'disk_eELzs_survey'+oversample_suffix,'wb') as f:
    pickle.dump(disk_eELzs,f)
##wi