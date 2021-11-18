#!/usr/bin/env python
# coding: utf-8

# ------------------------------------------------------------------------
#
# TITLE - make_solar_samples.py
# AUTHOR - James Lane
# PROJECT - mw-dfs
#
# ------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Draw samples from the DFs at the solar position. Creating the solar sample
from Lane+ 2021
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
data_dir = '/geir_data/scr/lane/projects/mw-dfs/data/data_Sept_2021/gaia_apogee_processed/'
out_dir = '/geir_data/scr/lane/projects/mw-dfs/data/data_Sept_2021/df_samples/'
# Ensure directories exist
if not os.path.exists(out_dir+'df_stash/'):
    os.makedirs(out_dir+'df_stash/')
##fi

# Sampling
n_samples = int(1e5)
betas = [0.,-0.5,0.5,0.7,0.9]

# Make which of the 4 datasets
make_solar_halo = True
make_solar_disk = False

# MW Potential and interpolated spherical MW Potential
rmin = 1/ro # 1 kpc
rmax = 60/ro # 60 kpc
rmin_interp = rmin/2.
rmax_interp = rmax*2.
ngrid = 401
mwpot = potential.MWPotential2014
print('Making interpolated Milky Way Potential')
interpot = project_potential.make_interpolated_mwpot(mwpot=mwpot,
    rmin=rmin_interp,rmax=rmax_interp, ngrid=ngrid, ro=ro, vo=vo, match_type='mass')
potential.turn_physical_on(mwpot,ro=ro,vo=vo)
potential.turn_physical_on(interpot,ro=ro,vo=vo)

# Halo density profile potential
alpha = 3.5 # halo density inner power law slope
rc = 30*apu.kpc
denspot = potential.PowerSphericalPotentialwCutoff(amp=1., r1=1.,
    alpha=alpha, rc=rc, ro=ro, vo=vo)
potential.turn_physical_on(denspot,ro=ro,vo=vo)

### Sample Halo

print('Making solar halo samples')
halo_orbs = []
for i in range(len(betas)):
    # Catch warnings otherwise the program is loud
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        print('Sampling beta: '+str(betas[i]))
        # Loading DFs takes awhile, so check if they're stashed
        df_name = out_dir+'df_stash/constantbetadf_beta_'+str(betas[i])+\
            '_rmax_'+str(round(rmax))+'_pot_'+interpot.__class__.__name__+\
            '_denspot_'+denspot.__class__.__name__+'.pkl'
        if os.path.exists(df_name):
            loaded_df = True
            print('Loading DF from '+df_name)
            with open(df_name,'rb') as f:
                dfm = pickle.load(f)
        else:
            loaded_df = False
            dfm = df.constantbetadf(pot=interpot, denspot=denspot, ro=ro, vo=vo, 
                beta=betas[i], rmax=rmax)
        ##ie
        o = dfm.sample(R=np.ones(n_samples)*ro*apu.kpc, 
            phi=np.zeros(n_samples), z=np.zeros(n_samples), rmin=rmin)
        halo_orbs.append(o)
        print('Target beta: '+str(betas[i]))
        print('Measured beta: '+str(round(project_df.calculate_beta_one_r(o),3))+'\n')
        if not loaded_df:
            print('Saving DF to '+df_name)
            with open(df_name,'wb') as f:
                pickle.dump(dfm,f)
            ##wi
        ##fi
    ##wi
###i

# Calculate the conserved quantities
print('Calculating kinematic quantities')
halo_actions, halo_eELzs = project_df.calculate_conserved_quantities(halo_orbs, 
    mwpot, aA_type='staeckel', ro=ro, vo=vo, zo=zo)

### Save Halo

with open(out_dir+'halo_orbits_solar.pkl','wb') as f:
    pickle.dump(halo_orbs,f)
##wi
with open(out_dir+'halo_actions_solar.pkl','wb') as f:
    pickle.dump(halo_actions,f)
##wi
with open(out_dir+'halo_eELzs_solar.pkl','wb') as f:
    pickle.dump(halo_eELzs,f)
##w

### Sample Disk

print('Making solar disk samples')
potential.turn_physical_off(mwpot)
# Using actionAngle adiabatic for stability, but produces same results as 
# Staeckel
aAA = aA.actionAngleAdiabatic(pot=mwpot, c=True, ro=ro, vo=vo, zo=zo, use_physical=True)

# Make the DFs
sigma_vR, sigma_vz  = 40*apu.km/apu.s, 30*apu.km/apu.s
thin_disk_Rscale, thin_disk_vRscale, thin_disk_vzscale = 3*apu.kpc, 8*apu.kpc, 8*apu.kpc
thin_disk_df = df.quasiisothermaldf(hr=thin_disk_Rscale,sr=sigma_vR,sz=sigma_vz,
    hsr=thin_disk_vRscale,hsz=thin_disk_vzscale, pot=mwpot, aA=aAA, ro=ro, vo=vo)
thick_disk_df = df.quasiisothermaldf(hr=thin_disk_Rscale*(2/3),sr=sigma_vR*1.5,
    sz=sigma_vz*1.5,hsr=thin_disk_vRscale, hsz=thin_disk_vzscale, 
    pot=mwpot,aA=aAA, ro=ro, vo=vo)

# Draw samples
thin_disk_samples = thin_disk_df.sampleV(R=ro/ro, z=zo/ro, n=n_samples)
thick_disk_samples = thick_disk_df.sampleV(R=ro/ro, z=zo/ro, n=n_samples)

# Make into orbits
vxvv_thin = np.array([np.ones(n_samples)*ro/ro, thin_disk_samples[:,0].value/vo,
    thin_disk_samples[:,1].value/vo, np.ones(n_samples)*zo/ro, 
    thin_disk_samples[:,2].value/vo, np.zeros(n_samples) ]).T
vxvv_thick = np.array([np.ones(n_samples)*ro/ro, thick_disk_samples[:,0].value/vo,
    thick_disk_samples[:,1].value/vo, np.ones(n_samples)*zo/ro, 
    thick_disk_samples[:,2].value/vo, np.zeros(n_samples) ]).T
os_thin = orbit.Orbit( vxvv=vxvv_thin, ro=ro, vo=vo, zo=zo )
os_thick = orbit.Orbit( vxvv=vxvv_thick, ro=ro, vo=vo, zo=zo  )
disk_orbs = [os_thin,os_thick]

# Calculating kinematic quantities
print('Calculating kinematic quantities')
disk_actions, disk_eELzs = project_df.calculate_conserved_quantities(disk_orbs, 
    mwpot, aA_type='staeckel', ro=ro, vo=vo, zo=zo)
 
### Save Disk

with open(out_dir+'disk_orbits_solar.pkl','wb') as f:
    pickle.dump(disk_orbs,f)
##wi
with open(out_dir+'disk_actions_solar.pkl','wb') as f:
    pickle.dump(disk_actions,f)
##wi
with open(out_dir+'disk_eELzs_solar.pkl','wb') as f:
    pickle.dump(disk_eELzs,f)
##wi