#!/usr/bin/env python
# coding: utf-8

# ------------------------------------------------------------------------
#
# TITLE - make_survey_halo_ig_samples.py
# AUTHOR - James Lane
# PROJECT - mw-dfs
#
# ------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Draw samples from the halo DFs at the positions of APOGEE stars in the 
inner Galaxy. Creating the IG halo portion of the survey sample from Lane+ 2021
'''

__author__ = "James Lane"

### Imports

## Basic
import numpy as np, pdb, sys, os, copy, warnings, dill as pickle
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
oversamp_fac = 100 # Factor by which to oversample the halo survey set
print('Oversample factor is: '+str(oversamp_fac))
betas = [0.,-0.5,0.5,0.7,0.9]

# Apply Gaia uncertainties to survey halo orbits?
survey_halo_gaia_errs = True

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


### Load APOGEE data

print('Loading Gaia-APOGEE halo data')
gaia2_halo = np.load(data_dir+'gaia_data_halo_ig.npy', allow_pickle=True)
allstar_halo = np.load(data_dir+'allstar_halo_ig.npy', allow_pickle=True)
with open(data_dir+'orbits_halo_ig.pkl','rb') as f:
    os_halo = pickle.load(f)
##wi
os_halo = orbit.Orbit(vxvv=np.array([os_halo.R().value/ro, 
    os_halo.vR().value/vo, os_halo.vT().value/vo, os_halo.z().value/ro, 
    os_halo.vz().value/vo, os_halo.phi().value]).T, ro=ro, vo=vo, zo=zo)

### Sampling

print('Making survey halo samples')
halo_orbs_noerr = []
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
        o = dfm.sample(R=os_halo.R().repeat(oversamp_fac),
            phi=os_halo.phi().repeat(oversamp_fac), 
            z=os_halo.z().repeat(oversamp_fac))
        halo_orbs_noerr.append(o)
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

# Now apply observational uncertainties to orbits
print('Applying observational uncertainties')
halo_orbs = []
for i in range(len(halo_orbs_noerr)):
    o = project_df.sample_orbits_with_uncertainties(halo_orbs_noerr[i], 
        gaia2_halo.repeat(oversamp_fac), allstar_halo.repeat(oversamp_fac),
        ro=ro, vo=vo, zo=zo, only_velocities=True)
    halo_orbs.append(o)
###i

# Calculating kinematic quantities
print('Calculating kinematic quantities')
halo_actions, halo_eELzs = project_df.calculate_conserved_quantities(halo_orbs, 
    mwpot, aA_type='staeckel', ro=ro, vo=vo, zo=zo)

### Save survey samples

oversample_suffix = '_os_factor'+str(oversamp_fac)+'.pkl'
with open(out_dir+'halo_orbits_ig_survey'+oversample_suffix,'wb') as f:
    pickle.dump(halo_orbs,f)
##wi
# noerr orbits just for validation of error addition if interested
with open(out_dir+'halo_orbits_ig_survey_noerr'+oversample_suffix,'wb') as f:
    pickle.dump(halo_orbs_noerr,f)
##wi
with open(out_dir+'halo_actions_ig_survey'+oversample_suffix,'wb') as f:
    pickle.dump(halo_actions,f)
##wi
with open(out_dir+'halo_eELzs_ig_survey'+oversample_suffix,'wb') as f:
    pickle.dump(halo_eELzs,f)
##wi
