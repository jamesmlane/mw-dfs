# ----------------------------------------------------------------------------
#
# TITLE - df.py
# PROJECT - mw-dfs
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Classes and functions for dealing with distribution functions
'''
__author__ = "James Lane"

### Imports
import numpy as np, sys, os, copy, warnings
import pdb
import astropy.units as apu
from galpy import potential
from galpy import orbit
from galpy import actionAngle as aA
import scipy.interpolate
import scipy.special
import scipy.integrate

# ----------------------------------------------------------------------------

def calculate_beta_one_r(os):
    '''calculate_beta_one_r:
    
    Calculate beta for orbits at a single radius
    
    Args:
        os (galpy.orbit.Orbit) - Orbit object
    '''
    vr = os.vr()
    vphi = os.vT()
    vtheta = os.vtheta()
    beta = 1 - (np.square(np.std(vphi)) + np.square(np.std(vtheta)))\
        /(2*np.square(np.std(vr)))
    if type(beta) is apu.quantity.Quantity:
        return beta.value
    else:
        return beta
    ##ie
#def

def calculate_conserved_quantities(orbs,pot,aA_type,ro,vo,zo):
    '''calculate_conserved_quantities:
    
    Calculate eccentricity, Energy, Angular momentum, along with 2 actions (+ Lz)
    
    Args:
        orbs (list) - list of galpy.orbit.Orbit objects
        pot (galpy.potential.Potential) - Potential for calculating Energy
        aA_type (string) - Type of actionAngle object to use, 'staeckel' or 
            'spherical' ['staeckel']
        ro,vo,zo (float) - galpy scale units
    '''
    # Make sure we are outputing in astropy units
    potential.turn_physical_on(pot,ro=ro,vo=vo)
    
    accs_out = []
    eELzs_out = []
    isUnbound_out = []
    
    if aA_type == 'staeckel':
        # Delta of 0.4 is temporary place-holder
        aA_use = aA.actionAngleStaeckel(pot=pot,c=True,delta=0.4,ro=ro,vo=vo,zo=zo)
    elif aA_type == 'spherical_staeckel':
        # Way faster than using Spherical. Must have delta>0 so 1e-8
        aA_use = aA.actionAngleStaeckel(pot=pot,c=True,delta=1e-8,ro=ro,vo=vo,zo=zo)
        
    # Ignore all the integration warnings.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        
        # Loop over each of the sets of orbits
        for i,os in enumerate(orbs):
            
            # Arrays for holding output of each set of orbit.Orbit objects
            eELz = np.zeros((3,len(os)))
            actions = np.zeros((3,len(os)))
            isUnbound = np.zeros(len(os),dtype=bool)
            
            # if we're using Staeckel then calculate deltas
            if aA_type == 'staeckel':
                print('\nComputing deltas..')
                # Add a small offset so delta estimations don't freak out when z=0
                staeckelOffset = (os.z().value==0).astype(int)*1e-8*apu.kpc
                deltas = aA.estimateDeltaStaeckel(pot, os.R(),
                    os.z()+staeckelOffset, no_median=True)
            ##fi
            
            # Do orbit kinematics individually to catch anything unbound
            for j,o in enumerate(os):
                
                # Print progress
                outstr = '\rOrbits '+str(i+1)+' / '+str(len(orbs))+\
                         ', Doing orbit '+str(j+1)+' / '+str(len(os))
                sys.stdout.write(outstr)
                sys.stdout.flush()
                
                # Try-catch the action / kinematic evaluation
                try:
                    if aA_type == 'staeckel':
                        actions_temp = aA_use(o,c=True,
                            use_physical=True,delta=deltas[j].value/ro, 
                            ro=ro,vo=vo,zo=zo)[:3]
                    elif aA_type == 'spherical_staeckel':
                        actions_temp = aA_use(o,c=True,
                            use_physical=True,ro=ro,vo=vo,zo=zo)[:3]
                        
                    actions[0,j] = actions_temp[0].value
                    actions[1,j] = actions_temp[1].value
                    actions[2,j] = actions_temp[2].value
                    eELz[0,j] = aA_use.EccZmaxRperiRap(o)[0]
                except aA.UnboundError:
                    print('Star index '+str(j)+' is unbound!')
                    actions[:,j] = (np.NaN,np.NaN,np.NaN)
                    eELz[0,j] = np.NaN
                    isUnbound[j] = True
                ##te
                eELz[1,j] = o.E(pot=pot,use_physical=True,ro=ro,vo=vo,zo=zo).value
                eELz[2,j] = o.Lz(use_physical=True,ro=ro,vo=vo,zo=zo).value
            ###o
            accs_out.append(actions)
            eELzs_out.append(eELz)
            isUnbound_out.append(isUnbound)
        ##os
    ##wi
    return accs_out,eELzs_out
#def

# def bootstrap_samples(n,gaia2,allstar,os):
#     assert len(gaia2) == len(allstar) == len(os)
#     n_sample = len(gaia2)
#     if n_sample >= n:
#         bootstrap_inds = np.random.choice(n_sample,size=n,replace=False)
#     else:
#         bootstrap_inds = np.random.choice(n_sample,size=n,replace=True)
#     ##ie
#     return gaia2[bootstrap_inds], allstar[bootstrap_inds], os[bootstrap_inds]
# #def

def sample_orbits_with_uncertainties(os, gaia2, allstar, ro, vo, zo, only_velocities=False):
    '''sample_orbits_with_uncertainties:
    
    Take os, gaia data, and allstar in order to make a new set of orbits 
    with uncertainties correctly handled os should correspond directly to 
    gaia2 and allstar. i.e. os.dist() should be the same as 
    allstar['weighted_dist'] etc...
    
    Args:
        os (galpy.orbit.Orbit object) - Original orbits
        gaia2 (array) - Gaia data
        allstar (array) - APOGEE data
        ro,vo,zo (float) - unit scales
        only_velocities (bool) - Only sample the velocities (i.e. if the positions are already impacted by uncertainties)
    
    Returns:
        os_sample (galpy.orbit.Orbit object) - New orbits with uncertainties
    '''
    if only_velocities:
        print('Warning, only perturbing velocities with uncertainties')
    n_samples = len(os)
    ra_sample,dec_sample,pmra_sample,pmdec_sample,dist_sample,rv_sample = np.zeros((6,n_samples))
    
    # Positions / velocities from the DF-sampled orbit
    ra,dec = os.ra().value, os.dec().value, 
    pmra,pmdec = os.pmra().value, os.pmdec().value
    dist,rv = os.dist().value, os.vlos().value 
    
    # Uncertainties. Note 'weighted_dist_error' is in pc (want kpc)
    ra_e,dec_e = gaia2['ra_error'], gaia2['dec_error']
    pmra_e,pmdec_e = gaia2['pmra_error'], gaia2['pmdec_error']
    dist_e,rv_e = allstar['weighted_dist_error']/1e3, allstar['VERR']
    
    if only_velocities:
        ra_e = dec_e = dist_e = np.zeros_like(ra_e)
    ##fi
    
    # Get covariances. Note that nothing covaries with radial velocity
    radec, radist, rapmra, rapmdec = gaia2['ra_dec_corr'], 0., gaia2['ra_pmra_corr'], gaia2['ra_pmdec_corr']
    decdist, decpmra, decpmdec = 0., gaia2['dec_pmra_corr'], gaia2['dec_pmdec_corr']
    distpmra, distpmdec = 0., 0.
    pmrapmdec = gaia2['pmra_pmdec_corr']
    
    if only_velocities:
        radec = radist = rapmra = rapmdec = np.zeros_like(radec)
        decdist = decpmra = decpmdec = np.zeros_like(decdist)
        distpmra = distpmdec = np.zeros_like(distpmra)
    ##fi
    
    # Covariance matrix
    cov = np.zeros((n_samples,6,6))
    rv_zeroed = np.zeros(n_samples)
    cov[:,0,:]  = np.dstack([ra_e**2, ra_e*dec_e*radec, ra_e*dist_e*radist, 
        ra_e*pmra_e*rapmra, ra_e*dec_e*rapmdec, rv_zeroed])[0] 
    cov[:,1,1:] = np.dstack([dec_e**2, dec_e*dist_e*decdist, 
        dec_e*pmra_e*decpmra, dec_e*pmdec_e*decpmdec, rv_zeroed])[0]
    cov[:,2,2:] = np.dstack([dist_e**2, dist_e*pmra_e*distpmra, 
        dist_e*pmdec_e*distpmdec, rv_zeroed])[0]
    cov[:,3,3:] = np.dstack([pmra_e**2, pmra_e*pmdec_e*pmrapmdec, 
        rv_zeroed])[0]
    cov[:,4,4:] = np.dstack([pmdec_e**2, rv_zeroed])[0]
    cov[:,5,5] = rv_e**2
    
    # Matrix is symmetric
    cov[:,:,0] = cov[:,0,:]
    cov[:,1:,1] = cov[:,1,1:]
    cov[:,2:,2] = cov[:,2,2:]
    cov[:,3:,3] = cov[:,3,3:]
    cov[:,4:,4] = cov[:,4,4:]

    # Mean vectors
    mean = np.dstack([ra,dec,dist,pmra,pmdec,rv])[0]
    
    vxvv_resample = np.empty((n_samples,6))
    for i in range(n_samples):
        try:
            vxvv_resample[i] = np.random.multivariate_normal(mean[i], cov[i], 1)
        except ValueError:
            print(mean[o_mask][i])
            print(cov[o_mask][i]) 
        ##te
    ###i  
    
    os_sample = orbit.Orbit(vxvv=vxvv_resample, radec=True, ro=ro, vo=vo, zo=zo)
    return os_sample
#def