# ----------------------------------------------------------------------------
#
# TITLE - potential.py
# PROJECT - mw-dfs
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Functions for dealing with potentials
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb, copy

## Astropy
from astropy import units as apu

## galpy
from galpy import orbit
from galpy import potential
from galpy.potential import DehnenSmoothWrapperPotential
from galpy import actionAngle as aA
from galpy import df
from galpy.util import conversion

# ----------------------------------------------------------------------------

def make_interpolated_mwpot(mwpot='MWPotential2014',rmin=1/8.,rmax=800.,
    ngrid=201,ro=8.,vo=220.,match_type='mass'):
    '''make_interpolated_mwpot:
    
    Make an interpolated version of MW Potential using 
    potential.interpSphericalPotential. Can either be made to match the 
    radial mass profile of the MW or the radial force profile. 
    
    Args:
        mwpot (string or galpy.potential.Potential) - Options are 'MWPotential2014'
        rmin (float) - Minimum radial position used to make the interpolated 
            grid
        rmax (float) - Maximum radial position used to make the interpolated
            grid
        ngrid (int) - Number of points in the radial interpolation grid
        ro (float) - galpy radial scale
        vo (float) - galpy velocity scale
        match_type (string) - Match the radial mass profile ('mass') or the 
            in-plane radial force profile ('force') ['mass'] 
    
    Returns:
        interpot (potential.interpSphericalPotential) - interpolated 
            MW Potential
    '''
    if isinstance(rmin,apu.quantity.Quantity):
        rmin = rmin.value/ro
    if isinstance(rmax,apu.quantity.Quantity):
        rmax = rmax.value/ro

    if isinstance(mwpot,potential.Potential) or isinstance(mwpot[0],potential.Potential):
        mwp = copy.deepcopy(mwpot)
    else:
        assert isinstance(mwpot,str), 'If not potential.Potential, mwpot must be str'
        if mwpot == 'MWPotential2014':
            mwp = potential.MWPotential2014
        elif mwpot == 'McMillan17':
            mwp = potential.mwpotentials.McMillan17
        ##ie
    ##ie
    potential.turn_physical_off(mwp)
    
    rgrid = np.geomspace(rmin,rmax,ngrid)
    assert match_type in ['mass','force'], 'match_type must be "mass" or "force"'
    if match_type == 'mass':
        rforce_fn = lambda r: -potential.mass(mwp,r,use_physical=False)/r**2
    elif match_type == 'force':
        rforce_fn = lambda r: potential.evaluaterforce(mwp,r,0,use_physical=False)
    ##ie
    interpot = potential.interpSphericalPotential(rforce=rforce_fn, 
        rgrid=rgrid, Phi0=potential.evaluatePotentials(mwp,rmin,0,use_physical=False), 
        ro=ro, vo=vo)
    return interpot
#defr

# ----------------------------------------------------------------------------

# def project_potential(ro,vo):
#     '''project_potential:
    
#     Return the potential which is used for this project. It's based on
#     MWPotential2014 but uses a Hernquist profile for the halo
    
#     Args:
#         ro (float) - galpy radial scale
#         vo (float) - galpy velocity scale
#     '''
#     bulge_pot,disk_pot,_ = potential.MWPotential2014
#     Mh = 1.2*10**12*apu.M_sun
#     ah = 42*apu.kpc
#     halo_pot = potential.HernquistPotential(amp=2*Mh,a=ah,ro=ro,vo=vo)
#     mwpot = [bulge_pot,disk_pot,halo_pot]
#     potential.turn_physical_on(mwpot,ro=ro,vo=vo)
#     return mwpot
# #def

# ----------------------------------------------------------------------------

# def powerlaw_expcut_density_profile(alpha,rc,ro,vo):
#     '''powerlaw_expcut_density_profile:
    
#     Make the powerlaw w/ exponentially truncated density profile for the 
#     paper DF
    
#     Args:
#         alpha (float) - Power law slope (should be positive)
#         rc (float, can be astropy quantity) - Exponential cutoff radius, 
#             can be 
#         ro (float) - galpy radial scale
#         vo (float) - galpy velocity scale
    
#     Returns:
#         denspot (galpy.Potential.Potential) - Density potential
#     '''
#     if not isinstance(rc,apu.quantity.Quantity):
#         rc *= apu.kpc
#     denspot = potential.PowerSphericalPotentialwCutoff(amp=1.,r1=1.,alpha=alpha,
#                                                        rc=rc,ro=ro,vo=vo)
#     return denspot
# #def

# ----------------------------------------------------------------------------

# def adiabatic_Hernquist_mwpot_disk_grow(Mh, ah, Mh_decay, tform, tsteady, ro, vo, zo):
#     '''adiabatic_Hernquist_mwpot_disk_gro:
    
#     Create a realization of MWPotential2014 with Hernquist halo that 
#     grows the disk and bulge component adiabatically. Useful for 
#     relaxing sampled halo orbits
    
#     Args:
#         Mh (float) - The mass of the final Hernquist halo in solar masses
#         ah (float) - The scale radius of the Hernquist halo in kpc
#         Mh_decay (float) - The mass of the halo that will decay while the 
#             disk and bulge are growing
#         tform (float) - Time at which adiabatic transformation begins
#         tsteady (float) - Duration of adiabatic transformation
#     '''
    
#     # Times in galpy units
#     tform = tform/conversion.time_in_Gyr(vo,ro)
#     tsteady = tsteady/conversion.time_in_Gyr(vo,ro)
    
#     # Fetch MWPotential2014
#     bulge, disk, _ = potential.MWPotential2014
#     potential.turn_physical_on(bulge,ro=ro,vo=vo)
#     potential.turn_physical_on(disk,ro=ro,vo=vo)
    
#     # The Hernquist halo
#     hernquist_halo = potential.HernquistPotential(amp=2*Mh, a=ah, ro=ro, vo=vo)
#     potential.turn_physical_on(hernquist_halo,ro=ro,vo=vo)
    
#     # Create a deficit halo that holds the mass of the disk and bulge while 
#     # they are growing
#     deficit_halo = potential.HernquistPotential(amp=2*Mh_decay*apu.M_sun, 
#                                                 a=ah, ro=ro, vo=vo)
#     potential.turn_physical_on(deficit_halo, ro=ro, vo=vo)
#     potential.turn_physical_off(deficit_halo)
    
#     # Now wrap all the potentials in Dehnen smooth wrappers
#     mwpot = [bulge, disk]
#     potential.turn_physical_off(mwpot)
#     potential.turn_physical_off(hernquist_halo)
#     mwpot_grow = DehnenSmoothWrapperPotential(pot=mwpot, tform=tform, 
#                                               tsteady=tsteady)
#     halo_decay = DehnenSmoothWrapperPotential(pot=deficit_halo, tform=tform, 
#                                               tsteady=tsteady, decay=True)
    
#     return potential.flatten([hernquist_halo,mwpot_grow,halo_decay])
# #def

# ----------------------------------------------------------------------------

# def adiabatic_mwpot_disk_grow(M_decay, tform, tsteady, ro, vo, zo):
#     '''adiabatic_mwpot_disk_gro:
    
#     Create a realization of MWPotential2014 that grows the disk and bulge 
#     component adiabatically. Useful for relaxing sampled halo orbits
    
#     Args:
#         Mh_decay (float) - The mass of the halo that will decay while the 
#             disk and bulge are growing in solar masses
#         tform (float) - Time at which adiabatic transformation begins
#         tsteady (float) - Duration of adiabatic transformation
#     '''
#     print('Warning! This function now returns a MWPotential with NFW halo instead of a Hernquist halo')
#     # Times in galpy units
#     tform = tform/conversion.time_in_Gyr(vo,ro)
#     tsteady = tsteady/conversion.time_in_Gyr(vo,ro)
    
#     # M_decay to unitless
#     if isinstance(M_decay,apu.quantity.Quantity): M_decay = M_decay.to(apu.M_sun).value
    
#     # Fetch MWPotential2014
#     bulge, disk, halo = potential.MWPotential2014
#     potential.turn_physical_on(bulge,ro=ro,vo=vo)
#     potential.turn_physical_on(disk,ro=ro,vo=vo)
#     potential.turn_physical_on(halo,ro=ro,vo=vo)
    
#     # Create a deficit halo that holds the mass of the disk and bulge while 
#     # they are growing
#     deficit_halo = potential.NFWPotential(mvir=M_decay/1e12, conc=halo.conc(), 
#                                           ro=ro, vo=vo)
#     potential.turn_physical_on(deficit_halo, ro=ro, vo=vo)
#     potential.turn_physical_off(deficit_halo)
    
#     # Now wrap all the potentials in Dehnen smooth wrappers
#     diskbulge = [bulge, disk]
#     potential.turn_physical_off(diskbulge)
#     potential.turn_physical_off(halo)
#     diskbulge_grow = DehnenSmoothWrapperPotential(pot=diskbulge, tform=tform, 
#                                               tsteady=tsteady)
#     halo_decay = DehnenSmoothWrapperPotential(pot=deficit_halo, tform=tform, 
#                                               tsteady=tsteady, decay=True)
    
#     return potential.flatten([halo,diskbulge_grow,halo_decay])
# #def

# ----------------------------------------------------------------------------

# def beta_in_radial_bins(os, rbin):
#     '''beta_in_radial_bins:
    
#     Calculate beta in radial bins defined by a set of bin edges
    
#     Args:
#         os (galpy.orbit.Orbit) - Orbits
#         rbin (array) - List of bin edges
#     '''
#     betas = np.zeros(len(rbin)-1)
#     vr = os.vr(use_physical=True).value
#     vphi = os.vT(use_physical=True).value
#     vtheta = os.vtheta(use_physical=True).value
#     for i in range(len(betas)):
#         where_in_bin = np.where(np.logical_and(os.r(use_physical=True).value>rbin[i], 
#                                                os.r(use_physical=True).value<rbin[i+1]
#                                                ))[0]
#         if len(where_in_bin) < 2:
#             betas[i] = np.NaN
#             continue
#         ##fi
#         svr = np.std( vr[where_in_bin] )
#         svphi = np.std( vphi[where_in_bin] )
#         svtheta = np.std( vtheta[where_in_bin] )
#         betas[i] = 1 - ( (svphi**2.+svtheta**2.)/(2*svr**2.) )
#     ###i
    
#     return betas
        

# ----------------------------------------------------------------------------

def get_ELz_boundary(pot,r):
    '''get_ELz_boundary:
    
    Get the energy - angular momentum boundary at a given radius for a given 
    potential. Defined as where the angular momentum is that of the escape 
    velocity. Will return the curve between the limits where E=0 
    
    Args:
        pot (galpy.potential.Potential) - Potential
        r (float) - radius at which to calculate the ELz boundary
        
    Returns:
        Lz_bound (array) - Lz of the boundary
        E_bound (array) - Energy of the boundary
    '''
    # First get the E=0 limits
    vesc = pot.vesc(r)
    Lz_lim = [-r*vesc,r*vesc]
    
    # Now calculate the curve
    Lz_bound = np.linspace(Lz_lim[0].value, Lz_lim[1].value, num=100)*apu.kpc*apu.km/apu.s
    E_bound = np.power(Lz_bound,2)/(2*r**2) + potential.evaluatePotentials(pot,r,0)
    
    return Lz_bound, E_bound
#def

# ----------------------------------------------------------------------------