# ----------------------------------------------------------------------------
#
# TITLE - plot.py
# AUTHOR - James Lane
# PROJECT - mw-dfs
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Plotting routines for the Milky Way DF project
'''
__author__ = "James Lane"

### Imports
import numpy as np
import numbers
import copy
import dill as pickle
from astropy import units as apu
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import patches
from galpy import potential
import scipy.stats
from scipy.special import erf
import warnings

# ----------------------------------------------------------------------------

def parse_mixture(orbs, mixture_arr, seed=None, absolute=True, return_inds=False):
    '''parse_mixture:
    
    Parses an array that denotes fractional numbers of orbits to be returned 
    as well as a list of orbit.Orbit objects and returns a list of orbit.Orbit 
    objects corresponding to the mixture
    
    In practice, the idea is to take many orbit.Orbit objects that have the 
    same number of orbits in them, then generate a new list of orbits which 
    has the same number of orbits as each of the original orbit.Orbit objects.
    
    if you had a list of orbit.Orbit objects, each of which 
    have 1000 orbits in them:
    
    [orb1,orb2]
    
    and a mixture_arr that looked like:
    
    [0.4,0.6]
    
    Then 400 orbits are taken from orb1 and 600 from orb2 and a list of 
    those orbitss is returned:
    
    [orb1_mixture,orb2_mixture]
    
    Note that if absolute=False, actual numbers are normalized, so a mixture 
    arr which is [4,6] Does the same thing as one which is [0.4,0.6]
    
    Args:
        orbs (list) - List of galpy.orbit.Orbit objects
        mixture_arr (list) - List of numbers corresponding to the fractional 
            amount of each element of orbs that should be returned
        seed (int) - Seed to use for sampling
        absolute (bool) - The fractions given in mixture_arr are absolute, 
            not relative. If True no element of absolute may be greater than 
            1.0 [True]
    
    Returns:
        orbs_mixture (list) - List of galpy.orbit.Orbit objects with fractional 
            amounts corresponding to mixture_arr
    '''
    if not absolute:
        warnings.warn('Warning, absolute=False. Not using absolute mixture array fractions!')
    assert len(mixture_arr) == len(orbs)
    assert isinstance(mixture_arr,np.ndarray)
    
    norm = np.sum(mixture_arr)
    
    orbs_mixture = []
    sample_inds_mixture = []
    
    for i in range(len(mixture_arr)):
        
        if mixture_arr[i] == 0:
            continue
        elif mixture_arr[i] == norm and not absolute:
            return [orbs[i],]
        else:
            if absolute:
                orb_frac = mixture_arr[i]
                assert orb_frac <= 1., 'absolute fractions cannot exceed 1'
                if orb_frac == 1.:
                    orbs_mixture.append(orbs[i])
                    continue
            else:
                orb_frac = mixture_arr[i]/norm
            try:
                n_orbs = orbs[i].shape[1]
                _isNumpy = True
            except IndexError:
                _isNumpy = False
            if _isNumpy:
                n_samples = int(n_orbs*orb_frac)
                if seed is not None:
                    np.random.seed(seed)
                sample_inds = np.random.choice(n_orbs,n_samples,replace=False)
                orbs_mixture.append(orbs[i][:,sample_inds])
                sample_inds_mixture.append(sample_inds)
            else:
                n_orbs = len(orbs[i])
                n_samples = int(n_orbs*orb_frac)
                if seed is not None:
                    np.random.seed(seed)
                sample_inds = np.random.choice(n_orbs,n_samples,replace=False)
                orbs_mixture.append(orbs[i][sample_inds])
                sample_inds_mixture.append(sample_inds)
            ##te
        ##ie
    ###i
    if return_inds:
        return orbs_mixture, sample_inds_mixture
    else:
        return orbs_mixture
    ##ie
#def

# def plot_velocity_ellipsoids_mixture(orbs, mixture_arr, mixture_text=None, coords='cylindrical', vlim=300, velocity_uncertainty=None):
#     '''plot_velocity_ellipsoids_mixture:
    
#     Plots the velocity ellipsoids for a list of orbit.Orbit objects. For each plot allow for a mixture from many orbits
    
#     Args:
#         orbs (array) - Array of different orbit.Orbit objects
#         mixture_arr (list) - list of arrays, each of size len(orbs), which 
#             denote the fraction of samples from each orbit.Orbit object in 
#             orbs that will be plotted
#         mixture_text (array) - Array of text for each mixture instance
#         coords (str) - Velocity coordinates to use: 'rectangular', 
#             'cylindrical', 'spherical'
#         vlim (float) - Velocity limits for the plots
#         velocity_uncertainty (float) - Standard deviation of gaussian errors 
#             artificially applied to velocities [0]
#     '''
#     n_orbs = len(orbs)
#     n_mix = len(mixture_arr)
    
#     fig = plt.figure( figsize=(15,5*n_mix) )
#     axs = fig.subplots(nrows=n_mix,ncols=3).reshape(n_mix*3)
#     assert len(mixture_text) == n_mix
#     if mixture_text == None:
#         mixture_text = np.empty(n_mix,dtype='string')
#         mixture_text[:] = ''
#     ##fi

#     for i in range(n_mix):
        
#         orbs_mixture = parse_mixture(orbs, mixture_arr[i])
        
#         for j in range(len(orbs_mixture)):
        
#             if coords == 'rectangular':
#                 v1,v2,v3 = orbs_mixture[j].vx(),orbs_mixture[j].vy(),orbs_mixture[j].vz()
#                 v1_label,v2_label,v3_label = r'$v_{x}$', r'$v_{y}$', r'$v_{z}$' 
#             ##fi
#             if coords == 'cylindrical':
#                 v1,v2,v3 = orbs_mixture[j].vR(),orbs_mixture[j].vT(),orbs_mixture[j].vz()
#                 v1_label,v2_label,v3_label = r'$v_{R}$', r'$v_{T}$', r'$v_{z}$' 
#             ##fi
#             if coords == 'spherical':
#                 orbs_mixture_r = np.sqrt( np.square(orbs_mixture[j].R) + np.square(orbs_mixture[j].z) )
#                 v1 = (orbs_mixture[j].vR()*orbs_mixture[j].R() + orbs_mixture[j].vz()*orbs_mixture[j].z()) / orbs_mixture_r
#                 v2 = orbs_mixture[j].vT()
#                 v3 = (orbs_mixture[j].R()*orbs_mixture[j].vz() - orbs_mixture[j].vR()*orbs_mixture[j].z()) / orbs_mixture_r
#                 v1_label,v2_label,v3_label = r'$v_{r}$', r'$v_{T}$', r'$v_{\theta}$' 
#             ##fi
            
#             # Apply Gaussian uncertainties to 
#             if velocity_uncertainty is not None:
#                 if isinstance(velocity_uncertainty,float) or isinstance(velocity_uncertainty,int):
#                     velocity_uncertainty = np.ones(3)*velocity_uncertainty
#                 n_os = len(orbs_mixture[j])
#                 v1 += np.random.normal(loc=0, scale=velocity_uncertainty[0], size=n_os)*apu.km/apu.s
#                 v2 += np.random.normal(loc=0, scale=velocity_uncertainty[1], size=n_os)*apu.km/apu.s
#                 v3 += np.random.normal(loc=0, scale=velocity_uncertainty[2], size=n_os)*apu.km/apu.s
#             ##fi
            
#             axs[3*i].scatter( v1, v2, s=0.5, alpha=0.25, color='Black' )
#             axs[3*i+1].scatter( v1, v3, s=0.5, alpha=0.25, color='Black' )
#             axs[3*i+2].scatter( v2, v3, s=0.5, alpha=0.25, color='Black' )

#             axs[3*i].set_xlabel(v1_label, fontsize=16)
#             axs[3*i].set_ylabel(v2_label, fontsize=16)
#             axs[3*i+1].set_xlabel(v1_label, fontsize=16)
#             axs[3*i+1].set_ylabel(v3_label, fontsize=16)
#             axs[3*i+2].set_xlabel(v2_label, fontsize=16)
#             axs[3*i+2].set_ylabel(v3_label, fontsize=16)

#             axs[3*i].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#             axs[3*i+1].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#             axs[3*i+2].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#             if velocity_uncertainty is not None:
#                 axs[3*i].annotate(r'$\sigma = $ '+str(round(velocity_uncertainty[0])), xy=(0.025,0.05), xycoords='axes fraction', fontsize=11)
#                 axs[3*i+1].annotate(r'$\sigma = $ '+str(round(velocity_uncertainty[1])), xy=(0.025,0.05), xycoords='axes fraction', fontsize=11)
#                 axs[3*i+2].annotate(r'$\sigma = $ '+str(round(velocity_uncertainty[2])), xy=(0.025,0.05), xycoords='axes fraction', fontsize=11)
#             ##fi
            
#             axs[3*i].set_xlim(-vlim,vlim)
#             axs[3*i].set_ylim(-vlim,vlim)
#             axs[3*i+1].set_xlim(-vlim,vlim)
#             axs[3*i+1].set_ylim(-vlim,vlim)
#             axs[3*i+2].set_xlim(-vlim,vlim)
#             axs[3*i+2].set_ylim(-vlim,vlim)
#         ###o
#     ###i
#     fig.subplots_adjust(wspace=0.3)
#     return fig, axs
# #def

# def plot_kinematics_mixture(eELzs, mixture_arr, mixture_text=None):
#     '''plot_kinematics_mixture:
    
#     Plots the conserved kinematic quantities e, Lz, and E
    
#     Args:
#         eELzs (array) - list of arrays (each shape 3xN) containing e, E, and Lz
#         pot (galpy.potential.Potential) - Potential for deriving kinematic quantities
#         mixture_arr (list) - list of arrays, each of size len(orbs), which denote the fraction
#             of samples from each orbit.Orbit object in orbs that will be plotted
#         mixture_text (array) - Array of text for each mixture instance
#         orbs_text (array) - Array of text for each orbit.Orbit object
#     '''
#     n_eELzs = len(eELzs)
#     n_mix = len(mixture_arr)
    
#     fig = plt.figure( figsize=(15,5*n_mix) )
#     axs = fig.subplots(nrows=n_mix,ncols=3).reshape(n_mix*3)
#     assert len(mixture_text) == n_mix
#     if mixture_text == None:
#         mixture_text = np.empty(n_mix,dtype='string')
#         mixture_text[:] = ''
#     ##fi

#     for i in range(n_mix):
                    
#         eELzs_mixture = parse_mixture(eELzs, mixture_arr[i])
        
#         for j in range(len(eELzs_mixture)):
            
#             e,E,Lz = eELzs_mixture[j]
        
#             axs[3*i].scatter( Lz, e, s=0.5, alpha=0.25, color='Black' )
#             axs[3*i+1].scatter( Lz, E, s=0.5, alpha=0.25, color='Black' )
#             axs[3*i+2].scatter( e, E, s=0.5, alpha=0.25, color='Black' )
#         ###j
        
#         axs[3*i].set_xlabel(r'$L_{z}$', fontsize=16)
#         axs[3*i].set_ylabel(r'e', fontsize=16)
#         axs[3*i+1].set_xlabel(r'$L_{z}$', fontsize=16)
#         axs[3*i+1].set_ylabel(r'$E$', fontsize=16)
#         axs[3*i+2].set_xlabel(r'e', fontsize=16)
#         axs[3*i+2].set_ylabel(r'$E$', fontsize=16)
        
#         axs[3*i].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#         axs[3*i+1].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#         axs[3*i+2].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)

#         e_lim = [-0.05,1.05]
#         Etot_lim = [-160000,-80000]
#         Lz_lim = [-3000,3000]
#         #axs[3*i].set_xlim(Lz_lim[0],Lz_lim[1])
#         #axs[3*i].set_ylim(e_lim[0],e_lim[1])
#         #axs[3*i+1].set_xlim(Lz_lim[0],Lz_lim[1])
#         #axs[3*i+1].set_ylim(Etot_lim[0],Etot_lim[1])
#         #axs[3*i+2].set_xlim(e_lim[0],e_lim[1])
#         #axs[3*i+2].set_ylim(Etot_lim[0],Etot_lim[1])
        
#     ###i
#     fig.subplots_adjust(wspace=0.3)
#     return fig, axs
# #def

# def plot_actions_mixture(accs, mixture_arr, mixture_text=None):
#     '''plot_actions_mixture:
    
#     Plots the actions: jR, jphi, jz
    
#     Args:
#         accs (array) - Array of different action arrays (array of 3-element lists of actions [jR,jphi,jz])
#         mixture_arr (list) - list of arrays, each of size len(orbs), which denote the fraction
#             of samples from each orbit.Orbit object in orbs that will be plotted
#         mixture_text (array) - Array of text for each mixture instance
#     '''
#     n_orbs = len(accs)
#     n_mix = len(mixture_arr)
    
#     fig = plt.figure( figsize=(15,5*n_mix) )
#     axs = fig.subplots(nrows=n_mix,ncols=3).reshape(n_mix*3)
#     assert len(mixture_text) == n_mix
#     if mixture_text == None:
#         mixture_text = np.empty(n_mix,dtype='string')
#         mixture_text[:] = ''
#     ##fi
    
#     jR_text, Lz_text, jz_text = r'$J_{R}$', r'$L_{z}$', r'$J_{z}$'

#     for i in range(n_mix):
        
#         accs_mixture = parse_mixture(accs, mixture_arr[i])
        
#         for j in range(len(accs_mixture)):
#             jR,Lz,jz = accs_mixture[j]
#             jR[jR > 1e4] = np.nan
#             Lz[jR > 1e4] = np.nan
#             jz[jR > 1e4] = np.nan
#             axs[3*i].scatter( jR, Lz, s=0.5, alpha=0.25, color='Black' )
#             axs[3*i+1].scatter( jR, jz, s=0.5, alpha=0.25, color='Black' )
#             axs[3*i+2].scatter( Lz, jz, s=0.5, alpha=0.25, color='Black' )
#         ###j
        
#         axs[3*i].set_xlabel(jR_text, fontsize=16)
#         axs[3*i].set_ylabel(Lz_text, fontsize=16)
#         axs[3*i+1].set_xlabel(jR_text, fontsize=16)
#         axs[3*i+1].set_ylabel(jz_text, fontsize=16)
#         axs[3*i+2].set_xlabel(Lz_text, fontsize=16)
#         axs[3*i+2].set_ylabel(jz_text, fontsize=16)

#         axs[3*i].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#         axs[3*i+1].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#         axs[3*i+2].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)

#         jR_lim = [0,5000]
#         Lz_lim = [-3000,3000]
#         jz_lim = [0,2500]
#     ###i
        
#     fig.subplots_adjust(wspace=0.3)
#     return fig, axs
# #def

# def plot_action_map_mixture(accs, mixture_arr, mixture_text=None):
#     '''plot_action_map_mixture:
    
#     Plot the wierd action map people so love
    
#     Args:
#         accs (array) - Array of different action arrays (array of 3-element lists of actions [jR,jphi,jz])
#         mixture_arr (list) - list of arrays, each of size len(orbs), which denote the fraction
#             of samples from each orbit.Orbit object in orbs that will be plotted
#         mixture_text (array) - Array of text for each mixture instance
#     '''
#     n_orbs = len(accs)
#     n_mix = len(mixture_arr)
#     n_rows = int(np.ceil(n_mix/3))
#     remove_axes_beyond = n_rows*3 - n_mix
    
#     fig = plt.figure( figsize=(15,5*n_rows) )
#     axs = fig.subplots(nrows=n_rows,ncols=3).flatten()
#     for i in range(remove_axes_beyond):
#         axs[-i-1].set_axis_off()
#     ###i

#     assert len(mixture_text) == n_mix
#     if mixture_text == None:
#         mixture_text = np.empty(n_mix,dtype='string')
#         mixture_text[:] = ''
#     ##fi

#     for i in range(n_mix):
        
#         accs_mixture = parse_mixture(accs, mixture_arr[i])
        
#         for j in range(len(accs_mixture)):
#             jR,Lz,jz = accs_mixture[j]
#             Jtot = np.abs(jR) + np.abs(Lz) + np.abs(jz)
#             Jz_JR_norm = (jz-jR) / Jtot
#             Jphi_norm = Lz / Jtot
#             axs[i].scatter(Jphi_norm, Jz_JR_norm, s=0.5, alpha=0.25, color='Black')
#         ###j
        
#         axs[i].plot([-1,0], [0,-1], linestyle='dashed', color='Black', linewidth=0.5)
#         axs[i].plot([0,1],  [-1,0], linestyle='dashed', color='Black', linewidth=0.5)
#         axs[i].plot([1,0],  [0,1],  linestyle='dashed', color='Black', linewidth=0.5)
#         axs[i].plot([0,-1], [1,0],  linestyle='dashed', color='Black', linewidth=0.5)
#         axs[i].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#         axs[i].set_xlabel(r'$J_{\phi}/J_{\rm tot}$', fontsize=16)
#         axs[i].set_ylabel(r'$(J_{\rm z}-J_{\rm R})/J_{\rm tot}$', fontsize=16)
        
#         axs[i].set_xlim(-1,1)
#         axs[i].set_ylim(-1,1)
#     ###i
#     fig.subplots_adjust(wspace=0.3)
#     return fig, axs
# #def

# def plot_velocity_ellipsoids_purity(orbs1, orbs2, mixture_arr1, mixture_arr2, coords='cylindrical', vlim=300, velocity_uncertainty=None):
#     '''plot_velocity_ellipsoids_purity:
    
#     Plots the purity of two samples of orbit.Orbit objects. For each plot allow for a mixture from many orbits. Purity of (1) will be 
#     calculated with respect to (1 and 2)
    
#     Args:
#         orbs(1/2) (array) - Array of different orbit.Orbit objects
#         mixture_arr(1/2) (list) - list of arrays, each of size len(orbs), which denote the fraction
#             of samples from each orbit.Orbit object in orbs that will primarily shown
#         coords (str) - Velocity coordinates to use: 'rectangular', 'cylindrical', 'spherical'
#         vlim (float) - Velocity limits for the plots
#         velocity_uncertainty (float) - Standard deviation of gaussian errors artificially applied to velocities [0]
#     '''
#     assert len(mixture_arr1) == len(mixture_arr2), 'Mixture array lengths must be the same'
#     n_mix = len(mixture_arr1)
    
#     assert coords == 'cylindrical', 'only cylindrical coordinates currently supported'
#     v1_label,v2_label,v3_label = r'$v_{R}$', r'$v_{T}$', r'$v_{z}$' 
    
#     fig = plt.figure( figsize=(18,5*n_mix) )
#     axs = fig.subplots(nrows=n_mix,ncols=3).reshape(n_mix*3)
    
#     h1_out = []
#     h2_out = []

#     for i in range(n_mix):
        
#         orbs_mixture1 = parse_mixture(orbs1, mixture_arr1[i])
#         orbs_mixture2 = parse_mixture(orbs2, mixture_arr2[i])
        
#         v1_1 = np.array([])
#         v2_1 = np.array([])
#         v3_1 = np.array([])
#         for j in range(len(orbs_mixture1)):
#             v1_temp,v2_temp,v3_temp = orbs_mixture1[j].vR(),orbs_mixture1[j].vT(),orbs_mixture1[j].vz() 
#             v1_1 = np.append(v1_1,v1_temp.value)
#             v2_1 = np.append(v2_1,v2_temp.value)
#             v3_1 = np.append(v3_1,v3_temp.value)
        
#         v1_2 = np.array([])
#         v2_2 = np.array([])
#         v3_2 = np.array([])
#         for j in range(len(orbs_mixture2)):    
#             v1_temp,v2_temp,v3_temp = orbs_mixture2[j].vR(),orbs_mixture2[j].vT(),orbs_mixture2[j].vz()
#             v1_2 = np.append(v1_2,v1_temp.value)
#             v2_2 = np.append(v2_2,v2_temp.value)
#             v3_2 = np.append(v3_2,v3_temp.value)
            
#         # Apply Gaussian uncertainties to 
#         if velocity_uncertainty is not None:
#             if isinstance(velocity_uncertainty,float) or isinstance(velocity_uncertainty,int):
#                 velocity_uncertainty = np.ones(3)*velocity_uncertainty
#             n_os1 = len(v1_1)
#             v1_1 += np.random.normal(loc=0, scale=velocity_uncertainty[0], size=n_os1)
#             v2_1 += np.random.normal(loc=0, scale=velocity_uncertainty[1], size=n_os1)
#             v3_1 += np.random.normal(loc=0, scale=velocity_uncertainty[2], size=n_os1)
#             n_os2 = len(v1_2)
#             v1_2 += np.random.normal(loc=0, scale=velocity_uncertainty[0], size=n_os2)
#             v2_2 += np.random.normal(loc=0, scale=velocity_uncertainty[1], size=n_os2)
#             v3_2 += np.random.normal(loc=0, scale=velocity_uncertainty[2], size=n_os2)
#         ##fi
                
#         h1_1 = np.histogram2d( v1_1, v2_1, bins=50, range=([-vlim,vlim],[-vlim,vlim]) )[0]
#         h2_1 = np.histogram2d( v1_1, v3_1, bins=50, range=([-vlim,vlim],[-vlim,vlim]) )[0]
#         h3_1 = np.histogram2d( v2_1, v3_1, bins=50, range=([-vlim,vlim],[-vlim,vlim]) )[0]
#         h1_2 = np.histogram2d( v1_2, v2_2, bins=50, range=([-vlim,vlim],[-vlim,vlim]) )[0]
#         h2_2 = np.histogram2d( v1_2, v3_2, bins=50, range=([-vlim,vlim],[-vlim,vlim]) )[0]
#         h3_2 = np.histogram2d( v2_2, v3_2, bins=50, range=([-vlim,vlim],[-vlim,vlim]) )[0]
        
#         p1 = np.divide( h1_1 , ( h1_1 + h1_2 ) )
#         p2 = np.divide( h2_1 , ( h2_1 + h2_2 ) )
#         p3 = np.divide( h3_1 , ( h3_1 + h3_2 ) )

#         i1 = axs[3*i].imshow( np.rot90(p1), cmap='rainbow', vmin=0, vmax=1, extent=(-vlim,vlim,-vlim,vlim))
#         i2 = axs[3*i+1].imshow( np.rot90(p2), cmap='rainbow', vmin=0, vmax=1, extent=(-vlim,vlim,-vlim,vlim))
#         i3 = axs[3*i+2].imshow( np.rot90(p3), cmap='rainbow', vmin=0, vmax=1, extent=(-vlim,vlim,-vlim,vlim))

#         cm1 = fig.colorbar(i1, ax=axs[3*i])
#         cm2 = fig.colorbar(i2, ax=axs[3*i+1])
#         cm3 = fig.colorbar(i3, ax=axs[3*i+2])
#         cm1.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')
#         cm2.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')
#         cm3.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')

#         axs[3*i].set_xlabel(v1_label, fontsize=16)
#         axs[3*i].set_ylabel(v2_label, fontsize=16)
#         axs[3*i+1].set_xlabel(v1_label, fontsize=16)
#         axs[3*i+1].set_ylabel(v3_label, fontsize=16)
#         axs[3*i+2].set_xlabel(v2_label, fontsize=16)
#         axs[3*i+2].set_ylabel(v3_label, fontsize=16)

#         #axs[3*i].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#         #axs[3*i+1].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#         #axs[3*i+2].annotate(mixture_text[i], xy=(0.025,0.95), xycoords='axes fraction', fontsize=11)
#         if velocity_uncertainty is not None:
#             axs[3*i].annotate(r'$\sigma = $ '+str(round(velocity_uncertainty[0])), xy=(0.025,0.05), xycoords='axes fraction', fontsize=11)
#             axs[3*i+1].annotate(r'$\sigma = $ '+str(round(velocity_uncertainty[1])), xy=(0.025,0.05), xycoords='axes fraction', fontsize=11)
#             axs[3*i+2].annotate(r'$\sigma = $ '+str(round(velocity_uncertainty[2])), xy=(0.025,0.05), xycoords='axes fraction', fontsize=11)
#         ##fi

#         axs[3*i].set_xlim(-vlim,vlim)
#         axs[3*i].set_ylim(-vlim,vlim)
#         axs[3*i+1].set_xlim(-vlim,vlim)
#         axs[3*i+1].set_ylim(-vlim,vlim)
#         axs[3*i+2].set_xlim(-vlim,vlim)
#         axs[3*i+2].set_ylim(-vlim,vlim)
        
#         h1_out.append([h1_1,h2_1,h3_1])
#         h2_out.append([h1_2,h2_2,h3_2])
#     ###i
#     fig.subplots_adjust(wspace=0.3)
#     return fig, axs, h1_out, h2_out
# #def

# def plot_kinematics_mixture_purity(eELzs1, eELzs2, mixture_arr1, mixture_arr2, mixture_text=None):
#     '''plot_kinematics_mixture_purity:
    
#     Plots the conserved kinematic quantities e, Lz, and E. Purity of (1) will be calculated with respect to (1 and 2)
    
#     Args:
#         eELzs(1/2) (array) - list of arrays (each shape 3xN) containing e, E, and Lz
#         mixture_arr(1/2) (list) - list of arrays, each of size len(orbs), which denote the fraction
#             of samples from each orbit.Orbit object in orbs that will primarily shown (i.e. the purity of this)
#         pot (galpy.potential.Potential) - Potential for deriving kinematic quantities
#         mixture_arr (list) - list of arrays, each of size len(orbs), which denote the fraction
#             of samples from each orbit.Orbit object in orbs that will be plotted
#         mixture_text (array) - Array of text for each mixture instance
#         orbs_text (array) - Array of text for each orbit.Orbit object
#     '''
#     assert len(mixture_arr1) == len(mixture_arr2), 'Mixture array lengths must be the same'
#     n_mix = len(mixture_arr1)
    
#     fig = plt.figure( figsize=(18,5*n_mix) )
#     axs = fig.subplots(nrows=n_mix,ncols=3).reshape(n_mix*3)
# #     assert len(mixture_text) == n_mix
# #     if mixture_text == None:
# #         mixture_text = np.empty(n_mix,dtype='string')
# #         mixture_text[:] = ''
# #     ##fi
    
#     h1_out = []
#     h2_out = []
    
#     for i in range(n_mix):
                    
#         eELzs_mixture1 = parse_mixture(eELzs1, mixture_arr1[i])
#         eELzs_mixture2 = parse_mixture(eELzs2, mixture_arr2[i])
        
#         e_1 = np.array([])
#         E_1 = np.array([])
#         Lz_1 = np.array([])
#         for j in range(len(eELzs_mixture1)):
#             e_temp,E_temp,Lz_temp = eELzs_mixture1[j]
#             if isinstance(E_temp,apu.quantity.Quantity):
#                 E_temp = E_temp.value
#             if isinstance(Lz_temp,apu.quantity.Quantity):
#                 Lz_temp = Lz_temp.value
#             e_1 = np.append(e_1,e_temp)
#             E_1 = np.append(E_1,E_temp)
#             Lz_1 = np.append(Lz_1,Lz_temp)
#         ###j
          
#         e_2 = np.array([])
#         E_2 = np.array([])
#         Lz_2 = np.array([])
#         for j in range(len(eELzs_mixture2)):
#             e_temp,E_temp,Lz_temp = eELzs_mixture2[j]
#             if isinstance(E_temp,apu.quantity.Quantity):
#                 E_temp = E_temp.value
#             if isinstance(Lz_temp,apu.quantity.Quantity):
#                 Lz_temp = Lz_temp.value
#             e_2 = np.append(e_2,e_temp)
#             E_2 = np.append(E_2,E_temp)
#             Lz_2 = np.append(Lz_2,Lz_temp)
#         ###j
        
#         Elim = [-70000,0]
#         Lzlim = [-2500,2500]
#         elim = [0,1]
        
#         h1_1 = np.histogram2d( Lz_1, e_1, bins=50, range=([Lzlim[0],Lzlim[1]],[elim[0],elim[1]]) )[0]
#         h2_1 = np.histogram2d( Lz_1, E_1, bins=50, range=([Lzlim[0],Lzlim[1]],[Elim[0],Elim[1]]) )[0]
#         h3_1 = np.histogram2d( e_1, E_1, bins=50, range=([elim[0],elim[1]],[Elim[0],Elim[1]]) )[0]
#         h1_2 = np.histogram2d( Lz_2, e_2, bins=50, range=([Lzlim[0],Lzlim[1]],[elim[0],elim[1]]) )[0]
#         h2_2 = np.histogram2d( Lz_2, E_2, bins=50, range=([Lzlim[0],Lzlim[1]],[Elim[0],Elim[1]]) )[0]
#         h3_2 = np.histogram2d( e_2, E_2, bins=50, range=([elim[0],elim[1]],[Elim[0],Elim[1]]) )[0]
        
#         p1 = np.divide( h1_1 , ( h1_1 + h1_2 ) )
#         p2 = np.divide( h2_1 , ( h2_1 + h2_2 ) )
#         p3 = np.divide( h3_1 , ( h3_1 + h3_2 ) )
                
#         i1 = axs[3*i].imshow( np.rot90(p1), cmap='rainbow', vmin=0, vmax=1, aspect='auto', extent=(Lzlim[0],Lzlim[1],elim[0],elim[1]))
#         i2 = axs[3*i+1].imshow( np.rot90(p2), cmap='rainbow', vmin=0, vmax=1, aspect='auto', extent=(Lzlim[0],Lzlim[1],Elim[0],Elim[1]))
#         i3 = axs[3*i+2].imshow( np.rot90(p3), cmap='rainbow', vmin=0, vmax=1, aspect='auto', extent=(elim[0],elim[1],Elim[0],Elim[1]))

#         cm1 = fig.colorbar(i1, ax=axs[3*i])
#         cm2 = fig.colorbar(i2, ax=axs[3*i+1])
#         cm3 = fig.colorbar(i3, ax=axs[3*i+2])
#         cm1.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')
#         cm2.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')
#         cm3.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')

#         axs[3*i].set_xlabel(r'$L_{z}$', fontsize=16)
#         axs[3*i].set_ylabel(r'$e$', fontsize=16)
#         axs[3*i+1].set_xlabel(r'$L_{z}$', fontsize=16)
#         axs[3*i+1].set_ylabel(r'$E$', fontsize=16)
#         axs[3*i+2].set_xlabel(r'$e$', fontsize=16)
#         axs[3*i+2].set_ylabel(r'$E$', fontsize=16)

#         axs[3*i].set_xlim(Lzlim[0],Lzlim[1])
#         axs[3*i].set_ylim(elim[0],elim[1])
#         axs[3*i+1].set_xlim(Lzlim[0],Lzlim[1])
#         axs[3*i+1].set_ylim(Elim[0],Elim[1])
#         axs[3*i+2].set_xlim(elim[0],elim[1])
#         axs[3*i+2].set_ylim(Elim[0],Elim[1])
        
#         h1_out.append([h1_1,h2_1,h3_1])
#         h2_out.append([h1_2,h2_2,h3_2])
#     ###i
#     fig.subplots_adjust(wspace=0.3)
#     return fig, axs, h1_out, h2_out
# #def

# def plot_actions_purity(accs1, accs2, mixture_arr1, mixture_arr2):
#     '''plot_actions_purity:
    
#     Plots the purity of actions for sample 1 with sample 2 being contaminent
    
#     Args:
#         accs(1/2) (array) - Array of different action arrays (array of 3-element lists of actions [jR,jphi,jz])
#         mixture_arr(1/2) (list) - list of arrays, each of size len(orbs), which denote the fraction
#             of samples from each orbit.Orbit object in orbs that will be plotted
#     '''
#     assert len(mixture_arr1) == len(mixture_arr2), 'Mixture array lengths must be the same'
#     n_mix = len(mixture_arr1)
    
#     fig = plt.figure( figsize=(15,5*n_mix) )
#     axs = fig.subplots(nrows=n_mix,ncols=3).reshape(n_mix*3)
    
#     jR_text, Lz_text, jz_text = r'$J_{R}$', r'$L_{z}$', r'$J_{z}$'
    
#     h1_out = []
#     h2_out = []

#     for i in range(n_mix):
        
#         accs_mixture1 = parse_mixture(accs1, mixture_arr1[i])
#         accs_mixture2 = parse_mixture(accs2, mixture_arr2[i])
        
#         jR_1 = np.array([])
#         Lz_1 = np.array([])
#         jz_1 = np.array([])
#         for j in range(len(accs_mixture1)):
#             jR_temp,Lz_temp,jz_temp = accs_mixture1[j]
#             jR_1 = np.append(jR_1,jR_temp)
#             Lz_1 = np.append(Lz_1,Lz_temp)
#             jz_1 = np.append(jz_1,jz_temp)
#         ###j
        
#         jR_2 = np.array([])
#         Lz_2 = np.array([])
#         jz_2 = np.array([])
#         for j in range(len(accs_mixture2)):
#             jR_temp,Lz_temp,jz_temp = accs_mixture2[j]
#             jR_2 = np.append(jR_2,jR_temp)
#             Lz_2 = np.append(Lz_2,Lz_temp)
#             jz_2 = np.append(jz_2,jz_temp)
#         ###j
        
#         jRlim = [0,2000]
#         Lzlim = [-2500,2500]
#         jzlim = [0,2000]
        
#         h1_1 = np.histogram2d( jR_1, Lz_1, bins=50, range=([jRlim[0],jRlim[1]],[Lzlim[0],Lzlim[1]]) )[0]
#         h2_1 = np.histogram2d( jR_1, jz_1, bins=50, range=([jRlim[0],jRlim[1]],[jzlim[0],jzlim[1]]) )[0]
#         h3_1 = np.histogram2d( Lz_1, jz_1, bins=50, range=([Lzlim[0],Lzlim[1]],[jzlim[0],jzlim[1]]) )[0]
#         h1_2 = np.histogram2d( jR_2, Lz_2, bins=50, range=([jRlim[0],jRlim[1]],[Lzlim[0],Lzlim[1]]) )[0]
#         h2_2 = np.histogram2d( jR_2, jz_2, bins=50, range=([jRlim[0],jRlim[1]],[jzlim[0],jzlim[1]]) )[0]
#         h3_2 = np.histogram2d( Lz_2, jz_2, bins=50, range=([Lzlim[0],Lzlim[1]],[jzlim[0],jzlim[1]]) )[0]
        
#         p1 = np.divide( h1_1 , ( h1_1 + h1_2 ) )
#         p2 = np.divide( h2_1 , ( h2_1 + h2_2 ) )
#         p3 = np.divide( h3_1 , ( h3_1 + h3_2 ) )
                
#         i1 = axs[3*i].imshow( np.rot90(p1), cmap='rainbow', vmin=0, vmax=1, aspect='auto', extent=(jRlim[0],jRlim[1],Lzlim[0],Lzlim[1]))
#         i2 = axs[3*i+1].imshow( np.rot90(p2), cmap='rainbow', vmin=0, vmax=1, aspect='auto', extent=(jRlim[0],jRlim[1],jzlim[0],jzlim[1]))
#         i3 = axs[3*i+2].imshow( np.rot90(p3), cmap='rainbow', vmin=0, vmax=1, aspect='auto', extent=(Lzlim[0],Lzlim[1],jzlim[0],jzlim[1]))

#         cm1 = fig.colorbar(i1, ax=axs[3*i])
#         cm2 = fig.colorbar(i2, ax=axs[3*i+1])
#         cm3 = fig.colorbar(i3, ax=axs[3*i+2])
#         cm1.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')
#         cm2.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')
#         cm3.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')
        
        
#         axs[3*i].set_xlabel(jR_text, fontsize=16)
#         axs[3*i].set_ylabel(Lz_text, fontsize=16)
#         axs[3*i+1].set_xlabel(jR_text, fontsize=16)
#         axs[3*i+1].set_ylabel(jz_text, fontsize=16)
#         axs[3*i+2].set_xlabel(Lz_text, fontsize=16)
#         axs[3*i+2].set_ylabel(jz_text, fontsize=16)

#         axs[3*i].set_xlim(jRlim[0],jRlim[1])
#         axs[3*i].set_ylim(Lzlim[0],Lzlim[1])
#         axs[3*i+1].set_xlim(jRlim[0],jRlim[1])
#         axs[3*i+1].set_ylim(jzlim[0],jzlim[1])
#         axs[3*i+2].set_xlim(Lzlim[0],Lzlim[1])
#         axs[3*i+2].set_ylim(jzlim[0],jzlim[1])
        
#         h1_out.append([h1_1,h2_1,h3_1])
#         h2_out.append([h1_2,h2_2,h3_2])
#     ###i
        
#     fig.subplots_adjust(wspace=0.3)
#     return fig, axs, h1_out, h2_out
# #def

# def plot_action_map_purity(accs1, accs2, mixture_arr1, mixture_arr2):
#     '''plot_action_map_purity:
    
#     Plot the wierd action map people so love
    
#     Args:
#         accs (array) - Array of different action arrays (array of 3-element lists of actions [jR,jphi,jz])
#         mixture_arr (list) - list of arrays, each of size len(orbs), which denote the fraction
#             of samples from each orbit.Orbit object in orbs that will be plotted
#         mixture_text (array) - Array of text for each mixture instance
#     '''
#     assert len(mixture_arr1) == len(mixture_arr2), 'Mixture array lengths must be the same'
#     n_mix = len(mixture_arr1)
    
#     n_rows = int(np.ceil(n_mix/3))
#     remove_axes_beyond = n_rows*3 - n_mix
    
#     fig = plt.figure( figsize=(18,5*n_rows) )
#     axs = fig.subplots(nrows=n_rows,ncols=3).flatten()
#     for i in range(remove_axes_beyond):
#         axs[-i-1].set_axis_off()
#     ###i
    
#     h1_out = []
#     h2_out = []
    
#     for i in range(n_mix):
        
#         accs_mixture1 = parse_mixture(accs1, mixture_arr1[i])
#         accs_mixture2 = parse_mixture(accs2, mixture_arr2[i])
        
#         jR_1 = np.array([])
#         Lz_1 = np.array([])
#         jz_1 = np.array([])
#         for j in range(len(accs_mixture1)):
#             jR_temp,Lz_temp,jz_temp = accs_mixture1[j]
#             jR_1 = np.append(jR_1,jR_temp)
#             Lz_1 = np.append(Lz_1,Lz_temp)
#             jz_1 = np.append(jz_1,jz_temp)
#         ###j
        
#         jR_2 = np.array([])
#         Lz_2 = np.array([])
#         jz_2 = np.array([])
#         for j in range(len(accs_mixture2)):
#             jR_temp,Lz_temp,jz_temp = accs_mixture2[j]
#             jR_2 = np.append(jR_2,jR_temp)
#             Lz_2 = np.append(Lz_2,Lz_temp)
#             jz_2 = np.append(jz_2,jz_temp)
#         ###j
        
#         Jtot_1 = np.abs(jR_1) + np.abs(Lz_1) + np.abs(jz_1)
#         Jz_JR_norm_1 = (jz_1-jR_1) / Jtot_1
#         Jphi_norm_1 = Lz_1 / Jtot_1
#         Jtot_2 = np.abs(jR_2) + np.abs(Lz_2) + np.abs(jz_2)
#         Jz_JR_norm_2 = (jz_2-jR_2) / Jtot_2
#         Jphi_norm_2 = Lz_2 / Jtot_2    
        
#         h1 = np.histogram2d( Jphi_norm_1, Jz_JR_norm_1, bins=50, range=([-1,1],[-1,1]) )[0]
#         h2 = np.histogram2d( Jphi_norm_2, Jz_JR_norm_2, bins=50, range=([-1,1],[-1,1]) )[0]
#         p1 = np.divide(h1,h1+h2)
        
#         i1 = axs[i].imshow( np.rot90(p1), cmap='rainbow', vmin=0, vmax=1, aspect='auto', extent=(-1,1,-1,1))
#         cm1 = fig.colorbar(i1, ax=axs[i])
#         cm1.set_label(r'Purity: $N_{T}/(N_{T}+N_{F})$')
        
#         axs[i].plot([-1,0], [0,-1], linestyle='solid', color='Black', linewidth=2.)
#         axs[i].plot([0,1],  [-1,0], linestyle='solid', color='Black', linewidth=2.)
#         axs[i].plot([1,0],  [0,1],  linestyle='solid', color='Black', linewidth=2.)
#         axs[i].plot([0,-1], [1,0],  linestyle='solid', color='Black', linewidth=2.)
#         axs[i].set_xlabel(r'$J_{\phi}/J_{\rm tot}$', fontsize=16)
#         axs[i].set_ylabel(r'$(J_{\rm z}-J_{\rm R})/J_{\rm tot}$', fontsize=16)
        
#         axs[i].set_xlim(-1,1)
#         axs[i].set_ylim(-1,1)
        
#         h1_out.append([h1,])
#         h2_out.append([h2,])
#     ###i
#     fig.subplots_adjust(wspace=0.3)
#     return fig, axs, h1_out, h2_out
# #def

# def plot_purity_histograms(hist1, hist2):
#     assert len(hist1) == len(hist2)
#     n_mix = len(hist1)
#     assert len(hist1[0]) == len(hist2[0])
#     n_panel = len(hist1[0])
    
#     fig = plt.figure( figsize=(5*n_panel,5*n_mix) )
#     axs = fig.subplots(nrows=n_mix,ncols=n_panel)
#     if n_panel == 1:
#         axs = axs.reshape((n_mix,1))
        
#     ymax = 0
    
#     for i in range(n_mix):
#         for j in range(n_panel):
            
#             h1 = hist1[i][j]
#             h2 = hist2[i][j]
#             pur1 = np.divide(h1,h1+h2)
#             pur2 = np.divide(h2,h1+h2)
            
#             stat1, be1, binn1 = scipy.stats.binned_statistic(pur1.flatten(), h1.flatten(), bins=10, statistic='sum', range=(0,1))
#             bmid1 = be1[:-1]+be1[0]/2
#             stat2, be2, binn2 = scipy.stats.binned_statistic(pur2.flatten(), h2.flatten(), bins=10, statistic='sum', range=(0,1))
#             bmid2 = be2[:-1]+be2[0]/2
            
#             stat1 = 100*stat1/np.sum(stat1)
#             stat2 = 100*stat2/np.sum(stat2)
            
#             axs[i,j].step(be1[:-1], stat1, where='post', color='DodgerBlue')
#             axs[i,j].step(be2[:-1], stat2, where='post', color='DarkOrange')
#             axs[i,j].plot([be1[-2],be1[-1]], [stat1[-1],stat1[-1]])
#             axs[i,j].plot([be2[-2],be2[-1]], [stat2[-1],stat2[-1]])
            
#             axs[i,j].set_xlabel('Purity')
#             axs[i,j].set_ylabel('Fraction [Percent]')
            
#             if np.max(stat1) > ymax or np.max(stat2) > ymax:
#                 ymax = np.max([np.max(stat1),np.max(stat2)])
#             ##fi
#         ###j
#     ###i
#     for i in range(n_mix):
#         for j in range(n_panel):
#             axs[i,j].set_ylim(0,ymax*1.1)
#             axs[i,j].set_xlim(0,1)
#         ###j
#     ###i
#     return fig,axs
# #def 

# ----------------------------------------------------------------------------

def get_plottable_data(orbs,eELz,accs,mixture_arr,plot_type,phi0=0,seed=0,
                        absolute=False):
    '''get_plottable_data:
    
    Take in all the arrays of data for the different beta models, and output 
    a mixture of the correct type of data to be plotted. plot type can be 
    any of ['vRvT','Toomre','ELz','JRLz','eLz','AD']. Can also be 'Rz' to 
    get radius and height above the disk.
    
    Args:
        orbs (list) - list of orbit.Orbit objects for each beta
        eELz (list) - list of np.array([e,E,Lz]) for each beta
        accs (list) - list of np.array([JR,Lz,Jz]) for each beta
        mixture_arr (list) - 
        plot_type (str) - Type of data to extract for plot
        phi0 (float) - Potential at infinity to subtract from energies [0]
        seed (int) - seed to use when randomly extracting data
        absolute (int) - Use absolute fractions of total amount of data in 
            mixture_arr [False]

    Returns:
        plot_x (np.array) - y coordinate
        plot_y (np.array) - x coordinates
    '''
    
    orbs_mix = parse_mixture(orbs, mixture_arr, seed=seed, absolute=absolute)
    eELz_mix = parse_mixture(eELz, mixture_arr, seed=seed, absolute=absolute)
    accs_mix = parse_mixture(accs, mixture_arr, seed=seed, absolute=absolute)
    n_in_mix = len(orbs_mix)
    
    if plot_type == 'vRvT':
        vR = np.array([])
        vT = np.array([])
        for i in range(n_in_mix):
            vR = np.append( vR, orbs_mix[i].vR().value )
            vT = np.append( vT, orbs_mix[i].vT().value )
        return vR,vT
    
    elif plot_type == 'Toomre':
        vT = np.array([])
        vperp = np.array([])
        for i in range(n_in_mix):
            vT = np.append( vT, orbs_mix[i].vT().value )
            this_vperp = np.sqrt( np.square( orbs_mix[i].vR().value ) +\
                                  np.square( orbs_mix[i].vz().value ) )
            vperp = np.append( vperp, this_vperp )
        return vT,vperp
    
    elif plot_type == 'ELz':
        E = np.array([])
        Lz = np.array([])
        for i in range(n_in_mix):
            E = np.append( E, (eELz_mix[i][1]-phi0)/1e5 )
            Lz = np.append( Lz, eELz_mix[i][2] )
        return Lz,E
    
    elif plot_type == 'JRLz':
        JR = np.array([])
        Lz = np.array([])
        for i in range(n_in_mix):
            JR = np.append( JR, np.sqrt(accs_mix[i][0]) )
            Lz = np.append( Lz, eELz_mix[i][2] )
        return Lz,JR
    
    elif plot_type == 'eLz':
        e = np.array([])
        Lz = np.array([])
        for i in range(n_in_mix):
            e = np.append( e, eELz_mix[i][0] )
            Lz = np.append( Lz, eELz_mix[i][2] )
        return Lz,e
    
    elif plot_type == 'AD':
        Jz_JR = np.array([])
        Jphi = np.array([])
        for i in range(n_in_mix):
            JR,Lz,Jz = accs_mix[i]
            Jtot = np.abs(JR) + np.abs(Lz) + np.abs(Jz)
            Jz_JR_norm = (Jz-JR) / Jtot
            Jphi_norm = Lz / Jtot
            Jz_JR = np.append(Jz_JR, Jz_JR_norm)
            Jphi = np.append(Jphi, Jphi_norm)
        return Jphi,Jz_JR
    
    elif plot_type == 'Rz':
        R = np.array([])
        z = np.array([])
        for i in range(n_in_mix):
            R = np.append( R, orbs_mix[i].R().value )
            z = np.append( z, orbs_mix[i].z().value )
        return R,z
    
    else:
        raise ValueError('plot_type not recognized')
    ##ie
#def

# ----------------------------------------------------------------------------

# Plot specific elements

def add_ELz_boundary(ax,pot,rs,Lz_facs,z=False,E_scale=1.,phi0=0,ro=8,
                     plot_kws={}):
    '''add_ELz_boundary:
    
    Add an ELz boundary to a figure at a number of different radii
    
    Args:
        ax (matplotlib axis object) - axis object
        rs (float or astropy object array) - radial positions
        Lz_facs (array) - Factors that set how wide to plot in Lz
    
    Returns:
        None
    '''
    nrs = len(rs)
    assert nrs == len(Lz_facs), 'Lz_facs must be same shape as rs'
    for i in range(nrs):
        Lz_bound, E_bound = get_ELz_boundary(pot,rs[i],z=z,Lz_fac=Lz_facs[i],
            ro=ro)
        ax.plot(Lz_bound, (E_bound-phi0)/E_scale, **plot_kws)
    ###i
#def

def get_ELz_boundary(pot,r,z=False,Lz_fac=1000,ro=8):
    '''get_ELz_boundary:
    
    Calculate the bounding E-Lz curve at a given radius
    
    Args:
        pot (galpy Potential) - potential
        r (float or astropy unit) - 
        z (bool) - r is in the vertical direction instead of radial
        Lz_fac (float) - Arbitrary factor to multiply radii by to give angular
            momentum end-points
    '''
    if isinstance(r,apu.quantity.Quantity):
        r = r.to(apu.kpc).value
    # Lz values that will be considered
    Lz_lim = [-r*Lz_fac,r*Lz_fac]
    Lz_bound = np.linspace(Lz_lim[0], Lz_lim[1], num=100)*apu.kpc*apu.km/apu.s
    
    if z:
        Pot_bound = potential.evaluatePotentials(pot,0,r/ro).value
    else:
        Pot_bound = potential.evaluatePotentials(pot,r/ro,0).value
    E_bound = (np.power(Lz_bound,2)/(2*r**2)).value + Pot_bound
    
    return Lz_bound, E_bound
#def
    

def add_ellipse(ax,cent,semi,plot_kwargs={}):
    '''add_ellipse:
    
    Add a bounding ellipse to an axis
    
    Args:
        ax (matplotlib axis object) - Axis
        cent (list) - ellipse center [x,y]
        semi (list) - ellipse semi-major axis [x,y]
        plot_kws (dict) - Plot keywords
    
    Returns:
        None
    '''
    ellipse = patches.Ellipse(cent, width=semi[0]*2, height=semi[1]*2, 
        **plot_kwargs)
    ax.add_artist(ellipse)
#def

def add_legend(ax, disk_selection_kws=None, halo_selection_kws=None):
    '''add_legend:
    
    Args:
        ax (matplotlib Axis object) - axis
        disk_selection_kws (dict) - dictionary of disk selection keywords
        halo_selection_kws (dict) - dictionary of halo selection keywords
        
    Returns:
        None
    '''
    cp_disk_selection_kws = copy.deepcopy(disk_selection_kws)
    cp_halo_selection_kws = copy.deepcopy(halo_selection_kws)
    cp_disk_selection_kws.pop('facecolor', None)
    cp_disk_selection_kws.pop('edgecolor', None)
    cp_halo_selection_kws.pop('facecolor', None)
    cp_halo_selection_kws.pop('edgecolor', None)
    ax.plot([],[], **cp_disk_selection_kws)
    ax.plot([],[], **cp_halo_selection_kws)
    legend = ax.legend(loc='lower left', fontsize=7, handlelength=2.2, 
                       frameon=False, ncol=1, handletextpad=0.2)
#def

def axis_limits_and_labels(ax,xlim,ylim,xlabel,ylabel,mixture_text,
                            is_left_edge=False,is_top_edge=False,
                            label_fontsize=8):
    '''axis_spacing_and_labels:
    
    Do the labels and limits for the axis
    
    Args:
        ax (matplotlib axis object) - axis
        xlim (list) - x limits
        ylim (list) - y limits
        xlabel (string) - x label
        ylabel (string) - y label
        is_edge (bool) - Is the left-most edge panel?
    
    Returns:
        None
    '''
    ax.set_xlabel(xlabel, fontsize=label_fontsize, labelpad=-0.1)
    ax.tick_params(axis='both',labelsize=label_fontsize)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if is_left_edge:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    else:
        ax.tick_params(labelleft=False)
    ##ie
    if is_top_edge:
        ax.set_title(mixture_text, fontsize=label_fontsize, loc='center')
    ##fi
#def

def add_selection_boundaries(ax,selection,plot_cent=False,plot_kwargs={}):
    '''add_selection_boundaries:
    
    Draw selection boundaries. Argument is a list of selection boundaries. 
    Each selection is a list of 3 elements, a string denoting the type of 
    boundary ('ellipse' or 'line'), then two arrays which each have two 
    elements. The two arrays are either [xcent,ycent] and [xsemi,ysemi]
    (for ellipse) or [x1,x2] and [y1,y2] (for line).
    
    Args:
        ax (matplotlib axis object) - axis
        selection (list) - list of selections
        plot_cent (bool) - Plot the center of the bounding region (ellipse only)
        plot_kwargs (dict) - keyword arguments for the plot
    
    Returns:
        None
    '''
    n_selec = len(selection)
    for i in range(n_selec):
        cp_plot_kwargs = copy.deepcopy(plot_kwargs)
        if selection[i][0] == 'ellipse':
            cp_plot_kwargs.pop('color', None)
            _,cent,semi = selection[i]
            add_ellipse(ax,cent,semi,plot_kwargs=cp_plot_kwargs)
            if plot_cent:
                ax.scatter(cent[0], cent[1], marker='x', s=10, color='Black')
        if selection[i][0] == 'line':
            cp_plot_kwargs = copy.deepcopy(plot_kwargs)
            cp_plot_kwargs.pop('facecolor', None)
            cp_plot_kwargs.pop('edgecolor', None)
            _,xs,ys = selection[i]
            ax.plot(xs,ys,**cp_plot_kwargs)
        ##fi
    ###i
#def

def add_diamond_boundary(ax,dedge=1.2):
    '''add_diamond_boundary:
    
    Add the diamond-shaped boundary and other things to the action diamond
    panels.
    
    Args:
        ax (matplotlib axis object) - Axis
    
    Returns:
        None
    '''
    
    # Keywords
    diamond_zorder = 5
    
    # Plot the white background to clean up the edges of the diamond?
    white_background = True
    if white_background:
        blank_contour_pts_R = [[0.,-1.],[1.,0.],[0,1.],[-0.5,dedge],
                                [dedge,dedge],[dedge,-dedge],[-0.5,-dedge],
                                [0.,-1.]]
        blank_contour_pts_L = [[0.,-1.],[-1.,0.],[0,1.],[0.5,dedge],
                                [-dedge,dedge],[-dedge,-dedge],[0.5,-dedge],
                                [0.,-1.]]
        blank_contour_poly_R = patches.Polygon(blank_contour_pts_R, 
            closed=True, fill=True, zorder=2, edgecolor='None', 
            facecolor='White')
        blank_contour_poly_L = patches.Polygon(blank_contour_pts_L, 
            closed=True, fill=True, zorder=2, edgecolor='None', 
            facecolor='White')
        ax.add_artist(blank_contour_poly_R)
        ax.add_artist(blank_contour_poly_L)
        
        # Do artificial ticks because the ticks get overwritten by the 
        # drawing of the white patch
        ticklen = 2*dedge/20
        xticks = [-1,0,1]
        yticks = [-1,0,1]
        tickwidth = 1.
        for i in range(len(xticks)):
            ax.plot( [xticks[i],xticks[i]], [-dedge,-dedge+ticklen], 
                linewidth=tickwidth, color='Black', zorder=diamond_zorder )
            ax.plot( [xticks[i],xticks[i]], [dedge,dedge-ticklen], 
                linewidth=tickwidth, color='Black', zorder=diamond_zorder )
            ax.plot( [-dedge,-dedge+ticklen], [yticks[i],yticks[i]], 
                linewidth=tickwidth, color='Black', zorder=diamond_zorder )
            ax.plot( [dedge,dedge-ticklen], [yticks[i],yticks[i]], 
                linewidth=tickwidth, color='Black', zorder=diamond_zorder )
        ###i
    #def
    
    # Plot the diamond itself
    ax.plot([-1,0], [0,-1], linestyle='solid', color='Black', linewidth=0.5, 
            zorder=diamond_zorder)
    ax.plot([0,1],  [-1,0], linestyle='solid', color='Black', linewidth=0.5, 
            zorder=diamond_zorder)
    ax.plot([1,0],  [0,1],  linestyle='solid', color='Black', linewidth=0.5, 
            zorder=diamond_zorder)
    ax.plot([0,-1], [1,0],  linestyle='solid', color='Black', linewidth=0.5, 
            zorder=diamond_zorder)
#def

# ----------------------------------------------------------------------------

# Check for points inside regions in kinematic spaces

def is_in_scaled_selection(x,y,selection,factor=1):
    '''is_in_scaled_selection:
    
    Lightweight wrapper of is_in_ellipse to allow scaling the ellipse by a 
    constant factor in either dimension or both 
    
    Args:
        x (array) - X coordinates
        y (array) - Y coordinates
        selection (list) - List of selections. Must be ellipses
        factor (float or array) - 2-element array of scaling factors for x 
            direction and y direction. If is float then will be cast as 
            np.array([factor,factor])
    
    Returns:
        is_inside (array) - Boolean array same size as x and y
    '''
    if isinstance(factor,numbers.Number):
        factor = np.array([factor,factor],dtype='float')
    assert len(x) == len(y), 'x and y must be same shape'
    is_in_ellipse_bools = np.zeros_like(x,dtype='bool')
    n_selec = len(selection)
    for i in range(n_selec):
        assert selection[i][0] == 'ellipse', 'selections must be ellipses'
        _,cent,semi = selection[i]
        semi_in = [semi[0]*factor[0],semi[1]*factor[1]]
        is_in_ellipse_bools = is_in_ellipse_bools |\
                              is_in_ellipse(x,y,cent,semi_in)
    ###i
    return is_in_ellipse_bools
#def

def is_in_ellipse(x,y,ellipse_cent,ellipse_ab):
    '''is_in_ellipse:
    
    Determine if points x,y are in an ellipse
    
    Args:
        x (array) - X coordinates
        y (array) - Y coordinates
        ellipse_cent (list) - 2-element list of ellipse center
        ellipse_ab (list) - 2-element list of semi-major axes
    
    Returns:
        is_inside (array) - Boolean array same size as x and y
    '''
    a,b = ellipse_ab
    xo,yo = ellipse_cent
    if isinstance(x,apu.quantity.Quantity):
        x = x.value
    if isinstance(y,apu.quantity.Quantity):
        y = y.value
    ##fi
    elliptical_dist = np.power(x-xo,2)/a**2 + np.power(y-yo,2)/b**2
    return elliptical_dist < 1
#def

# ----------------------------------------------------------------------------

# Lines equations

def line_equation(x,m,b):
    '''line_equation:
    
    Args:
        x (array) - x values
        m (float) - slope
        b (float) - y-intercept
    
    Returns:
        y (array) - y values
    '''
    return m*x + b
#def


def get_params_from_line(xs,ys):
    '''get_params_from_line:
    
    Get line parameters from 2 points
    
    Args:
        xs (list) - 2-element list of x points
        ys (list) - 2-element list of y points
    
    Returns:
        m (float) - slope
        b (float) - y-intercept
    '''
    m = (ys[1]-ys[0])/(xs[1]-xs[0])
    b = ys[0] - m*xs[0]
    return m,b
#def

# ----------------------------------------------------------------------------

# Colormaps

class colors(object):
    '''Stolen from Ted Mackereth halo-mass github'''
    def __init__(self):

        # colour table in HTML hex format
        self.hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', 
                   '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466',
                   '#4477AA']

        self.greysafecols = ['#809BC8', '#FF6666', '#FFCC66', '#64C204']

        self.xarr = [[12], 
                [12, 6], 
                [12, 6, 5], 
                [12, 6, 5, 3], 
                [0, 1, 3, 5, 6], 
                [0, 1, 3, 5, 6, 8], 
                [0, 1, 2, 3, 5, 6, 8], 
                [0, 1, 2, 3, 4, 5, 6, 8], 
                [0, 1, 2, 3, 4, 5, 6, 7, 8], 
                [0, 1, 2, 3, 4, 5, 9, 6, 7, 8], 
                [0, 10, 1, 2, 3, 4, 5, 9, 6, 7, 8], 
                [0, 10, 1, 2, 3, 4, 5, 9, 6, 11, 7, 8]]

    # get specified nr of distinct colours in HTML hex format.
    # usage: colour_list = safe_colours.distinct_list(num_colours_required)
    # returns: list of distinct colours in HTML hex
    def distinct_list(self,nr):

        # check if nr is in correct range
        if nr < 1 or nr > 12:
            print("wrong nr of distinct colours!")
            return

        # get list of indices
        lst = self.xarr[nr-1]
        
        # generate colour list by stepping through indices and looking them up
        # in the colour table
        i_col = 0
        col = [0] * nr
        for idx in lst:
            col[i_col] = self.hexcols[idx]
            i_col+=1
        return col

    # Generate a dictionary of all the safe colours which can be addressed by colour name.
    def distinct_named(self):
        cl = self.hexcols

        outdict = {'navy':cl[0],\
                   'cyan':cl[1],\
                   'turquoise':cl[2],\
                   'green':cl[3],\
                   'olive':cl[4],\
                   'sandstone':cl[5],\
                   'coral':cl[6],\
                   'maroon':cl[7],\
                   'magenta':cl[8],\
                   'brown':cl[9],\
                   'skyblue':cl[10],\
                   'pink':cl[11],\
                   'blue':cl[12]}

        return outdict


    # For making colourmaps.
    # Usage: cmap = safe_colours.colourmap('rainbow')
    def colourmap(self,maptype,invert=False):

        if maptype == 'diverging':
            # Deviation around zero colormap (blue--red)
            cols = []
            for x in np.linspace(0,1, 256):
                rcol = 0.237 - 2.13*x + 26.92*x**2 - 65.5*x**3 + 63.5*x**4 - 22.36*x**5
                gcol = ((0.572 + 1.524*x - 1.811*x**2)/(1 - 0.291*x + 0.1574*x**2))**2
                bcol = 1/(1.579 - 4.03*x + 12.92*x**2 - 31.4*x**3 + 48.6*x**4 - 23.36*x**5)
                cols.append((rcol, gcol, bcol))

            if invert==True:
                cols = cols[::-1]

            return plt.get_cmap(matplotlib.colors.LinearSegmentedColormap.from_list("PaulT_plusmin", cols))

        elif maptype == 'heat':
            # Linear colormap (white--red)
            cols = []
            for x in np.linspace(0,1, 256):
                rcol = (1 - 0.392*(1 + erf((x - 0.869)/ 0.255)))
                gcol = (1.021 - 0.456*(1 + erf((x - 0.527)/ 0.376)))
                bcol = (1 - 0.493*(1 + erf((x - 0.272)/ 0.309)))
                cols.append((rcol, gcol, bcol))

            if invert==True:
                cols = cols[::-1]

            return plt.get_cmap(matplotlib.colors.LinearSegmentedColormap.from_list("PaulT_linear", cols))

        elif maptype == 'rainbow':
            # Linear colormap (rainbow)
            cols = []
            for x in np.linspace(0,1, 254):
                rcol = (0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3)
                gcol = 0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6
                bcol = 1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5)
                cols.append((rcol, gcol, bcol))

            if invert==True:
                cols = cols[::-1]

            return plt.get_cmap(matplotlib.colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols))

        else:
            raise KeyError('Please pick a valid colourmap, options are "diverging", "heat" or "rainbow"')
#cls