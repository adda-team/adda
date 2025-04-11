#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Generates Figure 8 for the publication
# Virkki, A. and Yurkin, M. (2025), Microwave scattering by rough polyhedral particles on 
# a surface. Reads scattering matrix grid arrays from ADDA and plots a mosaic of the
# scattering-matrix elements (1,1), (1,2), (2,2), (3,3), and (4,4) in the scattering plane  
# for incidence angles 0°, 20°, 40°, and 60°.
# The code reads by default zenith angles up to 88, located at  
# 'SurfSphere/Sphere_m217i0004_x{size parameter}/Prop-or{incidence index}/mueller_scatgrid'.
# Free-space case can be excluded by setting freespace = False.

import numpy as np
from matplotlib import pyplot as plt

sp = 6.0           # size parameter
freespace = True   # include the free space case? Use True/False for yes/no
figname = 'fig8.pdf'
            
def draw(inc0, ax, mn, mx):
    # Draw reference lines 
    ax.plot([inc0,inc0], [mn,mx], '--', color=[0.4,0.4,0.4], lw=2)  
    ax.plot([-inc0,-inc0], [mn,mx], '--', color=[0.7,0.7,0.7], lw=2)

def main():

    Na = 89
    Nf = Na * 73

    angles = np.linspace(-Na+1,Na-1,177)  # From -88 to 88
    
    fig2, ax2 = plt.subplots(nrows=5,ncols=4,sharex=True,figsize=(16,9))
    fig2.subplots_adjust(top=0.95, right=0.96, left=0.08, bottom=0.12, hspace=0.2)
    ax2[0,0].set_title(r'$\theta_i$: 0°', fontsize=16)
    ax2[0,1].set_title(r'$\theta_i$: 20°', fontsize=16) 
    ax2[0,2].set_title(r'$\theta_i$: 40°', fontsize=16)
    ax2[0,3].set_title(r'$\theta_i$: 60°', fontsize=16)

    ax2[4,0].set_xlabel(r'$\theta_z$ [°]', fontsize=16)
    ax2[4,1].set_xlabel(r'$\theta_z$ [°]', fontsize=16) 
    ax2[4,2].set_xlabel(r'$\theta_z$ [°]', fontsize=16)
    ax2[4,3].set_xlabel(r'$\theta_z$ [°]', fontsize=16)

    ax2[0,0].set_ylabel(r'$F_{11}$', fontsize=16)
    ax2[1,0].set_ylabel(r'$-F_{12}/F_{11}$', fontsize=16) 
    ax2[2,0].set_ylabel(r'$F_{22}/F_{11}$', fontsize=16)
    ax2[3,0].set_ylabel(r'$F_{33}/F_{11}$', fontsize=16)
    ax2[4,0].set_ylabel(r'$F_{44}/F_{11}$', fontsize=16)
    
    for a in ax2.flat:
        a.tick_params(axis='both',labelsize=16)
        a.set_xlim([-89,89])
        
    zext = []
    
    # Free-space sphere case for reference
    if freespace == True:
        ScMb = np.zeros((Nf,5))
        ScMf = np.zeros((Nf,5))
        Nb = int(181*73 - Nf)
        data = np.loadtxt('NoSurfSphere/Sphere_m217i0004_x%d/mueller_scatgrid' % (10*sp), skiprows=1)
            
        ScMf[:,0] = data[:Nf,2]
        ScMf[:,1] = data[:Nf,3]    # F12
        ScMf[:,2] = data[:Nf,7]    # F22
        ScMf[:,3] = data[:Nf,12]    # F33
        ScMf[:,4] = data[:Nf,17]    # F44
    
        ScMb[:,0] = data[Nb:,2]
        ScMb[:,1] = data[Nb:,3]    # F12
        ScMb[:,2] = data[Nb:,7]    # F22
        ScMb[:,3] = data[Nb:,12]    # F33
        ScMb[:,4] = data[Nb:,17]    # F44
        
        valuesb = np.reshape(ScMb,(Na,73,5))
        valuesf = np.reshape(ScMf,(Na,73,5))
        valuesb = np.flip(valuesb,axis=0)
    
        ##-- Plot ------------------------------------------------
        # Polarizations
        for i in range(1,5):
            polb = valuesb[:,:,i] / valuesb[:,:,0]            
            back = np.hstack([np.flip(polb[:,54]),polb[1:,18]])
            polf = valuesf[:,:,i] / valuesf[:,:,0]            
            forw = np.hstack([np.flip(polf[:,54]),polf[1:,18]])
    
            if i == 1:                
                ax2[i,0].plot(angles, -forw, '--', color='r', lw=2)
                ax2[i,0].plot(angles, -back, '-.', color='b', lw=2)
            else:
                ax2[i,0].plot(angles, forw, '--', color='r', lw=2, label='free space (f)')
                ax2[i,0].plot(angles, back, '-.', color='b', lw=2, label='free space (b)')     
                       
        f11 = np.hstack([np.flip(valuesf[:,54,0]),valuesf[1:,18,0]])
        b11 = np.hstack([np.flip(valuesb[:,54,0]),valuesb[1:,18,0]])
        ax2[0,0].semilogy(angles, f11, 'r--', lw=2)
        ax2[0,0].semilogy(angles, b11, 'b-.', lw=2)
        
        mn, mx = min(b11), max(f11)    
        zext.append(0.8*mn)
        zext.append(1.3*mx)  

    # orientation indices and the orientations in radians
    iang = [0, 20, 40, 60]
    
    for ori in range(4):
        ScM0 = np.zeros((Nf,5)) #6497

        data = np.loadtxt('SurfSphere/Sphere_m217i0004_x%d/Prop-or%d/mueller_scatgrid' % (10*sp,20*ori), skiprows=1)
        
        ScM0[:,0] = data[:Nf,2]
        ScM0[:,1] = data[:Nf,3]    # F12
        ScM0[:,2] = data[:Nf,7]    # F22
        ScM0[:,3] = data[:Nf,12]    # F33
        ScM0[:,4] = data[:Nf,17]    # F44
    
        valuesz = np.reshape(ScM0,(Na,73,5))

#     #-- Plot ------------------------------------------------
        # Polarizations
        for i in range(1,5):
            polz = valuesz[:,:,i] / valuesz[:,:,0]            
            zen = np.hstack([np.flip(polz[:,54]),polz[1:,18]])
                        
            draw(iang[ori],ax2[i,ori], -1.05, 1.05)
            
            if i == 1:                
                ax2[i,ori].plot(angles, -zen, '-', color='k', lw=2)
            else:
                ax2[i,ori].plot(angles, zen, '-', color='k', lw=2, label='surface')
            ax2[i,ori].set_ylim([-1.01,1.02])
        
        # Intensity
        zen = np.hstack([np.flip(valuesz[:,54,0]),valuesz[1:,18,0]])
        
        mn, mx = min(zen), max(zen)    
        zext.append(0.9*mn)
        zext.append(1.1*mx)        
#         draw(iang[ori],ax2[0,ori], mn, mx)
        ax2[0,ori].semilogy(angles, zen, '-', color='k', lw=2)
    
    for k in range(3):
        ax2[0,k].set_ylim([min(zext), max(zext)])

    mn, mx = min(zext), max(zext)
    for k in range(4):           
        draw(iang[k],ax2[0,k], 1.1*mn, 0.9*mx)
        for i in range(1,5):
            draw(iang[k],ax2[i,k], -1.1, 1.1)
        ax2[0,k].set_ylim([mn, mx])
    
    if freespace == True:
        ax2[2,0].legend(loc=0, prop={'size':14})
        
    fig2.savefig(fname=figname,format='pdf')
        
    plt.show()    
	
if __name__ == '__main__':
    main()