#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Generates Figure 11A (sp = 1.0) or Figure 11C (sp = 6.0) for the publication
# Virkki, A. and Yurkin, M. A. (2025), Microwave scattering by rough polyhedral particles on a surface.
# Reads scattering matrix grid arrays computed using ADDA and plots the log10-scale
# intensity (1,1-element) in polar plots for the zenith and nadir views and one for both hemispheres
# in the scattering plane (along the horizontal axis of the polar plots) for incidence angles 
# 0°, 20°, 40°, and 60°. Other elements can be plotted but the interpretation
# of some elements could be ambiguous out of the scattering plane.
# The code assumes that the data with 181 zenith angles and 73 azimuth angles (Naz) are  
# located at 'SurfSphere/Sphere_m217i0004_x{size parameter}/Prop-or{incidence index}/mueller_scatgrid'.
# The code also works with Naz = 25, if the data are so computed.

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

sp = 1.0   # size parameter
Naz = 73   # number of azimuths 
figname = 'fig11A.pdf'
            
def draw(inc0, ax, mn, mx):
    # Draw reference lines 
    
    inc = np.pi*0.5 + inc0
    rfl = np.pi*0.5 - inc0
    rfr = 1.5*np.pi + np.arcsin(np.sin(inc0)/1.55)
        
    ax.plot([inc,inc], [mn,mx], 'k--', lw=1)  
    ax.plot([rfr,rfr], [mn,mx], 'k--', lw=1)
    ax.plot([rfl,rfl], [mn,mx], 'k--', lw=1)

def main():
    
    Na = 81
    Nf = int(Na * Naz)
    Nb = int(181*Naz - Nf)
    angles = np.linspace(0,6.283,360)
    pp2 = np.pi/2
    
#     #-- Contour plot ------------------------------------------------
    # zenith angles up to 80 degrees, nadir angles up to 40 degrees
    zeniths = np.arange(Na)
    azimuths = np.radians(np.arange(0, 361, 360/(Naz-1)))
    rz, thetaz = np.meshgrid(zeniths, azimuths)
    rn, thetan = np.meshgrid(zeniths[:41], azimuths)
    
    fig, ax = plt.subplots(nrows=3,ncols=4, subplot_kw=dict(projection='polar'), figsize=(7,6))
    fig.subplots_adjust(top=0.92, right=0.96, left=0.06, bottom=0.06, hspace=0.36)
    ax[0,0].set_title(r'$\theta_i$: 0°', fontsize=14)
    ax[0,1].set_title(r'$\theta_i$: 20°', fontsize=14)
    ax[0,2].set_title(r'$\theta_i$: 40°', fontsize=14)
    ax[0,3].set_title(r'$\theta_i$: 60°', fontsize=14)
    
    for k in range(4):
        for n in range(3):
            ax[n,k].set_xticklabels([])
            ax[n,k].set_rlabel_position(-30)
    
    # orientation indices and the orientations in radians
    iang = [0, np.pi/9.0, np.pi/4.5, np.pi/3.0]
    
    for ori in range(4):
        ScM0 = np.zeros((Nf,2))
        ScM180 = np.zeros((Nf,2))

        data = np.loadtxt('SurfSphere/Sphere_m217i0004_x%d/Prop-or%d/mueller_scatgrid' % (10*sp,20*ori), skiprows=1)
        ScM0[:,0] = data[:Nf,2]
        ScM180[:,0] = data[Nb:,2]
#         ScM0[:,1] = data[:Nf,17]    # choose a polarization element here
#         ScM180[:,1] = data[Nb:,17]  # from columns 3-17
    
        valuesz = np.reshape(ScM0,(Na,Naz,2))
        valuesn = np.reshape(ScM180,(Na,Naz,2))
        valuesn = np.flip(valuesn,axis=0)

#     #-- Plot ------------------------------------------------
        if sum(ScM0[:,1]) == 0:
            ax[0,ori].contourf(thetaz-pp2, rz, np.log10(valuesz[:,:,0].T), cmap=plt.cm.Greens_r)
            ax[1,ori].contourf(thetan-pp2, rn, np.log10(valuesn[:41,:,0].T), cmap=plt.cm.Blues_r)
            
            if Naz == 73:
                zen = np.hstack([np.flip(valuesz[:,54,0]),valuesz[1:,18,0]])
                nad = np.hstack([np.flip(valuesn[:41,54,0]),valuesn[1:41,18,0]]) 
            elif Naz == 25:
                zen = np.hstack([np.flip(valuesz[:,18,0]),valuesz[1:,6,0]])
                nad = np.hstack([np.flip(valuesn[:41,18,0]),valuesn[1:41,6,0]])             
            zen = np.log10(zen)
            nad = np.log10(nad)           

        else:
            polz = valuesz[:,:,1] / valuesz[:,:,0]
            poln = valuesn[:,:,1] / valuesn[:,:,0]
            ax[0,ori].contourf(thetaz-pp2, rz, polz.T, cmap=plt.cm.Greens_r)
            ax[1,ori].contourf(thetan-pp2, rn, poln[:41,:].T, cmap=plt.cm.Blues_r)
                        
            if Naz == 73:
                zen = np.hstack([np.flip(polz[:,54]),polz[1:,18]])
                nad = np.hstack([np.flip(poln[:41,54]),poln[1:41,18]])
            elif Naz == 25:
                zen = np.hstack([np.flip(polz[:,18]),polz[1:,6]])
                nad = np.hstack([np.flip(poln[:41,18]),poln[1:41,6]])
            
        mn, mx = min([min(zen),min(nad)]), max([max(zen),max(nad)])
            
        draw(iang[ori],ax[2,ori], mn, mx)
                        
        ax[2,ori].plot(np.flip(angles[10:171]), zen, 'g-', lw=2)
        ax[2,ori].plot(angles[230:311], nad, 'b-', lw=2)  
        
    fig.savefig(fname=figname,format='pdf')
    plt.show()    
	
if __name__ == '__main__':
    main()