#!/usr/bin/python
# -*- coding: utf-8 -*-

# The code sets up a file structure and creates a file with the run commands for ADDA. 
# The file can then be run in the same location or moved to another location for 
# running it. The current setup generates 8 commands for the surface case or 2 for
# the free-space case. The ADDA commands include (in given order):
# shape (sphere), dipoles per wavelength (dpl), radius (eq_rad), directory for saving
# the output (dir), incident propagation direction (prop), scattering grid input file
# (scat_grid_inp and store_scat_grid; make sure it's in the same directory), particle
# refractive index (m, real and imaginary parts), iteration precision limit (eps) and
# other computational definitions needed for a high refractive index (int, pol, iter),
# and in the surface case, particle center's height from the surface and the substrate's 
# refractive index (surf).   

from numpy import pi,cos, sin, linspace
import os

# size parameters
xs = [1.0,6.0]
# incident propagation angles
incidenceprop = linspace(0.0, 60.0, 4)   # (first, last, number of angles)
# include a substrate (0) or free-space case (1), or both (2)
mode = 2
# run file name
run_file_name = 'run_sphere'

def frontmatter(f):
    f.write('#!/bin/bash\n# Uses common script to find ADDA executable, look inside it for details\n')
    f.write('ADDA=$(../../find_adda)\nif [ $? -ne 0 ]; then\n  exit 1\nfi\n\n')

def main_surface():
    # write run commands into a file
    f = open(run_file_name, "w")
    frontmatter(f)
    for sp in xs:
        if not os.path.exists('SurfSphere/Sphere_m217i0004_x%d' % (int(sp*10))):
            os.makedirs('SurfSphere/Sphere_m217i0004_x%d' % (int(sp*10)))    
        for orient in incidenceprop: 
            inc = orient * pi / 180.0    # incidence angle in radians 
            f.write('"$ADDA" -shape sphere -dpl 32.0 -eq_rad %.1f ' % (sp))
            f.write('-dir SurfSphere/Sphere_m217i0004_x%d/Prop-or%d ' % (int(sp*10),orient))
            f.write('-prop 0.0 %.4f %.4f ' % (sin(inc), -cos(inc)))   # propagation
            f.write('-scat_grid_inp scat_grid_directions.dat -store_scat_grid ')
            f.write('-m 2.17 0.004 -eps 4 -int fcd -iter qmr2 -pol fcd ') # refr. index
            f.write('-surf %.1f 1.55 0.004 \n' % (sp)) # height, substrate refr. index
    f.close()

def main_freespace():
    if not os.path.exists('NoSurfSphere'):
        os.makedirs('NoSurfSphere')
        
    if mode == 1:
        f = open(run_file_name, "w")
        frontmatter(f)
    elif mode == 2:
        f = open(run_file_name, "a")
    
    for sp in xs:
        f.write('"$ADDA" -shape sphere -dpl 32.0 -eq_rad %.1f ' % (sp))
        f.write('-dir NoSurfSphere/Sphere_m217i0004_x%d ' % (int(sp*10)))
        f.write('-prop 0.0 0.0 -1.0 ')  # propagation
        f.write('-scat_grid_inp scat_grid_directions.dat -store_scat_grid ')
        f.write('-eps 4 -m 2.17 0.004 -int fcd -iter qmr2 -pol fcd \n') # refr. index 
    f.close()
            
if __name__ == '__main__':
    if mode == 0 or mode == 2:
        main_surface()
    if mode == 1 or mode == 2:
        main_freespace()
