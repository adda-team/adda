'''
# This code represents the extrapolation for the comparison of backscattering
# intensities calculated with GLMT (reference data) and DDA (with ADDA) for
# the scattering of Bessel beam (CS type) by a sphere.
# Reference data were kindly provided by Zhuyang Chen
# (DOI:10.1088/2040-8978/16/5/055701).
'''

import os
import matplotlib.pyplot as plt
import bb_module as bb

if not os.path.exists('dda'):
    os.makedirs('dda')
if not os.path.exists('dda/extrapolation'):
    os.makedirs('dda/extrapolation')

run = 1 # 0- do not run adda; 1- run adda
theta = 180
# Attention! Large grids require high computational power!!!
grids = [64,80,96,112,128,160,192,224,256]
# less accurate extrapolation
#grids = [32,40,48,56,64,80,96,112,128]

if run == 1:
    bb.adda_run_grids(grids)

# extraction of ADDA results
ig = bb.extractData_grids(grids,theta)
# extraction of reference data (GLMT)
ref_thetaper,ref_iper,ref_thetapar,ref_ipar = bb.extractData('glmt',0)

# extrapolation
nextr = 9 #number of points for extrapolation
Nf = len(grids) # max grid
N1 = Nf-nextr # from what point we start extrapolation
N0 = Nf-N1-nextr # to show some extra points (small grids)

fig = plt.figure()
ys = bb.axesGen(grids[N0:N1])
plt.plot(ys,ig[N0:N1], marker="o", linestyle="none", color = (121/256,251/256,186/256))
ind = ref_thetaper.index(theta)
bb.fitDotPlot(grids[N1:Nf:],ig[N1:Nf:],ref_iper[ind])

# data save
os.makedirs('saved', exist_ok=True)
plt.savefig('saved/fig14_extrapolation.pdf', bbox_inches='tight')

plt.show()
