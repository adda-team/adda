'''
# This code represents the comparison of scattering intensities calculated 
# with GLMT (reference data) and DDA (with ADDA) for the scattering of 
# Bessel beam (CS type) by a coated sphere.
# Reference data were kindly provided by Zhuyang Chen 
# (DOI:10.1088/2040-8978/16/5/055701).
'''

import os
import matplotlib.pyplot as plt
import bb_module as bb

option = 2
bb.adda_run(option)
theta1,iper1,ipar1 = bb.extractData('dda',option) 
option = 3
bb.adda_run(option)
theta2,iper2,ipar2 = bb.extractData('dda',option) 
option = 4
bb.adda_run(option)
theta3,iper3,ipar3 = bb.extractData('dda',option)   


# data visualisation
fig = plt.figure()

# Perpendicular scattering intensity (1)
ax = fig.add_subplot(121)
bb.plotData3(theta1,iper1,theta2,iper2,theta3,iper3,1,ax)

#image of a scattering particle
img = plt.imread('particles/cube.png')
newax = fig.add_axes([0.22,0.50,0.2,0.2], anchor='NE', zorder=1)
newax.imshow(img)
newax.axis('off')

# Parallel scattering intensity (2)
ax = fig.add_subplot(122)
bb.plotData3(theta1,ipar1,theta2,ipar2,theta3,ipar3,2,ax)

plt.tight_layout()

# data save
os.makedirs('saved', exist_ok=True)
plt.savefig('saved/fig15_cube.pdf', bbox_inches='tight')

plt.show()