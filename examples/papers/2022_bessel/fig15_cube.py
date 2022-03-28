'''
# The code represents the figure 15 with demonstation of incident field 
# intensity profiles for 4 types of Bessel beams (CS, CS', TEL, and TML) and 
# their scattering intensities on a cube.
'''

import os
import matplotlib.pyplot as plt
import bb_module as bb

if not os.path.exists('dda'):
    os.makedirs('dda')

option = 2
bb.adda_run(option)
theta1,iper1,ipar1 = bb.extractData('dda',option) 
xd1,yd1,ed1,minz1 = bb.extractField(option)
option = 3
bb.adda_run(option)
theta2,iper2,ipar2 = bb.extractData('dda',option) 
xd2,yd2,ed2,minz2 = bb.extractField(option)
option = 4
bb.adda_run(option)
theta3,iper3,ipar3 = bb.extractData('dda',option)   
xd3,yd3,ed3,minz3 = bb.extractField(option)
option = 5
bb.adda_run(option)
theta4,iper4,ipar4 = bb.extractData('dda',option) 
xd4,yd4,ed4,minz4 = bb.extractField(option)  

# data visualisation
fig = plt.figure(figsize=(9,5))
spec = fig.add_gridspec(ncols=4, nrows=2,height_ratios=[1, 1.4],left=0.04, right=0.99,bottom=0.1,top=0.99, wspace=0.09,hspace=0.1)

# Visualisation of the intensity of the incident electric field
ax = fig.add_subplot(spec[0],projection='3d')
bb.plotField(xd1,yd1,ed1,minz1,2,ax)
ax = fig.add_subplot(spec[1],projection='3d')
bb.plotField(xd2,yd2,ed2,minz2,3,ax)
ax = fig.add_subplot(spec[2],projection='3d')
bb.plotField(xd3,yd3,ed3,minz3,4,ax)
ax = fig.add_subplot(spec[3],projection='3d')
bb.plotField(xd4,yd4,ed4,minz4,5,ax)
# Parallel scattering intensity (1)
ax = fig.add_subplot(spec[1:,:-2])
bb.plotData4(theta1,ipar1,theta2,ipar2,theta3,ipar3,theta4,ipar4,1,ax)
#image of a scattering particle
img = plt.imread('particles/cube.png')
newax = fig.add_axes([0.08,0.13,0.16,0.16], anchor='NE', zorder=1)
newax.imshow(img)
newax.axis('off')
# Perpendicular scattering intensity (1)
ax = fig.add_subplot(spec[1:,2:])
bb.plotData4(theta1,iper1,theta2,iper2,theta3,iper3,theta4,iper4,2,ax)

# data save
os.makedirs('saved', exist_ok=True)
plt.savefig('saved/fig15_cube.pdf', bbox_inches='tight')

plt.show()