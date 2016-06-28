## ====================================================
## Purpose:  Animate DIMES floats in 3D using Mayavi
## Author :  Dhruv Balwada
## ====================================================

import os
# Go to appropriate data directory
os.chdir('/Users/dhruvbalwada/Work/Scripts_mat/3Dvis')
import numpy
import scipy
import mayavi
from mayavi import mlab
from scipy import io
import math

# Load files
bath = scipy.io.loadmat('bathylatlon.mat')
f    = scipy.io.loadmat('floatslatlon.mat')


# Scaling factor for bathymetry
scl=0.0015

# Setup initial scene
mayavi.mlab.figure(1, bgcolor=(0, 0, 0))
mayavi.mlab.surf((bath['bx']),bath['by'],bath['bz'],warp_scale=scl,opacity=1,colormap='gist_earth') 

fx=f['fx']
fy=f['fy']
fz=f['fz']

s=numpy.array([[row for col in range(114)] for row in range(31)])

@mlab.animate(delay=10)
def anim():
	
	# First set of floats in to scene
	x=fx[0:31,0:114]
	y=fy[0:31,0:114]
	z=fz[0:31,0:114]*scl

	pts=mayavi.mlab.points3d(x,y,z,s,colormap='Reds', resolution=16)
	
	srcpt=pts.mlab_source
	v1=90 # azimuthal angle
	v2=50  # elevation
	v3=40 # distance
	cx=-90
	cy=-55
	cz=-2.68
		
	for i in range(100,980):
			
			mayavi.mlab.view(v1, v2, v3, focalpoint=(cx,cy,cz))
			
			if i<=100 :
				cx=cx+0.1
			elif (i>100 and i<=300):
				v2=v2+42./200
				v1=v1-3./20
				v3=v3-20./200
				cx=cx-5./200
			elif (i>300 and i<=500): 
				v1=v1+35./200
				v3=v3+10./200
				cx=cx+10./200	
			elif (i>500 and i<=700): 
				cx=cx+10./200
				v3=v3+15./200
				v2=v2-47./200
			elif (i>700 and i<=900):
				v1=v1-5./200
				v2=v2-35./200
				v3=v3+20./200
			else:
				print(i)
				
			# Progressively move floats 
			x=fx[i:i+31,0:114]
			y=fy[i:i+31,0:114]
			z=fz[i:i+31,0:114]*scl
					
			srcpt.reset(x=x, y=y, z=z)
			
			# Save file to be later combined in mm_frame
			# Changed later to combine in imovie
			filename=("./mm_frames/mm_frame_%03i.png" % i)
			mayavi.mlab.savefig(filename)
	yield

# Add a slightly differnt view at the end.
@mlab.animate(delay=10)
def anim2():
	for i in range(980,994):
			for l in range(114):
				a=fx[:,l]
				b=fy[:,l]
				c=fz[:,l]*scl
				o=(i-978)/(994-978)
				mlab.plot3d(a,b,c,opacity=o)
				filename=("./mm_frames/mm_frame_%03i.png" % i)
				mayavi.mlab.savefig(filename)
		
	yield

a=anim()

