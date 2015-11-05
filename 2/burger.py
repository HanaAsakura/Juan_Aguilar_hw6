from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np
import PIL
from PIL import Image
import os

u = (np.loadtxt("burger_data.dat"))[:,0]
v = (np.loadtxt("burger_data.dat"))[:,1]

length = 41**2
nx = 41
ny = 41
x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)

for i in range (0,501):
    
    if (i % 10 == 0):   
 
    	uu = u[length*i:length*(i+1)]
    	vv = v[length*i:length*(i+1)]
    
    	imagenu = uu.reshape(41,41)
    	imagenv = vv.reshape(41,41)
    
    	#fig = plt.figure()
    	fig = plt.figure(figsize=(11,7), dpi=100)
    	ax = fig.gca(projection='3d')
    	X,Y = np.meshgrid(x,y)
    	wire1 = ax.plot_wireframe(X,Y,imagenu)
    	wire2 = ax.plot_wireframe(X,Y,imagenv)
    	if (i<10):
    		plt.savefig("burger00"+str((i))+".png")
	if (i>=10 and i<100):
    		plt.savefig("burger0"+str((i))+".png")
	if (i>=100):
    		plt.savefig("burger"+str((i))+".png")
	
    	plt.close(fig)
    
os.system("convert -delay 20 -loop 0 burger*.png burger.gif")

