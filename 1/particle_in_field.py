from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot
import numpy
import sys
from matplotlib.backends.backend_pdf import PdfPages

array = loadtxt(sys.argv[0]))

x = array[:,0]
y = array[:,1]
z = array[:,2]
t = array[:,3]

with PdfPages('particle_in_field.pdf') as pdf:
    
    plt.plot(x,y)
    plt.title('Trayectoria xy')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    fig = pyplot.figure(figsize=(11,7), dpi=100)
    ax = fig.gca(projection='3d')
    
    ax.plot(x,y,z)

    plt.title('Trayectoria xyz')
    
    pdf.savefig()
    plt.close()
