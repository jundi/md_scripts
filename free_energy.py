import numpy
import argparse
import matplotlib.pyplot as plt
from scipy import integrate


# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", help="list of windows and pullf.xvg -files")
args = parser.parse_args()
if args.f:
    files = args.f
else:
    exit()


# Reads one pullf.xvg -file and calculates mean force.
def mean_force(pullf_file_name):
    
    # Open file
    pullf_file=open(pullf_file_name, 'r')
    # Init array
    forces=numpy.array([])
    
    # Read file
    for line in pullf_file:
        # Skip comment/header lines
        if not (line.startswith('#') or line.startswith('@')):
            # Add force to array
            force=line.split()[1]
            forces=numpy.append(forces,float(force))
    
    # Close file
    pullf_file.close()
    # Return mean force
    return forces.mean()


# Reads list of z-coordinates and corresponding pullf.xvg -files.
# Returns z-zoordinates and mean forces.
def get_force(list_file_name):
    z=[]
    f=[]
    list_file=open(list_file_name, 'r')
    for line in list_file:
        pullf_file_name=line.split()[1]
        f.append(mean_force(pullf_file_name))
        z.append(float(line.split()[0]))
    list_file.close()
    return [z, f]
    

# Interpolate forces
def interp_force(zp,fp):
    step=0.01
    z=numpy.arange(min(zp), max(zp), step)
    f=numpy.interp(z, zp, fp)
    return [z, f]


# Integrate mean force
def integ_force(z,f):
    deltaG=integrate.cumtrapz(z,f,initial=0)
    return deltaG
    
    
# main function    
def pmf(list_file_name):
    
    # get forces
    zp=get_force(list_file_name)[0]
    fp=get_force(list_file_name)[1]
    
    # mirror z-axis to start integration from bulk water
    z=zp[::-1]
    f=fp[::-1]
    z=numpy.multiply(zp, -1)
    
    # integrate
    deltaG=integ_force(z,f)
    
    # plot
    plt.plot(z, f, 'ro', z, deltaG, 'b-')
    plt.show()


# run main function
pmf(files)