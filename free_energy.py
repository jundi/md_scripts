#!/usr/bin/python

import numpy
import argparse
import matplotlib.pyplot as plt
from scipy import integrate


# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", help="list of windows and pullf.xvg -files")
parser.add_argument("-b", help="first time step to be used")
args = parser.parse_args()
if args.f:
    files = args.f
else:
    exit()
if args.b:
    first_tstep = float(args.b)
else:
    first_tstep = 0



# Reads one pullf.xvg -file and calculates mean force.
def mean_force(pullf_file_name):
    
    # Open file
    pullf_file=open(pullf_file_name, 'r')
    # Init array
    forces=numpy.array([])
    
    # Read file
    for line in pullf_file:
        # Skip comment/header lines 
        if (line.startswith('#') or  line.startswith('@')):
            continue

	# Add force to array
        timestep=float(line.split()[0])
        if (timestep >= first_tstep):
            force=float(line.split()[1])
            forces=numpy.append(forces,force)
    
    # Close file
    pullf_file.close()
     
    # Calculate mean force
    forces_mean=forces.mean()

    # Print something...
    print pullf_file_name + ' ' + 'time: ' + str(timestep) + ' ps, ' + str(forces.size) + ' frames, f=' + str(forces_mean)

    # Return mean force
    return forces_mean


# Reads list of z-coordinates and corresponding pullf.xvg -files.
# Returns z-zoordinates and mean forces.
def get_force(list_file_name):

    z=[] # init z-coord vector
    f=[] # init force vector

    # read list
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
    deltaG=integrate.cumtrapz(f,z,initial=0)
    return deltaG
    
    
# main function    
def pmf(list_file_name):
    
    # get forces
    data=get_force(list_file_name)
    zp=data[0]
    fp=data[1]
    
    # mirror z-axis to start integration from bulk water
    z=numpy.multiply(zp[::-1], -1)
    f=numpy.multiply(fp[::-1], -1)

    # integrate
    deltaG=integ_force(z,f)
    

    # plot
    plt.plot(z, f, 'ro', z, deltaG, 'b-')
    plt.show()


# run main function
if __name__ == "__main__":
  pmf(files)
