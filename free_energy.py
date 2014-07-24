#!/usr/bin/python

import numpy
import argparse
import matplotlib.pyplot as plt
import scipy.integrate


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



# Reads one pullf.xvg -file.
def read_xvg(pullf_file_name):
    
    # Open file
    pullf_file=open(pullf_file_name, 'r')
    # Init array
    forces=[]
    
    # Read file
    for line in pullf_file:
        # Skip comment/header lines 
        if (line.startswith('#') or  line.startswith('@')):
            continue

	# Add force to array
        timestep=float(line.split()[0])
        force=float(line.split()[1])
        forces.append([timestep,force])
    
    # Close file
    pullf_file.close()
     
    # create numpy-array
    forces_np = numpy.array(forces)

    # Return mean force
    return forces_np




# Reads list of z-coordinates and corresponding pullf.xvg -files.
# Returns z-zoordinates and mean forces.
def get_force_data(list_file_name):

    z=[] # init z-coord vector
    tf_list=[] # list of time-force functions

    # read list
    list_file=open(list_file_name, 'r')
    for line in list_file:
        pullf_file_name=line.split()[1]
	pullf_data=read_xvg(pullf_file_name)
        tf_list.append(pullf_data)
        z.append(float(line.split()[0]))
    list_file.close()

    # mirror axis
    z = numpy.multiply(z[::-1], -1)
    tf_list = tf_list[::-1]
    for tf in tf_list:
	tf[:,1] = numpy.multiply(tf[:,1], -1)

    return z, tf_list
    

    
# calculate deltaG 
def calc_deltaG(z,tf_list,first_tstep):
    
    f_mean=[]
    for tf in tf_list:
	# skip non-equiliberated time steps
	tf_equil = tf[ numpy.where( tf[:,0] >= first_tstep ) ]
	# calculate mean force
	f_mean.append(tf_equil[:,1].mean())
    
    # integrate f_mean as function of z
    deltaG=scipy.integrate.cumtrapz(f_mean,z,initial=0)

    return deltaG, f_mean




# run code
if __name__ == "__main__":

    # get data
    z, tf_list = get_force_data(files)
    
    # calculate deltaG
    deltaG, f_mean = calc_deltaG(z,tf_list,first_tstep)
    
    # plot
    plt.plot(z, f_mean, 'ro', z, deltaG, 'b-')
    plt.show()

