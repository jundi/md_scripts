#!/usr/bin/python

import numpy
import argparse
import matplotlib.pyplot as plt
import scipy.integrate
import progressbar as pb


# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", help="list of windows and pullf.xvg -files")
parser.add_argument("-b", help="first time step to be used")
parser.add_argument("-e", help="last time step to be used")
parser.add_argument("-dt", help="use frame only when t MOD dt = first time (ps)")
args = parser.parse_args()
if args.f:
    files = args.f
else:
    if __name__ == "__main__":
        exit()
if args.b:
    first_tstep = float(args.b)
else:
    first_tstep = 0
if args.e:
    last_tstep = float(args.e)
else:
    last_tstep = 0
if args.dt:
    dt = float(args.dt)
else:
    dt = 1



def read_xvg(pullf_file_name):
    '''Reads one pullf.xvg -file.'''
    
    # Open file
    pullf_file=open(pullf_file_name, 'r')
    # Init array
    forces=[]
    

    # Read file
    for line in pullf_file:
        # Skip comment/header lines 
        if (line.startswith('#') or  line.startswith('@')):
            continue

        timestep=float(line.split()[0])

        # Add force to array
        force=float(line.split()[1])
        forces.append([timestep,force])
    
    # Close file
    pullf_file.close()
     
    # Return mean force
    return forces





def get_force_data(list_file_name):
    '''Reads list of z-coordinates and corresponding pullf.xvg -files.  Returns
    list of z-zoordinates and list of numpy-arrays of pull-forces as function
    of time.'''

    z=[] # init z-coord vector
    tf_list=[] # list of time-force functions

    # read list
    list_file=open(list_file_name, 'r')
    for line in list_file:
        pullf_file_name=line.split()[1]
        print('reading ' + pullf_file_name + '\r', end='')
        pullf_data = read_xvg(pullf_file_name)
        tf_list.append(numpy.array(pullf_data))
        z.append(float(line.split()[0]))
    list_file.close()

    # mirror axis
    z = numpy.multiply(z[::-1], -1)
    tf_list = tf_list[::-1]
    for tf in tf_list:
        tf[:,1] = numpy.multiply(tf[:,1], -1)

    return z, tf_list
    


def skip_frames(tf_list, first_tstep, last_tstep, dt):
    '''Skips frames before 'first_tstep'.'''

    new_tf_list = []
    for tf in tf_list:
        # skip non-equiliberated time steps
        new_tf_list.append(tf[ numpy.where( tf[:,0] >= first_tstep ) ])


    # Remove frames after last_tstep
    if last_tstep > 0:
        tf_list_last = []
        for tf in new_tf_list:
            tf_list_last.append(tf[ numpy.where( tf[:,0] <= last_tstep ) ])
        new_tf_list=tf_list_last


    # Only write frame when t MOD dt = first time (ps)
    if dt > 0:
        tf_list_dt = []
        dt = round(1000*float(dt))
        for tf in new_tf_list:
            tf_list_dt.append(tf[ numpy.where( (tf[:,0]*1000).round() % dt == 0) ])
        new_tf_list = tf_list_dt


    return new_tf_list



def calc_mean_forces(tf_list):
    '''Calculates mean force'''

    f_mean=[]
    for tf in tf_list:
        # calculate mean force
        f_mean.append(tf[:,1].mean())

    return numpy.array(f_mean)


    
def calc_deltaG(z, tf_list, first_tstep, last_tstep, dt):
    '''Calculates deltaG'''
    
    # skip frames before 'first_tstep'
    tf_list=skip_frames(tf_list, first_tstep, last_tstep, dt)

    # get mean force -vector
    f_mean = calc_mean_forces(tf_list)

    # integrate f_mean as function of z
    deltaG = scipy.integrate.cumtrapz(f_mean,z,initial=0)

    return deltaG, f_mean




def autocor(x):
    '''Calculates autocorrelation function of x'''
    result = numpy.correlate(x, x, mode='full') / len(x)
    return result[result.size/2:]




def calc_diffcoef(z_list,tf_list,first_tstep, last_tstep, dt):
    '''Calculates diffusion coefficient D'''
    
    # skip frames before 'first_tstep'
    tf_list=skip_frames(tf_list, first_tstep, last_tstep, dt)

    # mean forces
    f_mean_list = calc_mean_forces(tf_list)

    # random forces deltaF 
    deltaF_list = []
    for tf,f_mean in zip(tf_list, f_mean_list):
        deltaF_list.append( tf[:,1] - f_mean )

    # time
    time_list = []
    for tf in tf_list:
        time_list.append(tf[:,0])

    # Autocorrelation function <deltaF(z,t),deltaF(z,0)>
    acor_list=[]
    for deltaF,z in zip(deltaF_list,z_list):
        line = 'calculating autocorrelation functions ' + str(z)
        print(line + '\r', end='' )

        acor_list.append( autocor(deltaF) )

        print(' '*len(line) + '\r', end='')



    # Diffusion coefficient
    D_list = []
    for acor,time,z in zip(acor_list,time_list,z_list):
        line = 'calculating diffusion coefficient ' + str(z)
        print(line + '\r', end='' )

        D  = 1/numpy.trapz(acor,time) # some random units...
        D_list.append(D)

        print(' '*len(line) + '\r', end='')
    
    return D_list



# run code
if __name__ == "__main__":

    # get data
    z, tf_list = get_force_data(files)

    # get data from npy
    #z, tf_list = get_force_data_npy(files)

    # save data for reuse
    #numpy.save(files+'.npy',tf_list)

    # calculate deltaG
    deltaG, f_mean = calc_deltaG(z, tf_list, first_tstep, last_tstep, dt)

    # save deltaG
    dG_plot = numpy.array([z, deltaG])
    numpy.savetxt('deltaG', dG_plot.transpose())
    
    # plot dG
    #plt.plot(z, f_mean, 'ro', z, deltaG, 'b-')
    #plt.show()

    # calculate diffusion coefficient
    #D=calc_diffcoef(z, tf_list, first_tstep, last_tstep, dt)
    #print(D)

     #plot diffcoeff.
    #plt.plot(z, D)
    #plt.show()

    #numpy.savetxt('tf_list', tf_list[1])
    #numpy.savetxt('D', D[1])
    #numpy.savetxt('f_mean', f_mean)
