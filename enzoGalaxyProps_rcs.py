import sys
import os
import glob
import yt
import numpy as np
from numpy import *
import astropy
from astropy.cosmology import Planck13 as cosmo
import findGalaxyProps as fGP


if __name__=="__main__":
    snaps = np.sort(np.asarray(glob.glob("RD????/RD????")))  #ENZO format a list of snapshots in separate directories
    form='ENZO'

    assert snaps.shape[0] > 0

    print("Calculating Galaxy Props for "+form+": ", snaps)

    abssnap = os.path.abspath(snaps[0])
    dirname = os.path.dirname(os.path.dirname(abssnap))
    simname = os.path.basename(dirname) #assumes directory name for simulation name

    print( "Simulation name:  ", simname)

    particle_headers = []
    particle_data = []
    stars_data = []
    new_snapfiles = []
    for sn in snaps:
        aname=os.path.basename(sn)
        adir=os.path.abspath(os.path.dirname(sn))
        snap_dir = os.path.join(adir,simname+'_'+aname+'_sunrise')
        yt_fig_dir = snap_dir+'/yt_projections'
        print( "Sunrise directory: ", snap_dir)
        if not os.path.lexists(snap_dir):
            print ("Creating Sunrise directory:", snap_dir)
            os.mkdir(snap_dir)        
        if not os.path.lexists(yt_fig_dir):
            print ("Creating YT figure directory:", yt_fig_dir)
            os.mkdir(yt_fig_dir)        




