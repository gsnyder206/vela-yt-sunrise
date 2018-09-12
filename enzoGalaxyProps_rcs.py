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
    if len(sys.argv)==3:
        snaps = np.asarray([sys.argv[1]])+'/' +np.asarray([sys.argv[1]])
        form= sys.argv[2]
    else:
        snaps = np.ssort(np.asarray(glob.glob("RD????/RD????")))  #ENZO format a list of snapshots in separate directories
        form='ENZO'


    print snaps, form







