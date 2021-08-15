import os
import sys
import numpy as np
import glob
import shutil
import tarfile
import string
#import pandas
import astropy
import astropy.cosmology
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
import subprocess


if __name__=="__main__":

    #run from HLSP directory under simulation e.g. vela06

    dirname=os.path.basename(os.path.abspath('.'))  #e.g. vela06

    cams=np.sort(np.asarray(glob.glob('vela??_cam??.tar')))
    print(dirname)
    print(cams)

    for camst in cams:
        print('untarring ', camst)
        subprocess.call(['tar','xf',camst])
