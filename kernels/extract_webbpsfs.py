import math
import string
import sys
import struct
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.colors as pycolors
import matplotlib.cm as cm
import matplotlib.patches as patches
import numpy as np
import cPickle
import asciitable
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import scipy.odr as odr
import astropy.io.fits as pyfits
import glob
import os
import make_color_image
import make_fake_wht
import gzip
import tarfile
import shutil
import cosmocalc
import congrid
import astropy.io.ascii as ascii


if __name__=="__main__":
    ncpsfs = ['PSF_NIRCam_F070W_revV-1.fits','PSF_NIRCam_F090W_revV-1.fits','PSF_NIRCam_F115W_revV-1.fits','PSF_NIRCam_F150W_revV-1.fits',
              'PSF_NIRCam_F200W_revV-1.fits','PSF_NIRCam_F277W_revV-1.fits','PSF_NIRCam_F356W_revV-1.fits','PSF_NIRCam_F444W_revV-1.fits',
              'PSF_MIRI_F560W_revV-1.fits','PSF_MIRI_F770W_revV-1.fits','PSF_MIRI_F1000W_revV-1.fits','PSF_MIRI_F1130W_revV-1.fits',
              'PSF_MIRI_F1280W_revV-1.fits','PSF_MIRI_F1500W_revV-1.fits','PSF_MIRI_F1800W_revV-1.fits','PSF_MIRI_F2100W_revV-1.fits','PSF_MIRI_F2550W_revV-1.fits']

    outpsfs = ['WebbPSF_F070W_trunc.fits','WebbPSF_F090W_trunc.fits','WebbPSF_F115W_trunc.fits','WebbPSF_F150W_trunc.fits',
               'WebbPSF_F200W_trunc.fits','WebbPSF_F277W_trunc.fits','WebbPSF_F356W_trunc.fits','WebbPSF_F444W_trunc.fits',
               'WebbPSF_F560W_trunc.fits','WebbPSF_F770W_trunc.fits','WebbPSF_F1000W_trunc.fits','WebbPSF_F1130W_trunc.fits',
               'WebbPSF_F1280W_trunc.fits','WebbPSF_F1500W_trunc.fits','WebbPSF_F1800W_trunc.fits','WebbPSF_F2100W_trunc.fits','WebbPSF_F2550W_trunc.fits',]

    truncs = [40,40,40,40,40,30,30,30,
              40,40,50,50,60,60,80,80,80]

    for i,fits in enumerate(ncpsfs):
        hdulist = pyfits.open(fits)
        det_samp = hdulist['DET_SAMP'].data
        psfc = det_samp.shape[0]/2
        st = truncs[i] #limit to the center 80 pixels
        det_samp_center = det_samp[psfc-st:psfc+st,psfc-st:psfc+st]
        
        orig_header = hdulist['DET_SAMP'].header

        newhdu = pyfits.PrimaryHDU(det_samp_center, header=orig_header)
        newlist = pyfits.HDUList([newhdu])
        newlist.writeto(outpsfs[i], clobber=True)
        res = pyfits.open(outpsfs[i])
        res.info()

