#code to do source detection, segmentation, and morphology measurements from simple mock images
import astropy
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import numpy as np
import astropy.units as u
import glob
import sys
import os
import photutils
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
import scipy.ndimage
import scipy as sp
import statmorph


#Run basic source detection
def detect_sources(in_image,ext_name='MockImage_SB25',**kwargs):

    try:
        this_fo=fits.open(in_image,'append')
        this_hdu=this_fo[ext_name]
    except:
        print('HDU not found! ', in_image, ext_name)
        return
    
    rms=this_hdu.header['RMSNOISE']
    thr_rms = 1.3
    thresh = thr_rms * rms
    npixels = 10

    #build kernel for pre-filtering.  How big?
    #don't assume redshift knowledge here
    typical_kpc_per_arcsec = 8.0

    #abzp = image_hdu.header['ABZP']

    kernel_kpc_fwhm = 5.0
    kernel_arcsec_fwhm = kernel_kpc_fwhm/typical_kpc_per_arcsec
    kernel_pixel_fwhm = kernel_arcsec_fwhm/this_hdu.header['PIXSIZE']

    sigma = kernel_pixel_fwhm * gaussian_fwhm_to_sigma
    nsize = int(5*kernel_pixel_fwhm)
    kernel = Gaussian2DKernel(sigma, x_size=nsize, y_size=nsize)

    #may consider inputting an error image here -- can be computed with photutils plus a GAIN keyword -- ratio of flux units to counts

    #https://photutils.readthedocs.io/en/stable/segmentation.html
    
    #https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.detect_sources.html#photutils.segmentation.detect_sources
    segmap_obj = photutils.detect_sources(this_hdu.data, thresh, npixels=npixels, filter_kernel=kernel, **kwargs)
    segmap = segmap_obj.data

    #https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.source_properties.html#photutils.segmentation.source_properties
    props = photutils.source_properties(this_hdu.data, segmap)

    props_table=astropy.table.Table(props.to_table())
    #these give problems given their format/NoneType objects
    props_table.remove_columns(['sky_centroid','sky_centroid_icrs','source_sum_err','background_sum','background_mean','background_at_centroid'])
    
    #save segmap and info
    nhdu = fits.ImageHDU(np.int32(segmap))
    nhdu.header['EXTNAME']='SEGMAP'

    thdu = fits.BinTableHDU(props_table)
    thdu.header['EXTNAME']='SEGMAP_PROPS'
    
    this_fo.append(nhdu)
    this_fo.append(thdu)
    
    this_fo.flush()

    return segmap_obj, kernel


#Run PhotUtils Deblender
def deblend_sources(in_image,segm_obj,kernel,ext_name='MockImage_SB25',**kwargs):

    try:
        this_fo=fits.open(in_image,'append')
        this_hdu=this_fo[ext_name]
    except:
        print('HDU not found! ', in_image, ext_name)
        return


    #https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.deblend_sources.html#photutils.segmentation.deblend_sources
    segm_deblend = photutils.deblend_sources(this_hdu.data, segm_obj, npixels=10,
                                             filter_kernel=kernel, nlevels=32,
                                             contrast=0.001)
    segmap = segm_deblend.data

    #https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.source_properties.html#photutils.segmentation.source_properties    
    props = photutils.source_properties(this_hdu.data, segmap)

    props_table=props.to_table()
    #these give problems given their format/NoneType objects
    props_table.remove_columns(['sky_centroid','sky_centroid_icrs','source_sum_err','background_sum','background_mean','background_at_centroid'])
    
    #save segmap and info
    nhdu = fits.ImageHDU(np.int32(segmap))
    nhdu.header['EXTNAME']='DEBLEND'

    thdu = fits.BinTableHDU(props_table)
    thdu.header['EXTNAME']='DEBLEND_PROPS'
    
    
    this_fo.append(nhdu)
    this_fo.append(thdu)
    
    this_fo.flush()


    return


#Run morphology code
def run_statmorph(in_image,ext_name='MockImage_SB25',seg_name='DEBLEND'):

    #https://statmorph.readthedocs.io/en/latest/installation.html

    try:
        this_fo=fits.open(in_image,'append')
        this_hdu=this_fo[ext_name]
        segm_hdu=this_fo[seg_name]
    except:
        print('HDU not found! ', in_image, ext_name)
        return

    source_morph=statmorph.source_morphology(this_hdu.data,segm_hdu.data,gain=1000)

    return source_morph

