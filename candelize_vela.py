import cProfile
import pstats
import math
import string
import sys
import struct
import matplotlib
import numpy as np
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import scipy.odr as odr
import glob
import os
import gzip
import tarfile
import shutil
import congrid
import astropy.io.ascii as ascii
import warnings
import subprocess
import photutils
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import *
import astropy.io.fits as pyfits
import statmorph
import datetime
import setup_synthetic_images_mp as ssimp
import vela_extract_images as vei
import traceback
import visualize_vela_hlsp as vvh


def process_snapshot(subdirpath='.',mockimage_parameters=None,clobber=False, max=None, galaxy=None,seg_filter_label='NC-F200W',magsb_limits=[23.0,25.0,27.0,29.0],camindices='All',do_idl=False,analyze=True,use_nonscatter=True,Np=2,smc=False):

    cwd = os.path.abspath(os.curdir)

    os.chdir(subdirpath)

    if smc is True:
        smclab='smc'
    else:
        smclab=''
        
    bbfile_list = np.sort(np.asarray(glob.glob('broadbandz'+smclab+'.fits*')))   #enable reading .fits.gz files
    print(bbfile_list)

    if galaxy is not None:
        thisbb = np.where(bbfile_list==galaxy)[0]
        bbfile_list= bbfile_list[thisbb]

    test_file = bbfile_list[0]
    tf = pyfits.open(test_file)
    print(tf.info())
    print(tf['BROADBAND'].header.cards)
    print(tf['SFRHIST'].header.get('star_adaptive_smoothing'))
    print(tf['SFRHIST'].header.get('star_radius_factor'))


    if camindices=='All':
        N_cam=int(tf['MCRX'].header['N_CAMERA'])
        camindices=range(N_cam)
    
    #this is critical for later
    
    fils = tf['FILTERS'].data.field('filter')
    print(fils)



    '''
    filters_to_analyze = ['hst/acs_f435w','hst/acs_f606w','hst/acs_f775w','hst/acs_f850lp',
                          'hst/wfc3_f105w','hst/wfc3_f125w','hst/wfc3_f160w',
                          'jwst/nircam_f070w', 'jwst/nircam_f090w','jwst/nircam_f115w', 'jwst/nircam_f150w', 
                          'jwst/nircam_f200w', 'jwst/nircam_f277w', 'jwst/nircam_f356w', 'jwst/nircam_f444w', 
                          'hst/wfc3_f140w',
                          'hst/wfc3_f275w', 'hst/wfc3_f336w',
                          'hst/acs_f814w',
                          'jwst/miri_F560W','jwst/miri_F770W','jwst/miri_F1000W','jwst/miri_F1130W',
                          'jwst/miri_F1280W','jwst/miri_F1500W','jwst/miri_F1800W','jwst/miri_F2100W','jwst/miri_F2550W']


    '''

    filters_to_analyze=['hst/wfc3_f336w',
                        'hst/acs_f435w','hst/acs_f606w','hst/acs_f814w',
                        'hst/wfc3_f125w','hst/wfc3_f160w',
                        'wfirst/wfi_r062','wfirst/wfi_z087','wfirst/wfi_y106','wfirst/wfi_j129','wfirst/wfi_w146','wfirst/wfi_h158','wfirst/wfi_f184',
                        'jwst/nircam_f115w',
                        'jwst/nircam_f150w',
                        'jwst/nircam_f200w',
                        'jwst/nircam_f277w',
                        'jwst/nircam_f356w',
                        'jwst/nircam_f444w',
                        'jwst/miri_F770W',
                        'jwst/miri_F1500W']

    skip_filter_boolean=[False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False]

    
    
    print(filters_to_analyze)
    
    pixsize_arcsec = [0.03,
                      0.03,0.03,0.03,
                      0.06,0.06,
                      0.055,0.055,0.055,0.055,0.055,0.055,0.055,
                      0.032,0.032,0.032,
                      0.065,0.065,0.065,
                      0.11,0.11]

    filter_labels = ['WFC3-F336W',
                     'ACS-F435W','ACS-F606W','ACS-F814W',
                     'WFC3-F125W','WFC3-F160W',
                     'WFI-R062','WFI-Z087','WFI-Y106','WFI-J129','WFI-W146','WFI-H158','WFI-F184',
                     'NC-F115W','NC-F150W','NC-F200W',
                     'NC-F277W','NC-F356W','NC-F444W',
                     'MIRI-F770W','MIRI-F1500W']

    filter_indices = []

    print(len(filters_to_analyze), len(skip_filter_boolean), len(filter_labels))
    
    for i,f in enumerate(filters_to_analyze):
        fi = np.where(fils==f)
        print(fi[0][0], f, fils[fi[0][0]], filter_labels[i]) #, filters_to_analyze[fi]
        filter_indices.append(fi[0][0])


    filter_indices = np.asarray(filter_indices)

    print(filter_indices)

    #order of filter_labels in wavelength space (i.e, F435W is in the "2" position)
    filter_lambda_order=np.asarray([1,
                         2,3,5,
                         9,14,
                         4,6,7,10,11,13,15,
                         8,12,16,
                         17,18,19,
                         20,21])-1
    '''
    filter_lambda_order = [2,3,4,6,7,8,10,
                           11,12,13,14,15,16,17,18,
                           9,0,1,5,
                           19,20,21,22,
                           23,24,25,26,27]
    '''

    #photfnu units Jy; flux in 1 ct/s

    photfnu_Jy = [4.93e-7,
                  1.96e-7,9.17e-8,1.52e-7,
                  1.17e-7,1.52e-7,
                  3.17e-8,3.17e-8,3.17e-8,3.17e-8,3.17e-8,3.17e-8,3.17e-8,  #fix WFI photfnu to NIRCAM F115W -- similar instrument, middle of wl range, only matters for morph statistics
                  3.17e-8,2.68e-8,2.64e-8,
                  2.25e-8,2.57e-8,2.55e-8,
                  3.10e-8,4.48e-8]

    '''
    photfnu_Jy = [1.96e-7,9.17e-8,1.97e-7,4.14e-7,
                  1.13e-7,1.17e-7,1.52e-7,
                  5.09e-8,3.72e-8,3.17e-8,2.68e-8,2.64e-8,2.25e-8,2.57e-8,2.55e-8,
                  9.52e-8,8.08e-7,4.93e-7,1.52e-7,
                  5.75e-8,3.10e-8,4.21e-8,1.39e-7,
                  4.65e-8,4.48e-8,5.88e-8,4.98e-8,1.15e-7]
    '''
    
    morphcode_dir = "/Users/gsnyder/Documents/pro/morph_december2013/morph_pro/"
    morphcode_files = np.asarray(glob.glob(os.path.join(morphcode_dir,"*.*")))

    #se_dir = '/Users/gsnyder/Documents/Projects/Illustris_Morphology/Illustris-CANDELS/SE_scripts'
    #se_files = np.asarray(glob.glob(os.path.join(se_dir,"*.*")))

    psf_files = []
    psf_dir = os.path.expandvars('$GFS_PYTHON_CODE/vela-yt-sunrise/kernels')
    '''
    psf_names = ['TinyTim_IllustrisPSFs/F435W_rebin.fits','TinyTim_IllustrisPSFs/F606W_rebin.fits','TinyTim_IllustrisPSFs/F775W_rebin.fits','TinyTim_IllustrisPSFs/F850LP_rebin.fits',
                 'TinyTim_IllustrisPSFs/F105W_rebin.fits','TinyTim_IllustrisPSFs/F125W_rebin.fits','TinyTim_IllustrisPSFs/F160W_rebin.fits',
                 'WebbPSF_F070W_trunc.fits','WebbPSF_F090W_trunc.fits','WebbPSF_F115W_trunc.fits','WebbPSF_F150W_trunc.fits',
                 'WebbPSF_F200W_trunc.fits','WebbPSF_F277W_trunc.fits','WebbPSF_F356W_trunc.fits','WebbPSF_F444W_trunc.fits',
                 'TinyTim_IllustrisPSFs/F140W_rebin.fits','TinyTim_IllustrisPSFs/F275W_rebin.fits','TinyTim_IllustrisPSFs/F336W_rebin.fits','TinyTim_IllustrisPSFs/F814W_rebin.fits',
                 'WebbPSF_F560W_trunc.fits','WebbPSF_F770W_trunc.fits','WebbPSF_F1000W_trunc.fits','WebbPSF_F1130W_trunc.fits',
                 'WebbPSF_F1280W_trunc.fits','WebbPSF_F1500W_trunc.fits','WebbPSF_F1800W_trunc.fits','WebbPSF_F2100W_trunc.fits','WebbPSF_F2550W_trunc.fits']
    '''

    psf_names = ['TinyTim_IllustrisPSFs/F336W_rebin.fits',
                 'TinyTim_IllustrisPSFs/F435W_rebin.fits','TinyTim_IllustrisPSFs/F606W_rebin.fits','TinyTim_IllustrisPSFs/F814W_rebin.fits',
                 'TinyTim_IllustrisPSFs/F125W_rebin.fits','TinyTim_IllustrisPSFs/F160W_rebin.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_R062_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_Z087_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_Y106_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_J129_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_W149_SCA01_center_psf_phaseb.fits', #note filter name inconsistency with older version of WEBBPSF
                 'WebbPSF_WFIRST/PhaseB/WFI_H158_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_F184_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_F115W_trunc.fits','WebbPSF_F150W_trunc.fits','WebbPSF_F200W_trunc.fits',
                 'WebbPSF_F277W_trunc.fits','WebbPSF_F356W_trunc.fits','WebbPSF_F444W_trunc.fits',
                 'WebbPSF_F770W_trunc.fits','WebbPSF_F1500W_trunc.fits']
    
    #psf_pix_arcsec = [0.0125,0.0125,0.0125,0.0125,0.0325,0.0325,0.0325,0.007925,0.007925,0.007925,0.007925,0.007925,0.0162,0.0162,0.0162,0.0325,0.0100,0.0100,0.0125]
    #switch to JWST detector sampling for efficiency.  They're model psfs anyway, full accuracy not essential.  also JWST is better sampled than WFC3 or WFI -- smaller than det for HST/WFIRST is good

    psf_pix_arcsec = [0.03,
                      0.03,0.03,0.03,
                      0.06,0.06,
                      0.0275,0.0275,0.0275,0.0275,0.0275,0.0275,0.0275,
                      0.0317,0.0317,0.0317,
                      0.0648,0.0648,0.0648,
                      0.11,0.11]
    
    psf_truncate = [None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None]
    psf_hdu_num = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    
    psf_fwhm = [0.08,
                0.10,0.11,0.13,
                0.17,0.20,
                0.085*0.62,0.085*0.87,0.085*1.06,0.085*1.29,0.085*1.46,0.085*1.58,0.085*1.84,
                0.11,0.11,0.12,
                0.15,0.18,0.25,
                0.035*7.57,0.035*14.96]
    
    #these settings yield full subhalo (4 cams) convolution in 0.92s!  convolve_fft ftw!

    for pname in psf_names:
        psf_file = os.path.join(psf_dir,pname)
        psf_files.append(psf_file)
        print(psf_file, os.path.lexists(psf_file))

    ###  PSFSTD; WFC3 = 0.06 arcsec, ACS = 0.03 arcsec... I think
    ### NIRCAM in header with keyword 'PIXELSCL';  short 0.07925 long 0.0162
    ## acs wfc 0.05 arcsec pixels... PSFSTD x4 oversample?
    ## wfc3 ir 0.13 arcsec
    ## wfc3 uv 0.04 arcsec

    mockimage_parameters = ssimp.analysis_parameters('mockimage_default')
    mockimage_parameters.filter_indices = filter_indices
    mockimage_parameters.filter_labels = filter_labels
    mockimage_parameters.pixsize_arcsec = pixsize_arcsec
    mockimage_parameters.morphcode_base = morphcode_dir
    mockimage_parameters.morphcode_files = morphcode_files
    #mockimage_parameters.se_base = se_dir
    #mockimage_parameters.se_files = se_files
    mockimage_parameters.camera_indices = camindices #None #by default, do all
    mockimage_parameters.psf_files = psf_files
    mockimage_parameters.psf_pix_arcsec = psf_pix_arcsec
    mockimage_parameters.psf_truncate = psf_truncate
    mockimage_parameters.psf_hdu_num = psf_hdu_num
    mockimage_parameters.magsb_limits = magsb_limits
    mockimage_parameters.psf_fwhm_arcsec = psf_fwhm
    mockimage_parameters.photfnu_Jy = photfnu_Jy
    mockimage_parameters.filter_lambda_order = filter_lambda_order
    mockimage_parameters.skip_filters = skip_filter_boolean
    mockimage_parameters.use_nonscatter = use_nonscatter
    
    #use exactly one detection and segmentation per object, depending on redshift
    #enormous simplification
    #observationally, go w deepest filter.  here... ?

    mockimage_parameters.segment_filter_label = seg_filter_label
    mockimage_parameters.segment_filter_index = np.where(np.asarray(mockimage_parameters.filter_labels) == seg_filter_label)[0][0]

    print(mockimage_parameters.segment_filter_label)
    print(mockimage_parameters.segment_filter_index)
    
    assert(len(psf_pix_arcsec)==len(pixsize_arcsec))
    assert(len(filter_labels)==len(mockimage_parameters.psf_files))

    bbdirs = []
    
    for i,bbfile in enumerate(bbfile_list):

        try:
            bbdir = ssimp.process_single_broadband(bbfile,mockimage_parameters,clobber=clobber,do_idl=do_idl,analyze=analyze,bbase="broadbandz",Np=Np,zip_after=False,smc=smc)
            bbdirs.append(bbdir)
        except (KeyboardInterrupt,NameError,AttributeError,KeyError,TypeError,IndexError) as e:
            print(e)
            raise
        except:
            print("Exception while processing broadband: ", bbfile)
            print("Error:", sys.exc_info())
            print(traceback.print_tb(sys.exc_info()[2]))
        else:
            print("Successfully processed broadband: ", bbfile)

    os.chdir(cwd)

    return bbdirs




if __name__=="__main__":



    try:
        res1 = process_snapshot(subdirpath='.',clobber=False,seg_filter_label='NC-F200W',magsb_limits=[25.0,28.0],do_idl=False,analyze=False,use_nonscatter=False,smc=False,Np=4)
    except:
        print('failure in candelization of broadbandz.fits, scatter')
        continue
        
    try:
        res2 = process_snapshot(subdirpath='.',clobber=False,seg_filter_label='NC-F200W',magsb_limits=[25.0,28.0],do_idl=False,analyze=False,use_nonscatter=True,smc=False,Np=4)
    except:
        print('failure in candelization of broadbandz.fits, nonscatter')
        continue
        
    try:
        res3 = process_snapshot(subdirpath='.',clobber=False,seg_filter_label='NC-F200W',magsb_limits=[25.0,28.0],do_idl=False,analyze=False,use_nonscatter=False,smc=True,Np=4)
    except:
        print('failure in candelization of broadbandzsmc.fits, scatter')
        continue

    
        #extract hlsp stuff here

    vei.do_single_snap(genstr=sys.argv[1],output_dir='/nobackup/gfsnyder/VELA_sunrise/Outputs/HLSP')
    
    #tar and compress HLSP contents here?

    #now, tar ish up
    subprocess.call(['tar','cf',res1[0]+'.tar',res1[0]])
    subprocess.call(['rm','-rf',res1[0]])
    subprocess.call(['tar','cf',res2[0]+'.tar',res2[0]])
    subprocess.call(['rm','-rf',res2[0]])
    subprocess.call(['tar','cf',res3[0]+'.tar',res3[0]])
    subprocess.call(['rm','-rf',res3[0]])


    ss=os.path.abspath('.').split('/')[-3].lower()
    for i in range(25):
        camst='cam'+'{:02d}'.format(i)
        vvh.make_vela_stamps(bd='/nobackupp2/gfsnyder/VELA_sunrise/Outputs/HLSP/',sim=ss,cam=camst)
