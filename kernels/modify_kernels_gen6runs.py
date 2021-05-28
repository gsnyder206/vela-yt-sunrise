import astropy
import astropy.io.fits as fits

import os
import congrid

from PIL import Image
import numpy as np


psf_names = ['TinyTim_IllustrisPSFs/F336W_rebin.fits',
                 'TinyTim_IllustrisPSFs/F435W_rebin.fits','TinyTim_IllustrisPSFs/F606W_rebin.fits','TinyTim_IllustrisPSFs/F814W_rebin.fits',
                 'TinyTim_IllustrisPSFs/F125W_rebin.fits','TinyTim_IllustrisPSFs/F160W_rebin.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_R062_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_Z087_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_Y106_SCA01_center_psf_phaseb.fits',
                 #'WebbPSF_WFIRST/PhaseB/WFI_J129_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_W149_SCA01_center_psf_phaseb.fits', #note filter name inconsistency with older version of WEBBPSF
                 #'WebbPSF_WFIRST/PhaseB/WFI_H158_SCA01_center_psf_phaseb.fits',
                 'WebbPSF_WFIRST/PhaseB/WFI_F184_SCA01_center_psf_phaseb.fits',
                 'oversamp_WebbPSF_F115W_trunc.fits','oversamp_WebbPSF_F150W_trunc.fits','oversamp_WebbPSF_F200W_trunc.fits',
                 'oversamp_WebbPSF_F277W_trunc.fits','oversamp_WebbPSF_F356W_trunc.fits','oversamp_WebbPSF_F444W_trunc.fits',
                 'oversamp_WebbPSF_F770W_trunc.fits','oversamp_WebbPSF_F1500W_trunc.fits']

target_pix_sizes=[0.03,
                    0.03,0.03,0.03,
                    0.06,0.06,
                    0.06,0.06,0.06,0.06,0.06, #revisit?
                    0.03,0.03,0.03,
                    0.03,0.03,0.03,
                    0.12,0.12]


out_dir='gen6runs'
if not os.path.lexists(out_dir):
    os.makedirs(out_dir)

for pf,new_size in zip(psf_names,target_pix_sizes):
    out_name=os.path.basename(pf)
    print(out_name,new_size)
    fo=fits.open(pf)
    #print(fo.info())
    try:
        pix_arcsec=fo[0].header['PIXSCALE']
    except:
        pix_arcsec=fo[0].header['PIXELSCL']
    psf_use = fo[0].data
    np_orig = psf_use.shape[0]
    np_new = np_orig*pix_arcsec/new_size

    np_new_int=np.int64(np_new)
    orig_fov = np_orig*pix_arcsec
    new_fov = np_new_int*new_size
    diff = (orig_fov-new_fov)/2.0

    box_arcsec=(diff,diff,orig_fov-diff,orig_fov-diff)
    box=(diff/pix_arcsec,diff/pix_arcsec,(orig_fov-diff)/pix_arcsec,(orig_fov-diff)/pix_arcsec)

    print(pix_arcsec, np_orig, np_new,box)

    new_psf = Image.fromarray(psf_use).resize(size=(np_new_int, np_new_int),box=box)
    new_hdu = fits.PrimaryHDU(new_psf)

    new_hdu.header['PIXSCALE']=(new_size,'arcsec')
    new_hdu.header['ORIGSCL']=(pix_arcsec,'original PSF pixel size')


    print(new_hdu.data.shape)
    newfo=fits.HDUList(new_hdu)
    newfo.writeto(os.path.join(out_dir,out_name),overwrite=True)

    fo.close()
