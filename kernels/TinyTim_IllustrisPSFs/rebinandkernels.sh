#!/bin/bash


python ~/Documents/PythonCode/rebin_image.py result00_psf_F160W.fits F160W_rebin.fits 132 30
python ~/Documents/PythonCode/rebin_image.py result00_psf_F125W.fits F125W_rebin.fits 132 30
python ~/Documents/PythonCode/rebin_image.py result00_psf_F105W.fits F105W_rebin.fits 132 30
python ~/Documents/PythonCode/rebin_image.py result00_psf_F140W.fits F140W_rebin.fits 132 30

python ~/Documents/PythonCode/rebin_image.py result00_psf_F275W.fits F275W_rebin.fits 268 50
python ~/Documents/PythonCode/rebin_image.py result00_psf_F336W.fits F336W_rebin.fits 268 100
python ~/Documents/PythonCode/rebin_image.py result00_psf_F435W.fits F435W_rebin.fits 268 110
python ~/Documents/PythonCode/rebin_image.py result00_psf_F814W.fits F814W_rebin.fits 268 105
python ~/Documents/PythonCode/rebin_image.py result00_psf_F775W.fits F775W_rebin.fits 268 100
python ~/Documents/PythonCode/rebin_image.py result00_psf_F606W.fits F606W_rebin.fits 268 95
python ~/Documents/PythonCode/rebin_image.py result00_psf_F850LP.fits F850LP_rebin.fits 268 90

