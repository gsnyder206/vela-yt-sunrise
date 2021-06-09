import vela_extract_images as vei


if __name__=="__main__":

    vei.do_single_snap(camlist=['cam00','cam01'],aux_only=False,genstr='v6',output_dir='/Users/gsnyder/Documents/VELATEST/HLSP',
                        psf_file='/Users/gsnyder/Dropbox/Projects/PythonCode/vela-yt-sunrise/kernels/gen6runs/vela_gen6_psfs.fits',
                        filter_dir='/Users/gsnyder/Dropbox/Projects/FILTERS/sunrise_filters/')
