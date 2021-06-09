import vela_extract_images as vei


if __name__=="__main__":

    vei.do_single_snap(genstr=sys.argv[1],output_dir='/nobackup/gfsnyder/VELA_sunrise/Outputs/HLSP',
                        psf_file='/u/gfsnyder/PythonCode/vela-yt-sunrise/kernels/gen6runs/vela_gen6_psfs.fits',
                        filter_dir='/u/gfsnyder/sunrise_data/sunrise_filters/')
