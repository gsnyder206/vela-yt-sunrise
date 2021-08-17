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

import vela_extract_images as vei


def retar_vela_files_by_filter_and_dust(genstr='v6',duststr='mw'):

    #run from HLSP/vela/vela??/ subdirectory

    dirname=os.path.basename(os.path.abspath('.'))

    #hlsp_sim_dir=os.path.join('HLSP','vela',dirname.lower())
    #print(hlsp_sim_dir)

    #cwd=os.path.abspath(os.curdir)
    #os.chdir(hlsp_sim_dir)

    print('Tarring.. ',  dirname, genstr, duststr)

    obslist=['hst','jwst','roman','aux']
    for obs in obslist:
        for instrument in vei.insdict[obs]:
            for filname in vei.fildict[instrument]:
                lfil=filname.lower()
                print(obs,instrument,lfil)
                if lfil=='aux':
                    imagefiles=np.sort(np.asarray(glob.glob('cam*/*_'+genstr+'_sim.fits')))
                    tarfilename='hlsp_vela_none_none_'+dirname.lower()+'_'+lfil+'_'+genstr+'_sim.tar'
                else:
                    imagefiles=np.sort(np.asarray(glob.glob('cam*/'+obs+'/'+instrument+'/'+lfil+'/*_'+genstr+'_sim-'+duststr+'.fits')))
                    tarfilename='hlsp_vela_'+obs+'_'+instrument+'_'+dirname.lower()+'_'+lfil+'_'+genstr+'_sim-'+duststr+'.tar'

                print(tarfilename)
                print(imagefiles[0:25])
                print(imagefiles.shape)

                tfo=tarfile.open(tarfilename,mode='a')
                i=0
                for imf in imagefiles:
                    i=i+1
                    print(imf)
                    tfo.add(imf)
                    if i % 50==0:
                        sys.stdout.flush()

                tfo.close()


    return


if __name__=="__main__":
    retar_vela_files_by_filter_and_dust(genstr='v6',duststr='mw')
    retar_vela_files_by_filter_and_dust(genstr='v6',duststr='ns')
    retar_vela_files_by_filter_and_dust(genstr='v6',duststr='smc')

    try:
        retar_vela_files_by_filter_and_dust(genstr='v3-2',duststr='mw')
        retar_vela_files_by_filter_and_dust(genstr='v3-2',duststr='ns')
        retar_vela_files_by_filter_and_dust(genstr='v3-2',duststr='smc')
    except:
        pass
