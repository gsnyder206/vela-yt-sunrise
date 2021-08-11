import matplotlib
matplotlib.use('Agg')  #apparently needed for old matplotlib/python on pleiades (didn't want to risk updating)
import matplotlib.pyplot as pyplot
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
import visualize_vela_hlsp as vvh

from PIL import Image
from astropy.convolution import *


def parse_vela_files(dirname='VELA01',genstr='v6'):

    #update final bibcode reference?

    #create manifest of all FITS files /vela/vela??/etc and catalog of basic parameters
    #link with galprops info?? a la luvoir sims

    #file, sim, redshift, camera, telescope, instrument, filter, pristine apparent mag, stellar mass? , halo mass?
    #add the new mag to the pristine header?

    image_files=np.sort(np.asarray(glob.glob(os.path.join('HLSP','vela',dirname.lower(),'*/*/*/*/*_'+genstr+'_*.fits'))))

    aux_files=np.sort(np.asarray(glob.glob(os.path.join('HLSP','vela',dirname.lower(),'*/*_'+genstr+'_*.fits'))))


    catfile=os.path.abspath(os.path.join('HLSP','vela','catalogs','hlsp_vela_multi_multi_'+dirname.lower()+'_multi'+'_'+genstr+'_cat.txt'))
    #                        new_filename='hlsp_vela_'+obs+'_'+instrument+'_'+dirname.lower()+'-'+cam+'-'+target_dir[14:-8]+'_'+fil.lower()+'_v3'+'_sim.fits'
    auxcatfile=os.path.abspath(os.path.join('HLSP','vela','catalogs','hlsp_vela_multi_multi_'+dirname.lower()+'_multi'+'_'+genstr+'_auxcat.txt'))

    print(catfile)
    print(auxcatfile)

    if genstr=='v6':
        gendir='Gen6'
    elif genstr=='v3-2':
        gendir='Gen3'

    propsdir=os.path.join('/nobackup/gfsnyder/VELA_sunrise/Runs/', gendir,dirname.lower())

    datf=os.path.join(propsdir,dirname+'_galprops.npy')
    dat=np.load(datf,encoding='bytes').all()


    aux_tfo=open(auxcatfile,'w')
    aux_tfo.write('sim z scale cam mstar mgas mmet sfr mvir_dm path\n')


    #sfr_tau = 1.2e7 #12 Myr, from Ceverino et al. 2015.  SFR values are actually SFR*Tau, so must divide by this factor to get true SFR values

    for auxfile in aux_files:
        print(auxfile)
        auxfo=pyfits.open(auxfile,mode='readonly')
        auxpri=auxfo[0]
        #if 'HLSPID' in auxpri.header:
        #    print('Aux file already updated, skipping..', auxfile)
        #    continue
        #auxpri.header['HLSPID']='vela'
        #auxpri.header['HLSPLEAD']='Gregory F. Snyder'
        #auxpri.header['HLSPNAME']='Vela-Sunrise Mock Images'
        #auxpri.header['HLSPVER']='v3'
        #auxpri.header['LICENSE']='CC BY 4.0'
        #auxpri.header['LICENURL']='https://creativecommons.org/licenses/by/4.0/'
        #auxpri.header['PROPOSID']='HST-AR#13887'
        #auxpri.header['REFERENC']='Simons et al. 2019'

        auxmain=auxfo[1]
        kpc_per_pix=auxmain.header['CD1_1']
        stellar_mass_per_pix=auxmain.data[4,:,:]*(kpc_per_pix**2)
        gas_mass_per_pix=auxmain.data[0,:,:]*(kpc_per_pix**2)
        metal_mass_per_pix=auxmain.data[1,:,:]*(kpc_per_pix**2)

        #sfrdens_per_pix=auxmain.data[2,:,:]/sfr_tau
        #auxmain.data[2,:,:]=sfrdens_per_pix

        sfr_per_pix=auxmain.data[2,:,:]*(kpc_per_pix**2) #sfrdens_per_pix*(kpc_per_pix**2)

        total_mstar=np.sum(stellar_mass_per_pix)
        total_mgas=np.sum(gas_mass_per_pix)
        total_mmet=np.sum(metal_mass_per_pix)
        total_sfr=np.sum(sfr_per_pix)

        sim=dirname.lower()
        scalestr=auxfile.split('_')[4].split('-')[-1][-5:]
        camstr=auxfile.split('_')[4].split('-')[1]

        scalefloat=float(scalestr)
        zfloat=(1.0/scalefloat) - 1.0


        #this is actually problematic because apparently the filenames aren't exact?
        mvirdm=dat[b'Mvir_dm'][dat[b'scale']==scalefloat]

        auxfo.flush()

        try:
            writestr='{:10s} {:10s} {:15.8f} {:10s} {:10s} {:15.6e} {:15.6e} {:15.6e} {:15.6e} {:15.6e} {:75s}\n'.format(gendir, sim,zfloat,scalestr,camstr,total_mstar,total_mgas,total_mmet,total_sfr,mvirdm[0],auxfile[5:])
            aux_tfo.write(writestr)
        except:
            print('weird problem with mvirdm?')

        aux_tfo.flush()
        sys.stdout.flush()


    aux_tfo.close()

    aux_table=ascii.read(auxcatfile)

    cat_tfo=open(catfile,'w')
    cat_tfo.write('gen sim z scalestr cam mission instrument filter flux_njy abmag mstar mgas mmet sfr mvir_dm path\n')

    for imfile in image_files:
        print(imfile[5:])

        #update HLSP reference
        fo=pyfits.open(imfile,mode='update')
        fo[0].header['REFERENC']='Simons et al. 2019'

        #measure and add pristine apparent magnitude?
        try:
            imhdu=fo['IMAGE_PRISTINE']
        except:
            fo.flush()
            continue

        imdata=imhdu.data
        assert(imhdu.header['IMUNIT']=='nanoJanskies')
        total_flux_njy=np.sum(imdata)
        imhdu.header['FLUX_NJY']=total_flux_njy
        mag_val=np.nan

        if total_flux_njy > 0.0:
            pristine_ab_apparent_mag=-2.5*np.log10((total_flux_njy*1.0e-9)/3631.0)
            imhdu.header['ABMAG']=(pristine_ab_apparent_mag,'pristine AB apparent mag')
            mag_val=pristine_ab_apparent_mag


        sim=dirname.lower()
        bn=os.path.basename(imfile)
        scalestr=bn.split('_')[4].split('-')[-1][-5:]
        camstr=bn.split('_')[4].split('-')[1]
        scalefloat=float(scalestr)

        #zfloat=(1.0/scalefloat)-1.0
        zfloat = fo['BROADBAND'].header['redshift']

        filtername=imhdu.header['FILTER']
        '''
        if filtername=='hst/wfc3_f336w' or filtername=='hst/acs_f435w' or filtername=='hst/acs_f814w' or filtername=='hst/wfc3_f160w' or filtername=='jwst/nircam_f200w' or filtername=='jwst/nircam_f444w' or filtername=='jwst/miri_f770w':
            if int(camstr[-2:]) < 7:
                ns_hdu=get_nonscatter_hdu(sim,scalestr,camstr,filtername)
                if ns_hdu is not None:
                    fo.append(ns_hdu)
        '''
        #fo.flush()


        mvirdm=dat[b'Mvir_dm'][dat[b'scale']==scalefloat]

        try:
            thing=mvirdm[0]
        except:
            continue

        aux_index=np.logical_and(aux_table['scale']==scalefloat, aux_table['cam']==camstr)
        assert(np.sum(aux_index)==1)
        total_mstar=aux_table['mstar'][aux_index].data[0]
        total_mgas=aux_table['mgas'][aux_index].data[0]
        total_mmet=aux_table['mmet'][aux_index].data[0]
        total_sfr=aux_table['sfr'][aux_index].data[0]
        aux_mvirdm=aux_table['mvir_dm'][aux_index].data[0]

        assert(np.abs((aux_mvirdm-mvirdm[0])/mvirdm[0])<1.0e-4)

        mission=bn.split('_')[2]
        instr=bn.split('_')[3]
        filname=bn.split('_')[5]

        #output catalog value
        catstr='{:10s} {:10s} {:12.8f} {:10s} {:10s} {:10s} {:10s} {:10s} {:15.6e} {:12.6f} {:15.6e} {:15.6e} {:15.6e} {:15.6e} {:15.6e} {:75s}\n'.format(gendir, sim,zfloat,scalestr,camstr,mission,instr,filname,
                                                                                                                                                 total_flux_njy,mag_val,total_mstar,total_mgas,total_mmet,total_sfr,mvirdm[0],imfile[5:])

        cat_tfo.write(catstr)
        cat_tfo.flush()
        sys.stdout.flush()


    cat_tfo.close()

    return


if __name__=="__main__":
    simname=sys.argv[1]
    genname=sys.argv[2]

    parse_vela_files(dirname=simname.upper(),genstr=genname)
