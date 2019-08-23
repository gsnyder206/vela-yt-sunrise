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
#import subprocess

'''
['hst/wfc3_f275w', 'hst/wfc3_f336w', 'hst/acs_f435w',
       'hst/acs_f606w', 'hst/acs_f775w', 'hst/acs_f814w', 'hst/acs_f850lp',
       'hst/wfc3_f105w', 'hst/wfc3_f125w', 'hst/wfc3_f140w',
       'hst/wfc3_f160w', 'WFI_DRM15_Z087', 'WFI_DRM15_Y106',
       'WFI_DRM15_J129', 'WFI_DRM15_W149', 'WFI_DRM15_H158',
       'WFI_DRM15_F184', 'jwst/nircam_f070w', 'jwst/nircam_f090w',
       'jwst/nircam_f115w', 'jwst/nircam_f150w', 'jwst/nircam_f200w',
       'jwst/nircam_f277w', 'jwst/nircam_f356w', 'jwst/nircam_f444w',
       'jwst/miri_F560W', 'jwst/miri_F770W', 'jwst/miri_F1000W',
       'jwst/miri_F1130W', 'jwst/miri_F1280W', 'jwst/miri_F1500W',
       'jwst/miri_F1800W', 'jwst/miri_F2100W', 'jwst/miri_F2550W']
'''
'''
insdict={'hst':['wfc3','acs'],
         'jwst':['nircam','miri'],
         'wfirst':['wfidrm15'],
         'aux':['aux']}

cams=['cam00','cam01','cam02','cam03','cam04','cam05','cam06','cam07','cam08','cam09',
      'cam10','cam11','cam12','cam13','cam14','cam15','cam16','cam17','cam18']

fildict={'wfc3':['F275W','F336W','F105W','F125W','F140W','F160W'],
         'acs':['F435W','F606W','F775W','F814W','F850LP'],
         'nircam':['F070W','F090W','F115W','F150W','F200W','F277W','F356W','F444W'],
         'miri':['F560W','F770W','F1000W','F1130W','F1280W','F1500W','F1800W','F2100W','F2550W'],
         'wfidrm15':['Z087','Y106','J129','W149','H158','F184'],
         'aux':['aux']}

filfil={'F275W':'hst/wfc3_f275w', 'F336W':'hst/wfc3_f336w', 'F435W':'hst/acs_f435w',
       'F606W':'hst/acs_f606w', 'F775W':'hst/acs_f775w', 'F814W':'hst/acs_f814w', 'F850LP':'hst/acs_f850lp',
       'F105W':'hst/wfc3_f105w', 'F125W':'hst/wfc3_f125w', 'F140W':'hst/wfc3_f140w',
       'F160W':'hst/wfc3_f160w', 'Z087':'WFI_DRM15_Z087', 'Y106':'WFI_DRM15_Y106',
       'J129':'WFI_DRM15_J129', 'W149':'WFI_DRM15_W149', 'H158':'WFI_DRM15_H158',
       'F184':'WFI_DRM15_F184', 'F070W':'jwst/nircam_f070w', 'F090W':'jwst/nircam_f090w',
       'F115W':'jwst/nircam_f115w', 'F150W':'jwst/nircam_f150w', 'F200W':'jwst/nircam_f200w',
       'F277W':'jwst/nircam_f277w', 'F356W':'jwst/nircam_f356w', 'F444W':'jwst/nircam_f444w',
       'F560W':'jwst/miri_F560W', 'F770W':'jwst/miri_F770W', 'F1000W':'jwst/miri_F1000W',
       'F1130W':'jwst/miri_F1130W', 'F1280W':'jwst/miri_F1280W', 'F1500W':'jwst/miri_F1500W',
        'F1800W':'jwst/miri_F1800W', 'F2100W':'jwst/miri_F2100W', 'F2550W':'jwst/miri_F2550W','aux':'aux'}
'''

insdict={'hst':['wfc3','acs'],
         'jwst':['nircam','miri'],
         'wfirst':['wfi'],
         'aux':['aux']}

cams=['cam00','cam01','cam02','cam03','cam04','cam05','cam06','cam07','cam08','cam09',
      'cam10','cam11','cam12','cam13','cam14','cam15','cam16','cam17','cam18','cam19',
      'cam20','cam21','cam22','cam23','cam24']

fildict={'wfc3':['F336W','F125W','F160W'],
         'acs':['F435W','F606W','F814W'],
         'nircam':['F115W','F150W','F200W','F277W','F356W','F444W'],
         'miri':['F770W','F1500W'],
         'wfi':['Z087','Y106','J129','W146','H158','F184'],
         'aux':['aux']}

filfil={'F336W':'hst/wfc3_f336w', 'F435W':'hst/acs_f435w',
        'F606W':'hst/acs_f606w', 'F814W':'hst/acs_f814w',
        'F125W':'hst/wfc3_f125w', 'F160W':'hst/wfc3_f160w',
        'Z087':'wfirst/wfi_z087', 'Y106':'wfirst/wfi_y106',
        'J129':'wfirst/wfi_j129', 'W146':'wfirst/wfi_w146', 'H158':'wfirst/wfi_h158',
        'F184':'wfirst/wfi_f184','R062':'wfirst/wfi_r062',
        'F115W':'jwst/nircam_f115w', 'F150W':'jwst/nircam_f150w', 'F200W':'jwst/nircam_f200w',
        'F277W':'jwst/nircam_f277w', 'F356W':'jwst/nircam_f356w', 'F444W':'jwst/nircam_f444w',
        'F770W':'jwst/miri_F770W', 
        'F1500W':'jwst/miri_F1500W','aux':'aux'}


ilh = 0.704
illcos = astropy.cosmology.FlatLambdaCDM(H0=70.4,Om0=0.2726,Ob0=0.0456)


def vela_export_image(hdulist,camnum,filtername,label='',nonscatter=False):

    #get image in W/m/m^2/Sr

    #hdulist=fits.open(bbfile)
    fils=hdulist['FILTERS'].data['filter']
    fi=np.where(fils==filtername)[0][0]
    
    efl=hdulist['FILTERS'].data['lambda_eff']
    
    efl_microns=1.0e6 * efl[fi]

    if not np.isfinite(efl_microns):
        return None
    
    if nonscatter is False:
        key='CAMERA'+'{}'.format(camnum)+'-BROADBAND'
    else:
        key='CAMERA'+'{}'.format(camnum)+'-BROADBAND-NONSCATTER'

        

    camdata=hdulist[key].data[fi,:,:]
    
    #convert to nJy
    redshift=hdulist['BROADBAND'].header['redshift']
    pixsize_kpc=hdulist[key].header['CD1_1']
    pixsize_arcsec=pixsize_kpc/(illcos.kpc_proper_per_arcmin(redshift).value/60.0)
    
    sq_arcsec_per_sr = 42545170296.0
    c = 3.0e8

    pixel_Sr = (pixsize_arcsec**2)/sq_arcsec_per_sr  #pixel area in steradians:  Sr/pixel
    to_nJy_per_Sr = (1.0e9)*(1.0e14)*(efl_microns**2)/c   #((pixscale/206265.0)^2)*
    to_nJy_per_pix = to_nJy_per_Sr*pixel_Sr

    newdata=camdata*to_nJy_per_pix

    outhdu=pyfits.ImageHDU(newdata)

    outhdu.header['redshift']=hdulist['BROADBAND'].header['redshift']
    outhdu.header['PIXSIZE']=(pixsize_arcsec,'arcsec')
    outhdu.header['PIXKPC']=(pixsize_kpc,'kpc')
    outhdu.header['IMUNIT']=('nanoJanskies')
    outhdu.header['SBFACTOR']=(to_nJy_per_pix,'Image=SBFACTOR*Original')
    outhdu.header['EFLAMBDA']=(efl_microns,'microns')
    outhdu.header['CODE']=('Sunrise (Jonsson 06)')
    if nonscatter is False:
        outhdu.header['EXTNAME']='IMAGE_PRISTINE'
    else:
        outhdu.header['EXTNAME']='IMAGE_PRISTINE_NONSCATTER'

    outhdu.header['FILTER']=filtername
    outhdu.header['CAMERA']=(camnum,'Sunrise camera number')
    
    return outhdu


#utility to help with integrating pleiades pipeline
def do_single_snap(obslist=['hst','jwst','wfirst'],camlist=cams, aux_only=False, output_dir='/nobackup/gfsnyder/VELA_sunrise/Outputs/HLSP',genstr='v3-2'):

    #assume run within /images/ subdirectory?

    imagedir=glob.glob('images_*_sunrise_mw')[0]
    imagedir_ns=glob.glob('images_*_sunrise_nonscatter')[0]
    imagedir_smc=glob.glob('images_*_sunrise_smc')[0]
    
    bb_fits='broadbandz.fits'
    bb_fits_smc='broadbandzsmc.fits'

    hdulist=pyfits.open(bb_fits)
    hdulist_smc=pyfits.open(bb_fits_smc)

    target_dir=os.path.basename(imagedir)

    dirname=imagedir.split('_')[1]
    
    if True:
        for obs in obslist:
            for instrument in insdict[obs]:
                if instrument is 'wfc3':
                    instrumentfind='WFC3'
                elif instrument is 'acs':
                    instrumentfind='ACS'
                elif instrument is 'nircam':
                    instrumentfind='NC'
                elif instrument is 'miri':
                    instrumentfind='MIRI'
                elif instrument is 'wfi':
                    instrumentfind='WFI'
                for cam in camlist:
                    auxdir=os.path.join(output_dir,'vela',dirname.lower(),cam)
                    if not os.path.lexists(auxdir):
                        os.makedirs(auxdir)
                    auxhdu=hdulist['CAMERA'+str(int(cam[-2:]))+'-AUX']
                    auxoutfile=os.path.join(auxdir,'hlsp_vela_none_none_'+dirname.lower()+'-'+cam+'-'+target_dir[14:-8]+'_aux_'+genstr+'_sim.fits')

                    if not os.path.lexists(auxoutfile):
                        auxhdu.writeto(auxoutfile,overwrite=True)

                        auxfo=pyfits.open(auxoutfile,mode='update')
                        auxpri=auxfo[0]
                        auxpri.header['HLSPID']='vela'
                        auxpri.header['HLSPLEAD']='Gregory F. Snyder'
                        auxpri.header['HLSPNAME']='Vela-Sunrise Mock Images'
                        auxpri.header['HLSPVER']=genstr
                        auxpri.header['HLSPDOI']='10.17909/t9-ge0b-jm58'
                        auxpri.header['LICENSE']='CC BY 4.0'
                        auxpri.header['LICENURL']='https://creativecommons.org/licenses/by/4.0/'
                        auxpri.header['PROPOSID']='HST-AR#13887'
                        auxpri.header['REFERENC']='Simons et al. 2019'
                        auxpri.header['REFDOI']='10.3847/1538-4357/ab07c9'

                        sfr_tau = 1.2e7 #12 Myr, from Ceverino et al. 2015.  SFR values are actually SFR*Tau, so must divide by this factor to get true SFR values

                    
                        auxmain=auxfo[1]
                        kpc_per_pix=auxmain.header['CD1_1']
                        stellar_mass_per_pix=auxmain.data[4,:,:]*(kpc_per_pix**2)
                        gas_mass_per_pix=auxmain.data[0,:,:]*(kpc_per_pix**2)
                        metal_mass_per_pix=auxmain.data[1,:,:]*(kpc_per_pix**2)
                        
                        sfrdens_per_pix=auxmain.data[2,:,:]/sfr_tau
                        auxmain.data[2,:,:]=sfrdens_per_pix
                        
                        #sfr_per_pix=sfrdens_per_pix*(kpc_per_pix**2)
                        
                        #total_mstar=np.sum(stellar_mass_per_pix)
                        #total_mgas=np.sum(gas_mass_per_pix)
                        #total_mmet=np.sum(metal_mass_per_pix)
                        #total_sfr=np.sum(sfr_per_pix)
                        
                        #sim=dirname.lower()
                        #scalestr=auxfile.split('_')[4].split('-')[-1][-5:]
                        #camstr=auxfile.split('_')[4].split('-')[1]
                        #scalefloat=float(scalestr)
                        #zfloat=(1.0/scalefloat) - 1.0
                        
                        
                        #this is actually problematic because apparently the filenames aren't exact?
                        #mvirdm=dat[b'Mvir_dm'][dat[b'scale']==scalefloat]
                    
                        auxfo.flush()
                    else:
                        print(auxoutfile, ' exists, skipping..')

                        
                    for fil in fildict[instrument]:
                        if fil=='aux' or aux_only is True:
                            continue
                            
                        outdir=os.path.join(output_dir,'vela',dirname.lower(),cam,obs,instrument,fil.lower(),'')



                        new_filename='hlsp_vela_'+obs+'_'+instrument+'_'+dirname.lower()+'-'+cam+'-'+target_dir[14:-8]+'_'+fil.lower()+'_'+genstr+'_sim-mw.fits'
                        new_filename_ns='hlsp_vela_'+obs+'_'+instrument+'_'+dirname.lower()+'-'+cam+'-'+target_dir[14:-8]+'_'+fil.lower()+'_'+genstr+'_sim-ns.fits'
                        new_filename_smc='hlsp_vela_'+obs+'_'+instrument+'_'+dirname.lower()+'-'+cam+'-'+target_dir[14:-8]+'_'+fil.lower()+'_'+genstr+'_sim-smc.fits'

                        #bbhdu_ns=vela_export_image(hdulist,int(cam[-2:]),filfil[fil],nonscatter=True)

                        sys.stdout.flush()
                        
                        if os.path.lexists(os.path.join(outdir,new_filename)):
                            print(new_filename, ' exists, skipping.. ')
                            continue
                        else:
                            bbhdu=vela_export_image(hdulist,int(cam[-2:]),filfil[fil])
                            bbhdu_smc=vela_export_image(hdulist_smc,int(cam[-2:]),filfil[fil])
                            ns_hdu=vela_export_image(hdulist,int(cam[-2:]),filfil[fil],nonscatter=True)
                            print('saving.. ', new_filename)
                            
                        sys.stdout.flush()
                            
                        if not os.path.lexists(outdir):
                            os.makedirs(outdir)
                        
                        target_file=os.path.join(imagedir,target_dir[7:]+'_'+cam+'_'+instrumentfind+'-'+fil+'_SB00.fits')
                        target_file_ns=os.path.join(imagedir_ns,target_dir[7:]+'_nonscatter'+cam+'_'+instrumentfind+'-'+fil+'_SB00.fits')
                        target_file_smc=os.path.join(imagedir_smc,target_dir[7:]+'_'+cam+'_'+instrumentfind+'-'+fil+'_SB00.fits')
                        
                        if os.path.lexists(target_file) and os.path.lexists(target_file_ns) and os.path.lexists(target_file_smc):
                            ni=1
                        else:
                            ni=0
                            
                        #ni=np.where(np.asarray(names)==target_file)[0]

                        assert(ni==1)
                        
                        if ni==1:
                            #tfo.extractall(path=outdir,members=marr[ni])
                            shutil.copy2(target_file,outdir)
                            os.rename(os.path.join(outdir,os.path.basename(target_file)),os.path.join(outdir,new_filename))
                            #os.rmdir(os.path.join(outdir,os.path.dirname(target_file)))
                            shutil.copy2(target_file_ns,outdir)
                            os.rename(os.path.join(outdir,os.path.basename(target_file_ns)),os.path.join(outdir,new_filename_ns))
                            shutil.copy2(target_file_smc,outdir)
                            os.rename(os.path.join(outdir,os.path.basename(target_file_smc)),os.path.join(outdir,new_filename_smc))

                            for tfn,thdu,typename,hl in zip([new_filename,new_filename_ns,new_filename_smc],[bbhdu,ns_hdu,bbhdu_smc],['MW','NONSCATTER','SMCbar'],[hdulist,hdulist,hdulist_smc]):
                            
                                with pyfits.open(os.path.join(outdir,tfn),mode='update') as hdus:
                                    hdus.append(thdu)
                                    outhdu=hdus[0]
                                    outhdu.header['HLSPID']='vela'
                                    outhdu.header['HLSPLEAD']='Gregory F. Snyder'
                                    outhdu.header['HLSPNAME']='Vela-Sunrise Mock Images'
                                    outhdu.header['HLSPVER']=genstr
                                    outhdu.header['HLSPDOI']='10.17909/t9-ge0b-jm58'
                                    outhdu.header['LICENSE']='CC BY 4.0'
                                    outhdu.header['LICENURL']='https://creativecommons.org/licenses/by/4.0/'
                                    outhdu.header['PROPOSID']='HST-AR#13887'
                                    outhdu.header['REFERENC']='Simons et al. (2019)'
                                    outhdu.header['REFDOI']='10.3847/1538-4357/ab07c9'
                                    
                                    outhdu.header['MISSION']=obs.upper()
                                    outhdu.header['TELESCOP']=obs.upper()
                                    outhdu.header['INSTR']=instrument.upper()
                                    outhdu.header['INSTRUME']=instrument.upper()
                                    outhdu.header['EXTNAME']='IMAGE_PSF'
                                    outhdu.header['DUSTTYPE']=typename


    
    
    return


def extract_st_from_vela(dirname='VELA01',obslist=['hst','jwst','wfirst'],
                         camlist=cams, aux_only=False):

    imagedirs=np.sort(np.asarray(glob.glob(dirname+'/*_sunrise/images/images_*_sunrise')))


    
    for idir in imagedirs:
        print(idir)

        target_bb_fits=os.path.join(os.path.dirname(idir),'broadbandz.fits')
        if not os.path.lexists(target_bb_fits):
            print('unzipping.. ',target_bb_fits)
            subprocess.call(['gunzip', target_bb_fits])

        hdulist=pyfits.open(target_bb_fits)

        target_dir=os.path.basename(idir)
            
        for obs in obslist:
            for instrument in insdict[obs]:
                if instrument is 'wfc3':
                    instrumentfind='WFC3'
                elif instrument is 'acs':
                    instrumentfind='ACS'
                elif instrument is 'nircam':
                    instrumentfind='NC'
                elif instrument is 'miri':
                    instrumentfind='MIRI'
                elif instrument is 'wfidrm15':
                    instrumentfind='WFI_DRM15'
                for cam in camlist:
                    auxdir=os.path.join('HLSP','vela',dirname.lower(),cam)
                    if not os.path.lexists(auxdir):
                        os.makedirs(auxdir,exist_ok=True)
                    auxhdu=hdulist['CAMERA'+str(int(cam[-2:]))+'-AUX']
       
                    auxoutfile=os.path.join(auxdir,'hlsp_vela_none_none_'+dirname.lower()+'-'+cam+'-'+target_dir[14:-8]+'_aux_v3_sim.fits')
                    if aux_only is False:
                        if not os.path.lexists(auxoutfile):
                            auxhdu.writeto(auxoutfile,overwrite=False)
                    else:
                        auxhdu.writeto(auxoutfile,overwrite=True)

                    for fil in fildict[instrument]:
                        if fil=='aux' or aux_only is True:
                            continue
                            
                        outdir=os.path.join('HLSP','vela',dirname.lower(),cam,obs,instrument,fil.lower(),'')



                        new_filename='hlsp_vela_'+obs+'_'+instrument+'_'+dirname.lower()+'-'+cam+'-'+target_dir[14:-8]+'_'+fil.lower()+'_v3'+'_sim.fits'
                        #bbhdu_ns=vela_export_image(hdulist,int(cam[-2:]),filfil[fil],nonscatter=True)

                        sys.stdout.flush()
                        
                        if os.path.lexists(os.path.join(outdir,new_filename)):
                            print(new_filename, ' exists, skipping.. ')
                            continue
                        else:
                            bbhdu=vela_export_image(hdulist,int(cam[-2:]),filfil[fil])
                            print('saving.. ', new_filename)
                            
                        sys.stdout.flush()
                            
                        if not os.path.lexists(outdir):
                            os.makedirs(outdir,exist_ok=True)
                        
                        target_file=os.path.join(idir,target_dir[7:]+'_'+cam+'_'+instrumentfind+'-'+fil+'_SB00.fits')

                        if os.path.lexists(target_file):
                            ni=1
                        else:
                            ni=0
                            
                        #ni=np.where(np.asarray(names)==target_file)[0]

                        if ni==1:
                            #tfo.extractall(path=outdir,members=marr[ni])
                            shutil.copy2(target_file,outdir)
                            os.rename(os.path.join(outdir,os.path.basename(target_file)),os.path.join(outdir,new_filename))
                            #os.rmdir(os.path.join(outdir,os.path.dirname(target_file)))
                            with pyfits.open(os.path.join(outdir,new_filename),mode='update') as hdus:
                                hdus.append(bbhdu)
                                outhdu=hdus[0]
                                outhdu.header['HLSPID']='vela'
                                outhdu.header['HLSPLEAD']='Gregory F. Snyder'
                                outhdu.header['HLSPNAME']='Vela-Sunrise Mock Images'
                                outhdu.header['HLSPVER']='v3'
                                outhdu.header['LICENSE']='CC BY 4.0'
                                outhdu.header['LICENURL']='https://creativecommons.org/licenses/by/4.0/'
                                outhdu.header['PROPOSID']='HST-AR#13887'
                                outhdu.header['REFERENC']='Simons et al. (2019)'
                                outhdu.header['MISSION']=obs.upper()
                                outhdu.header['TELESCOP']=obs.upper()
                                outhdu.header['INSTR']=instrument.upper()
                                outhdu.header['INSTRUME']=instrument.upper()
                                outhdu.header['EXTNAME']='IMAGE_PSF'
                                    
                        elif ni==0 and obs=='wfirst':
                            bbhdu.writeto(os.path.join(outdir,new_filename),overwrite=True)
                            with pyfits.open(os.path.join(outdir,new_filename),mode='update') as hdus:
                                outhdu=hdus[0]
                                outhdu.header['HLSPID']='vela'
                                outhdu.header['HLSPLEAD']='Gregory F. Snyder'
                                outhdu.header['HLSPNAME']='Vela-Sunrise Mock Images'
                                outhdu.header['HLSPVER']='v3'
                                outhdu.header['LICENSE']='CC BY 4.0'
                                outhdu.header['LICENURL']='https://creativecommons.org/licenses/by/4.0/'
                                outhdu.header['PROPOSID']='HST-AR#13887'
                                outhdu.header['REFERENC']='Simons et al. (2019)'
                                outhdu.header['MISSION']=obs.upper()
                                outhdu.header['TELESCOP']=obs.upper()
                                outhdu.header['INSTR']=instrument.upper()
                                outhdu.header['INSTRUME']=instrument.upper()
                                
    return
    


def extract_jwst_from_vela(dirname='VELA01'):

    tarfiles=np.asarray(glob.glob(dirname+'/*_sunrise/images/images*.tar'))

    outdir='jwst_'+dirname
    if not os.path.lexists(outdir):
        os.mkdir(outdir)

    print(tarfiles)

    for tf in tarfiles:
        print(tf)
        with tarfile.TarFile(tf,mode='r') as tfo:
            names=tfo.getnames()
            members=tfo.getmembers()
            blist=[]
            for n,m in zip(names,members):
                if n.find('NC')>0 and n.find('SB00')>0 :
                    blist.append(True)
                elif n.find('MIRI')>0 and n.find('SB00')>0:
                    blist.append(True)
                else:
                    blist.append(False)

            barr=np.asarray(blist)
            print(barr.shape)
            narr=np.asarray(names)
            marr=np.asarray(members)

            print(narr[barr].shape)

            tfo.extractall(path=outdir,members=marr[barr])

    return outdir


def phot_from_vela(jd="stamps_VELA01"):

    fitsfiles=np.sort(np.asarray(glob.glob(jd+'/images_*_sunrise/*.fits')))
    print(fitsfiles)

    outf=open(jd+'_photometry.txt','w')
    outf.write('name,redshift,filter,apparentmag,distancemodulus\n')
    for fn in fitsfiles:
        header=pyfits.open(fn)[0].header
        mag=header['MAG']
        filname=header['FILTER']
        redshift=header['REDSHIFT']
        distmod=header['DISTMOD']
        print(fn,redshift,mag)
        outf.write('{:20s} {:8.4f} {:20s} {:8.4f} {:8.4f}\n'.format(jd,redshift,filname,mag,distmod))

    outf.close()

    return 

def untar_vela(dirname='VELA01'):
    tarfiles=np.sort(np.asarray(glob.glob(os.path.abspath(dirname+'/*_sunrise/images/images*.tar'))))

    cwd=os.path.abspath(os.curdir)
    
    for tf in tarfiles:
        dd=os.path.dirname(tf)
        os.chdir(dd)
        print('untarring.. ',tf,' into ', dd)
        sys.stdout.flush()
        subprocess.call(['tar', 'xf', os.path.basename(tf)])
        os.chdir(cwd)
        
    return



def get_nonscatter_hdu(sim,scalestr,camstr,filtername):

    dirname=sim.upper()
    astr='a'+scalestr
    bbf=os.path.join(dirname,dirname+'_'+astr+'_sunrise','images','broadbandz.fits')
    #print(bbf)
    #imagedirs=np.sort(np.asarray(glob.glob(dirname+'/*_sunrise/images/images_*_sunrise')))
    #print(int(camstr[-2:]))

    try:
        hdulist=pyfits.open(bbf,mode='readonly')
        #def vela_export_image(hdulist,camnum,filtername,label='',nonscatter=False):
        ns_hdu=vela_export_image(hdulist,int(camstr[-2:]),filtername,nonscatter=True)
        return ns_hdu
    except:
        print('Crap!')
        return None



def parse_vela_files(dirname='VELA01'):

    #update final bibcode reference?

    #create manifest of all FITS files /vela/vela??/etc and catalog of basic parameters
    #link with galprops info?? a la luvoir sims
    
    #file, sim, redshift, camera, telescope, instrument, filter, pristine apparent mag, stellar mass? , halo mass?
    #add the new mag to the pristine header?

    image_files=np.sort(np.asarray(glob.glob(os.path.join('HLSP','vela',dirname.lower(),'*/*/*/*/*.fits'))))

    aux_files=np.sort(np.asarray(glob.glob(os.path.join('HLSP','vela',dirname.lower(),'*/*.fits'))))
    
    catfile=os.path.abspath(os.path.join('HLSP','vela','catalogs','hlsp_vela_multi_multi_'+dirname.lower()+'_multi'+'_v3'+'_cat.txt'))
    #                        new_filename='hlsp_vela_'+obs+'_'+instrument+'_'+dirname.lower()+'-'+cam+'-'+target_dir[14:-8]+'_'+fil.lower()+'_v3'+'_sim.fits'
    auxcatfile=os.path.abspath(os.path.join('HLSP','vela','catalogs','hlsp_vela_multi_multi_'+dirname.lower()+'_multi'+'_v3'+'_auxcat.txt'))

    print(catfile)
    print(auxcatfile)
    
    datf=os.path.join(dirname,dirname+'_galprops.npy')
    dat=np.load(datf,encoding='bytes').all()

    
    aux_tfo=open(auxcatfile,'w')
    aux_tfo.write('sim z scale cam mstar mgas mmet sfr mvir_dm path\n')
    

    sfr_tau = 1.2e7 #12 Myr, from Ceverino et al. 2015.  SFR values are actually SFR*Tau, so must divide by this factor to get true SFR values

    for auxfile in aux_files:
        print(auxfile)
        auxfo=pyfits.open(auxfile,mode='update')
        auxpri=auxfo[0]
        if 'HLSPID' in auxpri.header:
            print('Aux file already updated, skipping..', auxfile)
            continue
        auxpri.header['HLSPID']='vela'
        auxpri.header['HLSPLEAD']='Gregory F. Snyder'
        auxpri.header['HLSPNAME']='Vela-Sunrise Mock Images'
        auxpri.header['HLSPVER']='v3'
        auxpri.header['LICENSE']='CC BY 4.0'
        auxpri.header['LICENURL']='https://creativecommons.org/licenses/by/4.0/'
        auxpri.header['PROPOSID']='HST-AR#13887'
        auxpri.header['REFERENC']='Simons et al. 2019'

        auxmain=auxfo[1]
        kpc_per_pix=auxmain.header['CD1_1']
        stellar_mass_per_pix=auxmain.data[4,:,:]*(kpc_per_pix**2)
        gas_mass_per_pix=auxmain.data[0,:,:]*(kpc_per_pix**2)
        metal_mass_per_pix=auxmain.data[1,:,:]*(kpc_per_pix**2)
        
        sfrdens_per_pix=auxmain.data[2,:,:]/sfr_tau
        auxmain.data[2,:,:]=sfrdens_per_pix
        
        sfr_per_pix=sfrdens_per_pix*(kpc_per_pix**2)
        
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
            writestr='{:10s} {:15.8f} {:10s} {:10s} {:15.6e} {:15.6e} {:15.6e} {:15.6e} {:15.6e} {:75s}\n'.format(sim,zfloat,scalestr,camstr,total_mstar,total_mgas,total_mmet,total_sfr,mvirdm[0],auxfile[5:])
            aux_tfo.write(writestr)
        except:
            print('weird problem with mvirdm?')
            
        aux_tfo.flush()
        sys.stdout.flush()
        
        
    aux_tfo.close()
    
    aux_table=ascii.read(auxcatfile)
    
    cat_tfo=open(catfile,'w')
    cat_tfo.write('sim z scale cam mission instrument filter flux_njy abmag mstar mgas mmet sfr mvir_dm path\n')

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
        zfloat=(1.0/scalefloat)-1.0
        
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
        catstr='{:10s} {:12.8f} {:10s} {:10s} {:10s} {:10s} {:10s} {:15.6e} {:12.6f} {:15.6e} {:15.6e} {:15.6e} {:15.6e} {:15.6e} {:75s}\n'.format(sim,zfloat,scalestr,camstr,mission,instr,filname,
                                                                                                                                                 total_flux_njy,mag_val,total_mstar,total_mgas,total_mmet,total_sfr,mvirdm[0],imfile[5:])
        
        cat_tfo.write(catstr)
        cat_tfo.flush()
        sys.stdout.flush()
        

    cat_tfo.close()
        
    return


def retar_vela_files(dirname='VELA01'):

    hlsp_sim_dir=os.path.join('HLSP','vela',dirname.lower())
    print(hlsp_sim_dir)
    
    cwd=os.path.abspath(os.curdir)
    os.chdir(hlsp_sim_dir)

    print('Tarring.. ',  hlsp_sim_dir)

    obslist=['hst','jwst','wfirst','aux']
    for obs in obslist:
        for instrument in insdict[obs]:
            for filname in fildict[instrument]:
                lfil=filname.lower()
                print(obs,instrument,lfil)
                if lfil=='aux':
                    imagefiles=np.sort(np.asarray(glob.glob('cam*/*.fits')))
                    tarfilename='hlsp_vela_none_none_'+dirname.lower()+'_'+lfil+'_v3_sim.tar'
                else:
                    imagefiles=np.sort(np.asarray(glob.glob('cam*/'+obs+'/'+instrument+'/'+lfil+'/*.fits')))
                    tarfilename='hlsp_vela_'+obs+'_'+instrument+'_'+dirname.lower()+'_'+lfil+'_v3_sim.tar'

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
                
                
    os.chdir(cwd)
    return    



def retar_vela_files_by_filter():

    hlsp_sim_dir=os.path.join('HLSP','vela')
    print(hlsp_sim_dir)
    
    cwd=os.path.abspath(os.curdir)
    os.chdir(hlsp_sim_dir)

    print('Tarring.. ',  hlsp_sim_dir)

    obslist=['hst','jwst','wfirst']
    for obs in obslist:
        for instrument in insdict[obs]:
            for filname in fildict[instrument]:
                lfil=filname.lower()
                print(obs,instrument,lfil)

                imagefiles=np.sort(np.asarray(glob.glob('vela??/cam*/'+obs+'/'+instrument+'/'+lfil+'/*.fits')))
                tarfilename='hlsp_vela_'+obs+'_'+instrument+'_vela_'+lfil+'_v3_sim.tar'

                print(tarfilename)
                print(imagefiles[0:25])
                print(imagefiles.shape)

                sys.stdout.flush()
                
                tfo=tarfile.open(tarfilename,mode='a')
                i=0
                for imf in imagefiles:
                    i=i+1
                    print(imf)
                    tfo.add(imf)
                    if i % 100==0:
                        sys.stdout.flush()
                        
                tfo.close()
                sys.stdout.flush()
                
                
    os.chdir(cwd)
    return    

if __name__=="__main__":

    vdir= np.sort(np.asarray(glob.glob('VELA??')))

    #untar_vela(dirname='VELA02')
    
    #extract_st_from_vela(dirname='VELA03',aux_only=True)
    
    #parse_vela_files(dirname='VELA03')
    
    #parse_vela_files(dirname='VELA02')
    #parse_vela_files(dirname='VELA03')

    
    
    #retar_vela_files(dirname='VELA01')

    #retar_vela_files_by_filter()


    #append_nonscatters(dirname='VELA01') ??
    
    
    for vd in vdir:
        retar_vela_files(dirname=vd)

#        if vd=='VELA01':
#            continue


#        parse_vela_files(dirname=vd)
        
        #untar_vela(dirname=vd)
        
        #extract_st_from_vela(dirname=vd)


