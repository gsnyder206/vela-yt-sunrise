import os
import sys
import numpy as np
import glob
import shutil
import tarfile
import string
import pandas
import astropy
import astropy.cosmology
import astropy.io.fits as pyfits
import subprocess

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
    outhdu.header['EXTNAME']='IMAGE_PRISTINE'
    outhdu.header['FILTER']=filtername
    outhdu.header['CAMERA']=(camnum,'Sunrise camera number')
    
    return outhdu

def extract_st_from_vela(dirname='VELA01',obslist=['hst','jwst','wfirst'],
                         camlist=cams):

    tarfiles=np.sort(np.asarray(glob.glob(dirname+'/*_sunrise/images/images*.tar')))


    
    for tf in tarfiles:
        print(tf)

        target_bb_fits=os.path.join(os.path.dirname(tf),'broadbandz.fits')
        if not os.path.lexists(target_bb_fits):
            print('unzipping.. ',target_bb_fits)
            subprocess.call(['gunzip', target_bb_fits])

        hdulist=pyfits.open(target_bb_fits)

        target_dir=os.path.basename(tf)[0:-4]
            
        with tarfile.TarFile(tf,mode='r') as tfo:
            names=tfo.getnames()
            members=tfo.getmembers()
            marr=np.asarray(members)

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
                        if not os.path.lexists(auxoutfile):
                            auxhdu.writeto(auxoutfile,overwrite=False)
                        
                        for fil in fildict[instrument]:
                            if fil=='aux':
                                continue
                            
                            outdir=os.path.join('HLSP','vela',dirname.lower(),cam,obs,instrument,fil.lower(),'')
                            if not os.path.lexists(outdir):
                                os.makedirs(outdir,exist_ok=True)


                                
                            
                            target_file=os.path.join(target_dir,target_dir[7:]+'_'+cam+'_'+instrumentfind+'-'+fil+'_SB00.fits')
                            
                            ni=np.where(np.asarray(names)==target_file)[0]

                            new_filename='hlsp_vela_'+obs+'_'+instrument+'_'+dirname.lower()+'-'+cam+'-'+target_dir[14:-8]+'_'+fil.lower()+'_v3'+'_sim.fits'
                            bbhdu=vela_export_image(hdulist,int(cam[-2:]),filfil[fil])

                            if os.path.lexists(os.path.join(outdir,new_filename)):
                                print(new_filename, ' exists, skipping.. ')
                                continue
                            
                            if ni.shape[0]==1:
                                tfo.extractall(path=outdir,members=marr[ni])
                                os.rename(os.path.join(outdir,target_file),os.path.join(outdir,new_filename))
                                os.rmdir(os.path.join(outdir,os.path.dirname(target_file)))
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
                                    outhdu.header['REFERENC']='in prep'
                                    outhdu.header['MISSION']=obs.upper()
                                    outhdu.header['TELESCOP']=obs.upper()
                                    outhdu.header['INSTR']=instrument.upper()
                                    outhdu.header['INSTRUME']=instrument.upper()
                                    outhdu.header['EXTNAME']='IMAGE_PSF'
                                    
                            elif ni.shape[0]==0 and obs=='wfirst':
                                bbhdu.writeto(os.path.join(outdir,new_filename))
                                with pyfits.open(os.path.join(outdir,new_filename),mode='update') as hdus:
                                    outhdu=hdus[0]
                                    outhdu.header['HLSPID']='vela'
                                    outhdu.header['HLSPLEAD']='Gregory F. Snyder'
                                    outhdu.header['HLSPNAME']='Vela-Sunrise Mock Images'
                                    outhdu.header['HLSPVER']='v3'
                                    outhdu.header['LICENSE']='CC BY 4.0'
                                    outhdu.header['LICENURL']='https://creativecommons.org/licenses/by/4.0/'
                                    outhdu.header['PROPOSID']='HST-AR#13887'
                                    outhdu.header['REFERENC']='in prep'
                                    outhdu.header['MISSION']=obs.upper()
                                    outhdu.header['TELESCOP']=obs.upper()
                                    outhdu.header['INSTR']=instrument.upper()
                                    outhdu.header['INSTRUME']=instrument.upper()
                                
    return outdir
    


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



def parse_vela_files(dirname='VELA01'):

    return

if __name__=="__main__":

    vdir= np.sort(np.asarray(glob.glob('VELA??')))
    for vd in vdir:
        if vd=='VELA01' or vd=='VELA02':
            continue
        extract_st_from_vela(dirname=vd)


        #parse_vela_files(dirname=vd)
