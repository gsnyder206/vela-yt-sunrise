import os
import sys
import numpy as np
import glob
import shutil
import tarfile
import string
import matplotlib
matplotlib.use('Agg')  #apparently needed for old matplotlib/python on pleiades (didn't want to risk updating)
import matplotlib.pyplot as pyplot
#import pandas
import astropy
import astropy.cosmology
import astropy.io.fits as fits
import astropy.io.ascii as ascii

import vela_extract_images as vei
import visualize_vela_hlsp as vvh


#vvh runs as vvh.make_vela_stamps(bd=output_dir,sim=dirname.lower(),cam=camst,single_aname=aname)

obstrs=['hst','hst','hst','hst','hst','hst','jwst','jwst','jwst','roman']
instrs=['wfc3','acs','acs','acs','wfc3','wfc3','nircam','nircam','nircam','wfi']
fstrs=['f336w','f435w','f606w','f814w','f125w','f160w','f277w','f200w','f444w','w149']

def run_from_aux(genname,simname,aux_tfo,tarfodd,limit=20):

    i=0
    tarinfo=1

    while tarinfo is not None:
        tarinfo=aux_tfo.next()

        auxn=tarinfo.name
        aname = auxn.split('_')[-4].split('-')[-1]
        camname = auxn.split('_')[-4].split('-')[-2]

        print(aname, camname, auxn)

        aux_i = aux_tfo.extractfile(tarinfo)
        aux_fo_i = fits.open(aux_i,mode='readonly')

        make_vela_stamps_hlsp(aux_fo_i,tarfodd,sim=simname,cam=camname,single_aname=aname,genname=genname)


        i=i+1
        if limit != None:
            if i>=limit:
                break


    return



#bd='/nobackup/gfsnyder/VELA_sunrise/Outputs/HLSP/',sim='vela06',cam='cam00',single_aname=None

def make_vela_stamps_hlsp(aux_fo_i,tarfodd,sim='vela06',cam='cam00',single_aname=None,genname='v6'):


    aname=single_aname
    genstr=genname
    auxn=str(int(cam[-2:]))

    stardens=aux_fo_i['CAMERA'+auxn+'-AUX'].data[4,:,:]
    gasdens=aux_fo_i['CAMERA'+auxn+'-AUX'].data[0,:,:]
    sfrdens=aux_fo_i['CAMERA'+auxn+'-AUX'].data[1,:,:]

    bd='.'

    dd={}
    dd['mw']={}
    dd['ns']={}
    dd['smc']={}



    preview_dir=os.path.join(cam,'previews')
    if not os.path.lexists(preview_dir):
        os.makedirs(preview_dir)

    preview_fn=os.path.join(preview_dir,'hlsp_vela_multi_multi_'+sim+'-'+cam+'-'+aname+'_multi_'+genstr+'_preview.png')

    fig=pyplot.figure(figsize=(12.0,8.0),dpi=600)
    fig.subplots_adjust(left=0.0,bottom=0.0,top=1.0,right=1.0,hspace=0.0,wspace=0.0)

    ax1=fig.add_subplot(2,3,1)
    vvh.make_quant_stamp(gasdens,ax1,text='gas mass')
    ax4=fig.add_subplot(2,3,4)
    vvh.make_quant_stamp(stardens+5000.0,ax4,text='star mass')

    try:
        for obs,ins,f in zip(obstrs,instrs,fstrs):

            #tar_file_mw='hlsp_vela_'+obs+'_'+ins+'_'+sim+'_'+f+'_'+genstr+'_sim-mw.tar.gz'
            #tar_file_ns='hlsp_vela_'+obs+'_'+ins+'_'+sim+'_'+f+'_'+genstr+'_sim-ns.tar.gz'
            #tar_file_smc='hlsp_vela_'+obs+'_'+ins+'_'+sim+'_'+f+'_'+genstr+'_sim-smc.tar.gz'

            #tarfo_mw=tarfile.open(tar_file_mw,'r')
            #tarfo_ns=tarfile.open(tar_file_ns,'r')
            #tarfo_smc=tarfile.open(tar_file_smc,'r')


            tfn=os.path.join(cam,obs,ins,f,'hlsp_vela_'+obs+'_'+ins+'_'+sim+'-'+cam+'-'+aname+'_'+f+'_'+genstr+'_sim-mw.fits')
            tfn_ns=os.path.join(cam,obs,ins,f,'hlsp_vela_'+obs+'_'+ins+'_'+sim+'-'+cam+'-'+aname+'_'+f+'_'+genstr+'_sim-ns.fits')
            tfn_smc=os.path.join(cam,obs,ins,f,'hlsp_vela_'+obs+'_'+ins+'_'+sim+'-'+cam+'-'+aname+'_'+f+'_'+genstr+'_sim-smc.fits')

            print(tfn,tfn_ns,tfn_smc)

            etfn=tarfodd['mw'][f].extractfile(tfn)
            etfn_ns=tarfodd['ns'][f].extractfile(tfn_ns)
            etfn_smc=tarfodd['smc'][f].extractfile(tfn_smc)

            tfo=fits.open(etfn,'readonly')
            tfo_ns=fits.open(etfn_ns,'readonly')
            tfo_smc=fits.open(etfn_smc,'readonly')


            dd['mw'][f]=tfo
            dd['ns'][f]=tfo_ns
            dd['smc'][f]=tfo_smc

    except:
        print("Could not open all expected files, skipping image previews, ", tfn)
        fig.savefig(preview_fn,dpi=600)
        pyplot.close(fig)
        return

    dd['redshift']=dd['mw']['f200w'][0].header['REDSHIFT']

    ax2=fig.add_subplot(2,3,2)
    vvh.make_combo_stamp(bd,sim,cam,aname,ax2,dd,key='IMAGE_PRISTINE',sigma_tuple=[0.0,0.0,0.0],dusttype='mw',text='starlight and dust')

    ax3=fig.add_subplot(2,3,3)
    vvh.make_combo_stamp(bd,sim,cam,aname,ax3,dd,key='IMAGE_PSF',sigma_tuple=[0.05,0.05,0.05],dusttype='mw',text='mock data')



    ax5=fig.add_subplot(2,3,5)
    vvh.make_combo_stamp(bd,sim,cam,aname,ax5,dd,key='IMAGE_PRISTINE',sigma_tuple=[0.0,0.0,0.0],dusttype='smc',text='SMCbar dust')

    ax6=fig.add_subplot(2,3,6)
    vvh.make_combo_stamp(bd,sim,cam,aname,ax6,dd,key='IMAGE_PRISTINE_NONSCATTER',sigma_tuple=[0.0,0.0,0.0],dusttype='ns',text='no dust RT')


    fig.savefig(preview_fn,dpi=600)
    pyplot.close(fig)



    return



if __name__=="__main__":
    #assume run within folder containing all hlsp*.tar.gz files for a sim, folder assumed to be called "vela??"

    #need to specify genstring
    genname=sys.argv[1]

    dir=os.path.abspath('.')
    simname=os.path.basename(dir)

    aux_file = 'hlsp_vela_none_none_'+simname+'_aux_'+genname+'_sim.tar.gz'

    #get list of cams and snaps?

    aux_tfo = tarfile.open(aux_file,'r')

    #SLOW
    #aux_names = aux_tfo.getnames()

    tarfodd={}
    tarfodd['mw']={}
    tarfodd['ns']={}
    tarfodd['smc']={}

    for obs,ins,f in zip(obstrs,instrs,fstrs):

        tar_file_mw='hlsp_vela_'+obs+'_'+ins+'_'+simname+'_'+f+'_'+genname+'_sim-mw.tar.gz'
        tar_file_ns='hlsp_vela_'+obs+'_'+ins+'_'+simname+'_'+f+'_'+genname+'_sim-ns.tar.gz'
        tar_file_smc='hlsp_vela_'+obs+'_'+ins+'_'+simname+'_'+f+'_'+genname+'_sim-smc.tar.gz'

        tarfo_mw=tarfile.open(tar_file_mw,'r')
        tarfo_ns=tarfile.open(tar_file_ns,'r')
        tarfo_smc=tarfile.open(tar_file_smc,'r')

        tarfodd['mw'][f]=tarfo_mw
        tarfodd['ns'][f]=tarfo_ns
        tarfodd['smc'][f]=tarfo_smc


    run_from_aux(genname,simname,aux_tfo,tarfodd)

    aux_tfo.close()

    for f in fstrs:
        tarfodd['mw'][f].close()
        tarfodd['ns'][f].close()
        tarfodd['smc'][f].close()


    #hlsp_files = np.sort(np.asarray(glob.glob('hlsp*_'+genname+'_sim-*.tar.gz')))
