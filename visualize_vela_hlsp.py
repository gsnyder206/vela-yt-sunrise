

import astropy.io.fits as fits
import os
import matplotlib
matplotlib.use('Agg')  #apparently needed for old matplotlib/python on pleiades (didn't want to risk updating)
import matplotlib.pyplot as pyplot
import make_color_image
import congrid
import glob
import numpy as np


def make_combo_stamp(bd,sim,cam,aname,ax,dd,key='IMAGE_PRISTINE',sigma_tuple=[0.0,0.0,0.0],dusttype='mw',text=''):

    #check for different shapes in case key='IMAGE_PSF'

    r1=dd[dusttype]['f277w'][key].data
    r2=dd[dusttype]['f200w'][key].data
    r3=dd[dusttype]['f444w'][key].data

    g1=dd[dusttype]['w146'][key].data

    b1=dd[dusttype]['f336w'][key].data
    b2=dd[dusttype]['f435w'][key].data
    b3=dd[dusttype]['f606w'][key].data

    iml=[r1,r2,r3,g1,b1,b2,b3]
    for ima in iml:
        print(np.sum(ima))
        
    if key=='IMAGE_PSF':
        mp=np.min([r1.shape[0],r2.shape[0],r3.shape[0],g1.shape[0],b1.shape[0],b2.shape[0],b3.shape[0]])
        for i,ima in enumerate(iml):
            ip=ima.shape[0]
            os=np.sum(ima)
            iml[i]=congrid.congrid(ima,(mp,mp))
            iml[i]=os*iml[i]/np.sum(iml[i])/5.0/(ip/mp)**2
            print(os,np.sum(iml[i]))

    
    rr=(iml[0]+iml[1]+iml[2])*0.75
    gg=(iml[3])*3.0
    bb=(iml[4]+iml[5]+iml[6])*2.0

    
    scalefactor=(1.0/(1.0+dd['redshift']))
    alph=10.0*(0.33/scalefactor)**5 ; Q=7.0

    rgbthing=make_color_image.make_interactive_nasa(bb,gg,rr,alph,Q,sigma_tuple=sigma_tuple)
    ax.set_axis_off()
    ax.imshow(rgbthing,interpolation='nearest',aspect='auto',origin='lower')

    ax.annotate(text,(0.5,0.15),xycoords='axes fraction',ha='center',va='center',color='white',size=14)
    
    return ax

def make_quant_stamp(quant,ax,text=''):

    ax.imshow(np.log10(np.flipud(np.transpose(quant))),interpolation='nearest',aspect='auto',cmap='magma')
    ax.annotate(text,(0.5,0.15),xycoords='axes fraction',ha='center',va='center',color='white',size=14)

    return ax

def make_vela_stamps(bd='/nobackup/gfsnyder/VELA_sunrise/Outputs/HLSP/',sim='vela06',cam='cam00',single_aname=None):



    auxes = np.sort(np.asarray(glob.glob(os.path.join(bd,'vela',sim,cam,'hlsp*aux*.fits'))))
    obstrs=['hst','hst','hst','hst','hst','hst','jwst','jwst','jwst','wfirst']
    instrs=['wfc3','acs','acs','acs','wfc3','wfc3','nircam','nircam','nircam','wfi']
    fstrs=['f336w','f435w','f606w','f814w','f125w','f160w','f277w','f200w','f444w','w146']

    auxn=str(int(cam[-2:]))
    
    for afn in auxes:
        aname=os.path.basename(afn).split('_')[4].split('-')[2]

        if single_aname is not None:
            if aname==single_aname:
                print('processing single aname: ', single_aname)
                pass
            else:
                continue
        
        #genstr based on aux file name -- should be one for each, so this is OK
        genstr=os.path.basename(afn).split('_')[6]
        
        stardens=fits.open(afn)['CAMERA'+auxn+'-AUX'].data[4,:,:]
        gasdens=fits.open(afn)['CAMERA'+auxn+'-AUX'].data[0,:,:]
        sfrdens=fits.open(afn)['CAMERA'+auxn+'-AUX'].data[1,:,:]

        dd={}
        dd['mw']={}
        dd['ns']={}
        dd['smc']={}

        
        
        preview_dir=os.path.join(bd,'vela',sim,cam,'previews')
        if not os.path.lexists(preview_dir):
            os.makedirs(preview_dir)
        preview_fn=os.path.join(preview_dir,'hlsp_vela_multi_multi_'+sim+'-'+cam+'-'+aname+'_multi_'+genstr+'_preview.png')

        fig=pyplot.figure(figsize=(12.0,8.0),dpi=600)
        fig.subplots_adjust(left=0.0,bottom=0.0,top=1.0,right=1.0,hspace=0.0,wspace=0.0)

        ax1=fig.add_subplot(2,3,1)
        make_quant_stamp(gasdens,ax1,text='gas mass')
        ax4=fig.add_subplot(2,3,4)
        make_quant_stamp(stardens+5000.0,ax4,text='star mass')
        
        try:
            for obs,ins,f in zip(obstrs,instrs,fstrs):
                tfn=os.path.join(bd,'vela',sim,cam,obs,ins,f,'hlsp_vela_'+obs+'_'+ins+'_'+sim+'-'+cam+'-'+aname+'_'+f+'_'+genstr+'_sim-mw.fits')
                tfn_ns=os.path.join(bd,'vela',sim,cam,obs,ins,f,'hlsp_vela_'+obs+'_'+ins+'_'+sim+'-'+cam+'-'+aname+'_'+f+'_'+genstr+'_sim-ns.fits')
                tfn_smc=os.path.join(bd,'vela',sim,cam,obs,ins,f,'hlsp_vela_'+obs+'_'+ins+'_'+sim+'-'+cam+'-'+aname+'_'+f+'_'+genstr+'_sim-smc.fits')

                print(tfn,tfn_ns,tfn_smc)
            
                tfo=fits.open(tfn)
                tfo_ns=fits.open(tfn_ns)
                tfo_smc=fits.open(tfn_smc)
                dd['mw'][f]=tfo
                dd['ns'][f]=tfo_ns
                dd['smc'][f]=tfo_smc
        except:
            print("Could not open all expected files, skipping image previews, ", afn)
            fig.savefig(preview_fn,dpi=600)
            pyplot.close(fig)
            continue

        dd['redshift']=dd['mw']['f200w'][0].header['REDSHIFT']
        
        ax2=fig.add_subplot(2,3,2)
        make_combo_stamp(bd,sim,cam,aname,ax2,dd,key='IMAGE_PRISTINE',sigma_tuple=[0.0,0.0,0.0],dusttype='mw',text='starlight and dust')
        
        ax3=fig.add_subplot(2,3,3)
        make_combo_stamp(bd,sim,cam,aname,ax3,dd,key='IMAGE_PSF',sigma_tuple=[0.05,0.05,0.05],dusttype='mw',text='mock data')


        
        ax5=fig.add_subplot(2,3,5)
        make_combo_stamp(bd,sim,cam,aname,ax5,dd,key='IMAGE_PRISTINE',sigma_tuple=[0.0,0.0,0.0],dusttype='smc',text='SMCbar dust')
        
        ax6=fig.add_subplot(2,3,6)
        make_combo_stamp(bd,sim,cam,aname,ax6,dd,key='IMAGE_PRISTINE_NONSCATTER',sigma_tuple=[0.0,0.0,0.0],dusttype='ns',text='no dust RT')

        
        #make_hst_stamp(bd,sim,cam,aname)

        #make_jwst_stamp(bd,sim,cam,aname)

        #make_wfirst_stamp(bd,sim,cam,aname)


        fig.savefig(preview_fn,dpi=600)
        pyplot.close(fig)

    

    return


