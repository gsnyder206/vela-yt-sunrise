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
import visualize_vela_hlsp as vvh


#vvh runs as vvh.make_vela_stamps(bd=output_dir,sim=dirname.lower(),cam=camst,single_aname=aname)

#obstrs=['hst','hst','hst','hst','hst','hst','jwst','jwst','jwst','roman']
#instrs=['wfc3','acs','acs','acs','wfc3','wfc3','nircam','nircam','nircam','wfi']
#fstrs=['f336w','f435w','f606w','f814w','f125w','f160w','f277w','f200w','f444w','w149']

def run_from_aux(aux_tfo,limit=5):

    i=0
    tarinfo=1

    while tarinfo is not None:
        tarinfo=aux_tfo.next()

        auxn=tarinfo.name
        aname = auxn.split('_')[-4].split('-')[-1]
        camname = auxn.split('_')[-4].split('-')[-2]

        print(aname, camname, auxn)

        i=i+1
        if i>=5:
            break


    return




if __name__=="__main__":
    #assume run within folder containing all hlsp*.tar.gz files for a sim, folder assumed to be called "vela??"

    #need to specify genstring
    genname=sys.argv[1]

    dir=os.path.abspath('.')
    simname=os.path.basename(dir)

    aux_file = 'hlsp_vela_none_none_'+simname+'_aux_'+genname+'_sim.tar*'

    #get list of cams and snaps?

    aux_tfo = tarfile.open(aux_file,'r')

    #SLOW
    #aux_names = aux_tfo.getnames()

    run_from_aux(aux_tfo)

    aux_tfo.close()

    #hlsp_files = np.sort(np.asarray(glob.glob('hlsp*_'+genname+'_sim-*.tar.gz')))
