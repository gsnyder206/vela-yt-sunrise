import os
import sys
import numpy as np
import glob
import shutil
import tarfile
import string


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


def phot_from_vela(jd="jwst_VELA01"):

    fitsfiles=np.sort(np.asarray(glob.glob(jd+'/images_*_sunrise/*.fits')))
    print(fitsfiles)

    return 


if __name__=="__main__":

    vdir= np.asarray(glob.glob('VELA??'))
    for vd in vdir:
        extract_jwst_from_vela(vd)
