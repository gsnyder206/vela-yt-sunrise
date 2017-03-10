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
                if n.find('NC')>0 and n.find('SB00')>0 and (n.find('cam17')>0 or n.find('cam18')>0 or n.find('cam16')>0 ):
                    blist.append(True)
                elif n.find('MIRI')>0 and n.find('SB00')>0:
                    blist.append(False)
                else:
                    blist.append(False)

            barr=np.asarray(blist)
            print(barr.shape)
            narr=np.asarray(names)
            marr=np.asarray(members)

            print(narr[barr].shape)

            tfo.extractall(path=outdir,members=marr[barr])

    return
