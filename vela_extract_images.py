import os
import sys
import numpy as np
import glob
import shutil
import tarfile
import string
import pandas
import astropy
import astropy.io.fits as pyfits

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

    outf=open(jd+'_photometry.txt','w')
    outf.write('name,redshift,filter,apparentmag,distancemodulus')
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


if __name__=="__main__":

    #vdir= np.asarray(glob.glob('VELA??'))
    #for vd in vdir:
    #    extract_jwst_from_vela(vd)

    phot_from_vela()
