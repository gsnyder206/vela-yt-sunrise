import os
import glob
import numpy as np


if __name__=="__main__":
    #identify galprops file, read it
    #identify simulation
    pf=glob.glob('*_galprops.npy')[0]
    simname=pf[0:6]
    print(simname)

    pd=np.load(pf).all()
    tc=pd['true_center']  #in physical kpc
    print(tc)

    #read ceverino centers file and parse into arrays

    #convert galprops centers to comoving

    #measure distance with ceverino centers


    #print out distances versus scalefactor
