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
    tc=pd['true_center']  #in physical kpc, a list of arrays
    print(tc)
    a=pd['scale']
    x=np.zeros_like(a)
    y=np.zeros_like(a)
    z=np.zeros_like(a)

    #convert galprops centers to comoving
    for i,scale,cen in enumerate(zip(a,tc)):
        print('{:5.2f} {:10.3f}'.format(scale, cen[0]*(0.70)/scale))
        x[i]=cen[0]*0.70/scale
        y[i]=cen[1]*0.70/scale
        z[i]=cen[2]*0.70/scale

    #read ceverino centers file and parse into arrays
    dcf='~/vela_data/ceverino_centers_gen6.txt'
    #VELA01/profileSpGe_Reca0.200.dat:# (xc,yc,zc)[Mpc h-1 comoving]=    5.1572    5.0681    4.3615
    dcfo=open(dcf)
    lines=dcfo.readlines()
    this_lines=[]
    for l in lines:
        line_sim=l[0:6]
        if line_sim==simname:
            this_lines.append(l)

    print(this_lines)

    dcfo.close()


    #measure distance with ceverino centers


    #print out distances versus scalefactor
