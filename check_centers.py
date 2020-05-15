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

    #convert galprops centers to comoving kpc
    for i,(scale,cen) in enumerate(zip(a,tc)):
        print('{:5.2f} {:10.3f}'.format(scale, cen[0]*(0.70)/scale))
        x[i]=cen[0]*0.70/scale
        y[i]=cen[1]*0.70/scale
        z[i]=cen[2]*0.70/scale

    #read ceverino centers file and parse into arrays
    dcf='/u/gfsnyder/vela_data/ceverino_centers_gen6.txt'
    #VELA01/profileSpGe_Reca0.200.dat:# (xc,yc,zc)[Mpc h-1 comoving]=    5.1572    5.0681    4.3615
    dcfo=open(dcf)
    lines=dcfo.readlines()
    this_lines=[]
    dc_x=[]
    dc_y=[]
    dc_z=[]
    dc_a=[]
    for l in lines:
        line_sim=l[0:6]
        if line_sim==simname:
            this_lines.append(l)
            ss=l[23:28]
            dc_a.append(ss)
            nums=l[-28:-1]
            dc_x.append(1000*np.float64(nums.split('   ')[0]))
            dc_y.append(1000*np.float64(nums.split('   ')[1]))
            dc_z.append(1000*np.float64(nums.split('   ')[2]))
    dcfo.close()

    dc_x=np.asarray(dc_x)
    dc_y=np.asarray(dc_y)
    dc_z=np.asarray(dc_z)
    dc_a=np.asarray(dc_a)

    #for a,dx in zip(dc_a,dc_x): print(a,dx)

    #measure distance with ceverino centers

    for scale,truecen in zip(z,tc):
        #match to dc_a
        dc_i=dc_a==scale
        print(np.sum(dc_i))



    #print out distances versus scalefactor
