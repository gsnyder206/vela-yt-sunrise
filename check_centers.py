import os
import glob
import numpy as np
import astropy
import astropy.io.ascii as ascii

if __name__=="__main__":
    #identify galprops file, read it
    #identify simulation
    pf=glob.glob('*_galprops.npy')[0]
    simname=pf[0:6]
    genname=os.path.basename(os.path.dirname(os.path.abspath('.')))
    print(simname)
    print(genname)

    pd=np.load(pf).all()
    tc=pd['true_center']  #in physical kpc, a list of arrays
    scale_exact=pd['scale_exact']  #vector of exact scale factor floats
    #print(tc)
    try:
        print(pd['centerfile'])
    except:
        print('No centerfile key in this galprops file.')

    aname=pd['scale_string']
    x=np.ndarray(shape=aname.shape, dtype=np.float64)
    y=np.ndarray(shape=aname.shape, dtype=np.float64)
    z=np.ndarray(shape=aname.shape, dtype=np.float64)

    #convert galprops centers to comoving kpc/h
    for i,(aname,cen,scale) in enumerate(zip(aname,tc,scale_exact)):
        print('{:5.2f} {:10.3f}'.format(np.float64(scale), cen[0]*(0.70)/np.float64(scale)))
        x[i]=cen[0]*0.70/np.float64(scale)
        y[i]=cen[1]*0.70/np.float64(scale)
        z[i]=cen[2]*0.70/np.float64(scale)



    dcf='/u/gfsnyder/PythonCode/vela-yt-sunrise/Ceverino_centers_'+genname.lower()+'_formatted.txt'
    cdata=ascii.read(dcf)
    si=np.where(cdata['col2']==simname)
    dc_x=np.asarray(cdata['col3'][si]*1000)  #in comoving kpc/h
    dc_y=np.asarray(cdata['col4'][si]*1000)
    dc_z=np.asarray(cdata['col5'][si]*1000)
    dc_a=np.asarray(cdata['col1'][si])

    print(dc_a)

    ''' old version
    #read ceverino centers file and parse into arrays
    dcf='ceverino_centers_gen6.txt'
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
            dc_a.append(np.float64(ss))
            nums=l[-28:-1]
            dc_x.append(1000*np.float64(nums.split('   ')[0]))
            dc_y.append(1000*np.float64(nums.split('   ')[1]))
            dc_z.append(1000*np.float64(nums.split('   ')[2]))
    dcfo.close()

    dc_x=np.asarray(dc_x)
    dc_y=np.asarray(dc_y)
    dc_z=np.asarray(dc_z)
    dc_a=np.asarray(dc_a)
    '''

    #for a,dx in zip(dc_a,dc_x): print(a,dx)

    #measure distance with ceverino centers
    d_ckpch=np.ndarray(shape=aname.shape, dtype=np.float64)
    print('#scale     3D distance (kpc)')
    for i,(aname,scale,tc_x,tc_y,tc_z) in enumerate(zip(aname,scale_exact,x,y,z)):
        #match to dc_a
        dc_i=np.where(dc_a==np.float64(aname[1:]))[0]
        #print(scale, dc_i, tc_x, dc_x[dc_i])
        if len(dc_i)==1:
            #print(scale,dc_x[dc_i],tc_x)
            d_ckpch[i]=((tc_x-dc_x[dc_i[0]])**2 + (tc_y-dc_y[dc_i[0]])**2 + (tc_z-dc_z[dc_i[0]])**2)**(0.5)
            #print out distances versus scalefactor
            print('{:8.3f}     {:12.4f}'.format(np.float64(scale),d_ckpch[i]*np.float64(scale)/0.70))
