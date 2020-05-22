import os
import glob
import numpy as np



if __name__=="__main__":


    dcf='Ceverino_centers_gen6.dat'
    #VELA01/profileSpGe_Reca0.200.dat:# (xc,yc,zc)[Mpc h-1 comoving]=    5.1572    5.0681    4.3615
    dcfo=open(dcf)
    lines=dcfo.readlines()
    this_lines=[]

    for l in lines:
        line_sim=l[0:6]
        ss=l[23:28]
        nums=l[-28:-1]
        numsplit=np.float64(np.asarray(nums.split('   ')))
        print('{:4.2f}  VELA{:2s}  {:12.6f} {:12.6f} {:12.6f}'.format(np.float64(ss),line_sim[-2:],numsplit[0],numsplit[1],numsplit[2]))

    dcfo.close()
