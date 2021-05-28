import os
import glob
import numpy as np
import astropy.io.ascii as ascii


if __name__=="__main__":

    input_folders=np.sort(np.asarray(glob.glob('/Users/gsnyder/Dropbox/AR-Disks/Vela/gen6_centers/*')))

    nmfo=open('Mandelker_centers_gen6_formatted.txt','w')

    for indir in input_folders:
        simname='VELA'+os.path.basename(indir)[-2:]
        print(simname)

        this_data=ascii.read(os.path.join(indir,'centers.txt'),data_start=1)
        this_scales=this_data['col1']
        this_x=this_data['col2']  #cMpc/h
        this_y=this_data['col3']
        this_z=this_data['col4']

        for scale,x,y,z in zip(this_scales,this_x,this_y,this_z):
            nmfo.write('{:4.2f}  {}  {}  {}  {}\n'.format(scale,simname,x,y,z))

    nmfo.close()
