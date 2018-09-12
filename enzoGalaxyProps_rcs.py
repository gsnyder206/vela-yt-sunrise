import sys
import os
import glob
import yt
import numpy as np
from numpy import *
import astropy
from astropy.cosmology import Planck13 as cosmo
import findGalaxyProps as fGP





if __name__=="__main__":
    snaps = np.sort(np.asarray(glob.glob("RD????/RD????")))  #ENZO format a list of snapshots in separate directories
    form='ENZO'

    assert snaps.shape[0] > 0

    print("Calculating Galaxy Props for "+form+": ", snaps)

    abssnap = os.path.abspath(snaps[0])
    dirname = os.path.dirname(os.path.dirname(abssnap))
    simname = os.path.basename(dirname) #assumes directory name for simulation name

    print( "Simulation name:  ", simname)

    particle_headers = []
    particle_data = []
    stars_data = []
    new_snapfiles = []
    for sn in snaps:
        aname=os.path.basename(sn)
        adir=os.path.abspath(os.path.dirname(sn))
        snap_dir = os.path.join(adir,simname+'_'+aname+'_sunrise')
        yt_fig_dir = snap_dir+'/yt_projections'
        print( "Sunrise directory: ", snap_dir)
        if not os.path.lexists(snap_dir):
            print ("Creating Sunrise directory:", snap_dir)
            os.mkdir(snap_dir)        
        if not os.path.lexists(yt_fig_dir):
            print ("Creating YT figure directory:", yt_fig_dir)
            os.mkdir(yt_fig_dir)        


        new_snapfiles.append(os.path.abspath(sn))


    new_snapfiles = np.asarray(new_snapfiles)
    galaxy_props = {}
    fields = ['scale', 'stars_total_mass', 'stars_com', 'stars_maxdens', 'stars_maxndens', 'stars_hist_center',
              'stars_rhalf', 'stars_mass_profile', 'stars_L','gas_total_mass', 'gas_maxdens', 'gas_L', 'rvir', 
              'Mvir_dm', 'stars_center','snap_files']
    for field in fields: 
        if field in ['scale', 'stars_total_mass', 'stars_rhalf', 'gas_total_mass' ]:
            galaxy_props[field] = np.array([])                
        else:
            galaxy_props[field] = []



    def _stars(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 2

    #this gets dark matter particles in zoom region only
    def _darkmatter(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 4

    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    yt.add_particle_filter("darkmatter",function=_darkmatter, filtered_type='all',requires=["particle_type"])

    ts = yt.DatasetSeries(new_snapfiles)


    for ds,snap_dir in zip(reversed(ts),np.flipud(new_snapfiles)):
        print( "Getting galaxy props: ",  snap_dir)

        ds.add_particle_filter('stars')
        ds.add_particle_filter('darkmatter')

        dd = ds.all_data()
        ds.domain_right_edge = ds.arr(ds.domain_right_edge,'code_length')
        ds.domain_left_edge  = ds.arr(ds.domain_left_edge,'code_length')
        print(ds.index.get_smallest_dx())



























