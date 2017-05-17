import sys
import os
import glob
import yt
import numpy as np
from numpy import *
import astropy
from astropy.cosmology import Planck13 as cosmo
reload(yt)
import findGalaxyProps as fGP

def _stars(pfilter, data):
    return data[(pfilter.filtered_type, "particle_type")] == 2


if __name__=="__main__":
    
    if len(sys.argv)==3:
        snaps = np.asarray([sys.argv[1]])
        form= sys.argv[2]
    else:
        try:
            snaps = np.sort(np.asarray(glob.glob("*.d")))  #VELA format
            form='VELA'
            assert snaps.shape[0] > 0
        except AssertionError as e:
            snaps = np.sort(np.asarray(glob.glob("RD????/RD????")))  #ENZO format a list of snapshots in separate directories
            form='ENZO'

    assert snaps.shape[0] > 0

    print "Calculating Galaxy Props for "+form+": ", snaps

    abssnap = os.path.abspath(snaps[0])
    assert os.path.lexists(abssnap)

    if form=='VELA':
        dirname = os.path.dirname(abssnap)
    elif form=='ENZO':
        dirname = os.path.dirname(os.path.dirname(abssnap))
        
    simname = os.path.basename(dirname) #assumes directory name for simulation name
    
    print "Simulation name:  ", simname

    particle_headers = []
    particle_data = []
    stars_data = []
    new_snapfiles = []
    for sn in snaps:
        if form=='VELA':
            aname = sn.split('_')[-1].rstrip('.d')
            particle_headers.append('PMcrd'+aname+'.DAT')
            particle_data.append('PMcrs0'+aname+'.DAT')
            stars_data.append('stars_'+aname+'.dat')
            snap_dir = os.path.join(simname+'_'+aname+'_sunrise')
        elif form=='ENZO':
            aname=os.path.basename(sn)
            adir=os.path.abspath(os.path.dirname(sn))
            
            snap_dir = os.path.join(adir,simname+'_'+aname+'_sunrise')
            
        yt_fig_dir = snap_dir+'/yt_projections'
            
        print "Sunrise directory: ", snap_dir
        if not os.path.lexists(snap_dir):
            os.mkdir(snap_dir)        
        if not os.path.lexists(yt_fig_dir):
            os.mkdir(yt_fig_dir)        


        if form=='VELA':
            #art snaps need to be isolated from others for some reason
            newf = os.path.join(snap_dir,sn)
            new_snapfiles.append(newf)
            if not os.path.lexists(newf):
                os.symlink(os.path.abspath(sn),newf)
                os.symlink(os.path.abspath(particle_headers[-1]),os.path.join(snap_dir,particle_headers[-1]))
                os.symlink(os.path.abspath(particle_data[-1]),os.path.join(snap_dir,particle_data[-1]))
                os.symlink(os.path.abspath(stars_data[-1]),os.path.join(snap_dir,stars_data[-1]))
        elif form=='ENZO':
            #can probably just work from snap directories?
            new_snapfiles.append(os.path.abspath(sn))
            
    
    new_snapfiles = np.asarray(new_snapfiles)
    galaxy_props = {}
    fields = ['scale', 'stars_total_mass', 'stars_com', 'stars_maxdens', 'stars_maxndens', 'stars_hist_center',
	      'stars_rhalf', 'stars_mass_profile', 'stars_L',
	      'gas_total_mass', 'gas_maxdens', 'gas_L', 'rvir', 'Mvir_dm', 'stars_center','snap_files']
    for field in fields: 
	if field in ['scale', 'stars_total_mass', 'stars_rhalf', 'gas_total_mass' ]:
	    galaxy_props[field] = np.array([])                
	else :
	    galaxy_props[field] = []


    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])

    ts = yt.DatasetSeries(new_snapfiles)

    for ds,snap_dir in zip(reversed(ts),np.flipud(new_snapfiles)):
        print "Getting galaxy props: ",  snap_dir

        ds.add_particle_filter('stars')
        
        dd = ds.all_data()
        ds.domain_right_edge = ds.arr(ds.domain_right_edge,'code_length')
        ds.domain_left_edge  = ds.arr(ds.domain_left_edge,'code_length')
        print ds.index.get_smallest_dx()

        #need to exit gracefully here if there's no stars.
        try:
            stars_pos_x = dd['stars', 'particle_position_x'].in_units('kpc')
            assert stars_pos_x.shape > 5
        except AttributeError,AssertionError:
            print "No star particles found, skipping: ", snap_dir
            continue


        scale = round(1.0/(ds.current_redshift+1.0),3)
        galaxy_props['scale'] = np.append(galaxy_props['scale'], scale)
    
        galaxy_props['snap_files'] = np.append(galaxy_props['snap_files'],snap_dir)


        print 'Determining center...'
        max_ndens_arr = find_center(dd, ds, cen_pos = ds.domain_center.in_units('kpc')[0].value[()], units = 'kpc')
        print '\tCenter = ', max_ndens_arr

        #Generate Sphere Selection
        print 'Determining virial radius...'
        rvir = find_rvirial(dd, ds, max_ndens_arr)
        print '\tRvir = ', rvir

        hc_sphere = ds.sphere(max_ndens_arr, rvir)

 
        galaxy_props['stars_maxndens'].append(max_ndens_arr.value)
        galaxy_props['rvir'] = np.append(galaxy_props['rvir'], rvir.value[()])
        galaxy_props['Mvir_dm'] = np.append(galaxy_props['Mvir_dm'], hc_sphere[('darkmatter', 'particle_mass')].in_units('Msun').sum().value[()])

		
        #Find Galaxy Properties
        galaxy_props = find_galaxyprops(galaxy_props, ds, hc_sphere, max_ndens_arr)


        del (hc_sphere)
        sys.stdout.flush()


    # Save galaxy props file
    galaxy_props_file = simname+'_galprops.npy'
    print '\nSuccessfully computed galaxy properties'
    print 'Saving galaxy properties to ', galaxy_props_file
    print
    np.save(galaxy_props_file, galaxy_props)  
