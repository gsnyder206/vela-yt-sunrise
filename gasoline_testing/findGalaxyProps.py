'''
Working from Miguel Rocha's script: findGalaxyProps.py. Find the center of the galaxy
at the peak in the stellar number density. Generate galaxy properties.


'''

import sys
import os
import glob
import yt
import numpy as np
from numpy import *
import astropy
from astropy.cosmology import Planck13 as cosmo
#reload(yt)

def find_center(dd, ds, units = 'kpc', cen_pos = 10.e3, bin_width = 2.e4, del_pos = 200):
    '''
    find the center using the number density
    all lengths are in kpc
    returns ndarray of max_ndens_arr = ([cenx, ceny, cenz])
    '''
    units = 'kpc' 
    #stars_pos_x = dd['stars', 'particle_position_x'].in_units(units)
    #stars_pos_y = dd['stars', 'particle_position_y'].in_units(units)
    #stars_pos_z = dd['stars', 'particle_position_z'].in_units(units)

    stars_pos_x = dd['Stars', 'Coordinates'][:,0].in_units('kpc')
    stars_pos_y = dd['Stars', 'Coordinates'][:,1].in_units('kpc')
    stars_pos_z = dd['Stars', 'Coordinates'][:,2].in_units('kpc')




    star_pos = [stars_pos_x.value, stars_pos_y.value, stars_pos_z.value]

    min_pos = cen_pos - bin_width
    max_pos = cen_pos + bin_width
    bins = [arange(min_pos,max_pos,del_pos), arange(min_pos,max_pos,del_pos), arange(min_pos,max_pos,del_pos)]



    H, edges = histogramdd(star_pos, bins = bins)
    max_ndens_index = unravel_index(H.argmax(), H.shape)


    print min_pos, max_pos, del_pos


    max_ndens_loc = array([(edges[0][max_ndens_index[0]] + edges[0][max_ndens_index[0]+1])/2., 
    					   (edges[1][max_ndens_index[1]] + edges[1][max_ndens_index[1]+1])/2.,
    					   (edges[2][max_ndens_index[2]] + edges[2][max_ndens_index[2]+1])/2.])




    max_ndens_arr = ds.arr([max_ndens_loc[0], max_ndens_loc[1], max_ndens_loc[2]], units)



    #end of First pass
    print('\tDone with coarse pass searching for center, moving to fine pass')
    	

    bin_width = 200
    del_pos = 0.5

    min_pos_x = float(max_ndens_arr[0]) - bin_width
    max_pos_x = float(max_ndens_arr[0]) + bin_width

    min_pos_y = float(max_ndens_arr[1]) - bin_width
    max_pos_y = float(max_ndens_arr[1]) + bin_width

    min_pos_z = float(max_ndens_arr[2]) - bin_width
    max_pos_z = float(max_ndens_arr[2]) + bin_width


    bins = [arange(min_pos_x,max_pos_x,del_pos), arange(min_pos_y,max_pos_y,del_pos), arange(min_pos_z,max_pos_z,del_pos)]

    H, edges = histogramdd(star_pos, bins = bins)
    max_ndens_index = unravel_index(H.argmax(), H.shape)

    max_ndens_loc = array([(edges[0][max_ndens_index[0]] + edges[0][max_ndens_index[0]+1])/2., 
    					   (edges[1][max_ndens_index[1]] + edges[1][max_ndens_index[1]+1])/2.,
    					   (edges[2][max_ndens_index[2]] + edges[2][max_ndens_index[2]+1])/2.])

    max_ndens_arr = ds.arr([max_ndens_loc[0], max_ndens_loc[1], max_ndens_loc[2]], units)


    return max_ndens_arr

def find_rvirial(dd, ds, center, start_rad = 0, delta_rad_coarse = 20, delta_rad_fine = 1, rad_units = 'kpc'):
    vir_check = 0
    r0 = ds.arr(start_rad, rad_units)
    critical_density = cosmo.critical_density(ds.current_redshift).value   #is in g/cm^3
    max_ndens_arr=center

    while True:
        r0_prev = r0
        r0 = r0_prev + ds.arr(delta_rad_coarse, rad_units)
        v_sphere = ds.sphere(max_ndens_arr, r0)
        dark_mass 	= v_sphere[('DarkMatter', 'Mass')].in_units('Msun').sum()	
        rho_internal = dark_mass.in_units('g')/((r0.in_units('cm'))**3.*(pi*4/3.))
        print rho_internal, r0, dark_mass
        if rho_internal < 200*ds.arr(critical_density,'g')/ds.arr(1.,'cm')**3.:
            #now run fine test
            print('\tNow running fine search on the virial radius')
            r0 = r0_prev
            while True:
                r0 += ds.arr(delta_rad_fine, rad_units)
                v_sphere = ds.sphere(max_ndens_arr, r0)
                dark_mass 	= v_sphere[('DarkMatter', 'Mass')].in_units('Msun').sum()	
                rho_internal = dark_mass.in_units('g')/((r0.in_units('cm'))**3.*(pi*4/3.))
                print rho_internal, r0
                if rho_internal < 200*ds.arr(critical_density,'g')/ds.arr(1.,'cm')**3.:
                    rvir = r0
                    print dark_mass.in_units('Msun'), rvir.in_units('kpc')
                    return rvir

def find_hist_center(positions, masses):
    '''
    Find the center of a particle distribution by interactively refining 
    a mass weighted histogram
    '''
    pos = np.array(positions)
    masses = np.array(masses)
    if len(pos) == 0: 
        return None
    mass_current = masses
    old_center = np.array([0,0,0])
    refined_pos = pos.copy()
    refined_mas = mass_current.copy()
    refined_dist = 1e20
    nbins=3
    center = None

    dist = lambda x,y:np.sqrt(np.sum((x-y)**2.0))
    dist2 = lambda x,y:np.sqrt(np.sum((x-y)**2.0,axis=1))

    j=0
    while len(refined_pos)>1e1 or j==0: 
        table,bins=np.histogramdd(refined_pos, bins=nbins, weights=refined_mas)
        bin_size = min((np.max(bins,axis=1)-np.min(bins,axis=1))/nbins)
        centeridx = np.where(table==table.max())
        le = np.array([bins[0][centeridx[0][0]],
                       bins[1][centeridx[1][0]],
                       bins[2][centeridx[2][0]]])
        re = np.array([bins[0][centeridx[0][0]+1],
                       bins[1][centeridx[1][0]+1],
                       bins[2][centeridx[2][0]+1]])
        center = 0.5*(le+re)
        refined_dist = dist(old_center,center)
        old_center = center.copy()
        idx = dist2(refined_pos,center)<bin_size
        refined_pos = refined_pos[idx]
        refined_mas = refined_mas[idx]
        j+=1    

    return center

def find_shapes(center, pos, ds, nrad=10, rmax=None):
    '''
    Find the shape of the given particle distribution at nrad different 
    radii, spanning from 0.1*rmax to rmax. 
    rmax = max(r(pos)) if not given.
    '''

    print('Starting shape calculation')

    units = center.units
    center = center.value

    try:
        pos = np.array([pos[:,0] - center[0],
                        pos[:,1] - center[1],
                        pos[:,2] - center[2]]).transpose()
        pos = ds.arr(pos, units)
        pos = pos.in_units(units).value
        r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
    except IndexError: # no stars found
        pos = np.array([])   

    if len(pos) > 1: 
        if not rmax: rmax = r.max()
        radii = np.linspace(0.1*rmax, rmax, nrad)
    else:
        radii = np.array([])
    
    c_to_a = np.empty(radii.size)     
    b_to_a = np.empty(radii.size)
    axes = []
        
    for i,r in enumerate(radii):
        # get shapes
        try:
            axis_out = axis_ratios(pos, r, axes_out=True, fix_volume = False)
            c_to_a[i] = axis_out[0][0]
            b_to_a[i] = axis_out[0][1]
            axes.append(axis_out[1])
        except UnboundLocalError:
            print( 'Not enough particles to find shapes at r = %g in snapshot %s'%(r, ds.parameter_filename ))
            b_to_a[i] = c_to_a[i] = None
            axes.append([])
    
    return radii, c_to_a, b_to_a, axes        
       
def L_crossing(x, y, z, vx, vy, vz, weight, center):
    x, y, z = x-center[0], y-center[1],z-center[2]
    cx, cy, cz = y*vz - z*vy, z*vx - x*vz, x*vy - y*vx
    lx, ly, lz = [np.sum(l * weight) for l in [cx, cy, cz]]
    L = np.array([lx, ly, lz])
    L /= np.sqrt(np.sum(L*L))
    return L

def find_galaxyprops(galaxy_props, ds, hc_sphere, max_ndens_arr):

        print( 'Determining stellar and gas mass...')
        # Get total stellar mass 
        stars_mass = hc_sphere[('Stars', 'Mass')].in_units('Msun')
        stars_total_mass = stars_mass.sum().value[()]
        galaxy_props['stars_total_mass'] = np.append(galaxy_props['stars_total_mass'], stars_total_mass)

        # Get total mass of gas
        gas_mass = hc_sphere[('Gas', 'Mass')].in_units('Msun')
        gas_total_mass = gas_mass.sum().value[()]
        galaxy_props['gas_total_mass'] = np.append(galaxy_props['gas_total_mass'], 
                                                   gas_total_mass)
        print( '\tlog Mgas/Msun = ', log10(gas_total_mass))
        print( '\tlog M*/Msun = ', log10(stars_total_mass))

        '''
        print( 'Determining location of max stellar density...')
        # Get max density of stars (value, location)
        stars_maxdens = hc_sphere.quantities.max_location(('deposit', 'stars_cic'))
        stars_maxdens_val = stars_maxdens[0].in_units('Msun/kpc**3').value[()]

        print( stars_maxdens)
        #difference bt yt-3.2.3 and yt-3.3dev: stars_maxdens has different # elements; this works for both
        stars_maxdens_loc = np.array([stars_maxdens[-3].in_units('kpc').value[()], 
                                      stars_maxdens[-2].in_units('kpc').value[()], 
                                      stars_maxdens[-1].in_units('kpc').value[()]])
        galaxy_props['stars_maxdens'].append((stars_maxdens_val, stars_maxdens_loc))
        print( '\t Max Stellar Density = ', stars_maxdens_loc)
        '''


        print( 'Determining location of max gas density...')
        # Get max density of gas
        gas_maxdens = hc_sphere.quantities.max_location(('Gas', 'Density'))
        gas_maxdens_val = gas_maxdens[0].in_units('Msun/kpc**3').value[()]
        gas_maxdens_loc = np.array([gas_maxdens[-3].in_units('kpc').value[()], 
                                    gas_maxdens[-2].in_units('kpc').value[()], 
                                    gas_maxdens[-1].in_units('kpc').value[()]])
        galaxy_props['gas_maxdens'].append((gas_maxdens_val, gas_maxdens_loc)) 
        print( '\t Max Gas Density = ', gas_maxdens_loc)



        print( 'Determining refined histogram center of stars...')
        #---Need to Check these--#
        # Get refined histogram center of stars
        stars_pos_x = dd['Stars', 'Coordinates'][:,0].in_units('kpc')
        stars_pos_y = dd['Stars', 'Coordinates'][:,1].in_units('kpc')
        stars_pos_z = dd['Stars', 'Coordinates'][:,2].in_units('kpc')





        '''
        stars_pos = np.array([stars_pos_x, stars_pos_y, stars_pos_z]).transpose()
        stars_hist_center = find_hist_center(stars_pos, stars_mass)
        galaxy_props['stars_hist_center'].append(stars_hist_center)
        print( '\t Refined histogram center of stars = ', stars_hist_center)
        '''
        
        print( 'Computing stellar density profile...')

        # Get stellar density profile
        sc_sphere_r = 0.1
        ssphere_r = sc_sphere_r*hc_sphere.radius
        while ssphere_r < ds.index.get_smallest_dx():
                ssphere_r = 2.0*ssphere_r
        sc_sphere =  ds.sphere(max_ndens_arr, ssphere_r)
        '''
        try:
                p_plot = yt.ProfilePlot(sc_sphere, 'radius', 'Stars_Mass', n_bins=50, weight_field=None, accumulation=True)
                p_plot.set_unit('radius', 'kpc')
                p_plot.set_unit('Stars_Mass', 'Msun')
                p = p_plot.profiles[0]

                radii, smass = p.x.value, p['Stars_Mass'].value 
                rhalf = radii[smass >= 0.5*smass.max()][0]
        except (IndexError, ValueError): # not enough stars found
                radii, smass = None, None 
                rhalf = None
        galaxy_props['stars_rhalf'] = np.append(galaxy_props['stars_rhalf'], rhalf)
        galaxy_props['stars_mass_profile'].append((radii, smass))       

        print( '\tStars half-light radius = ', rhalf)

        '''
        print( 'Determining center of mass within 15 kpc of the galaxy...')

        '''
        # Get center of mass of stars
        gal_sphere = ds.sphere(max_ndens_arr, (15, 'kpc'))
        stars_pos_x = dd['Stars', 'Coordinates'][:,0].in_units('kpc')
        stars_pos_y = dd['Stars', 'Coordinates'][:,1].in_units('kpc')
        stars_pos_z = dd['Stars', 'Coordinates'][:,2].in_units('kpc')

        gal_stars_mass = gal_sphere[('Stars', 'Mass')].in_units('Msun')
        gal_total_mass = gal_stars_mass.sum().value[()]

        stars_com = np.array([np.dot(stars_pos_x, gal_stars_mass)/gal_total_mass, 
                              np.dot(stars_pos_y, gal_stars_mass)/gal_total_mass, 
                              np.dot(stars_pos_z, gal_stars_mass)/gal_total_mass])
        galaxy_props['stars_com'].append(stars_com)
        print( '\tCenter of mass = ', stars_com)
        '''





        print( 'Setting stars center...')
        # Define center of stars
        center = 'maxndens'
        if center == 'max_dens':
                stars_center = stars_maxdens_loc
        elif center == 'com':
                stars_center = stars_com
        elif center == 'maxndens':
                stars_center = max_ndens_arr
        else: 
                stars_center = stars_hist_center

        stars_center = ds.arr(stars_center, 'kpc')
        galaxy_props['stars_center'].append(stars_center)
        print( '\tStars Center = ', stars_center)



        # Get angular momentum of stars
        try:
                x, y, z = [sc_sphere[('Stars', 'particle_position_%s'%s)] for s in 'xyz'] 
                vx, vy, vz = [sc_sphere[('Stars', 'particle_velocity_%s'%s)] for s in 'xyz'] 
                mass = sc_sphere[('Stars', 'particle_mass')]
                try:
                        metals = sc_sphere[('Stars', 'metallicity')]
                        stars_L = L_crossing(x, y, z, vx, vy, vz, mass*metals, sc_sphere.center)
                except:
                        stars_L = L_crossing(x, y, z, vx, vy, vz, mass, sc_sphere.center)
                
        except IndexError: # no stars found
                stars_L = [None, None, None]
                print("No stars exception")

        galaxy_props['stars_L'].append(stars_L)
        del(sc_sphere)


        # Get angular momentum of gas
        gas_center = ds.arr(gas_maxdens_loc, 'kpc')
        print ssphere_r, gas_center
        gc_sphere =  ds.sphere(gas_center, ssphere_r)
        x = gc_sphere[('index', 'x')]
        x, y, z = [gc_sphere[('index', '%s'%s)] for s in 'xyz'] 
        cell_volume = gc_sphere[('index', 'cell_volume')]
        density=gc_sphere[('gas', 'density')]
        vx, vy, vz = [gc_sphere[('gas', 'velocity_%s'%s)] for s in 'xyz'] 
        gas_L = L_crossing(x, y, z, vx, vy, vz, density*cell_volume**3., gc_sphere.center)

        '''
        try:
                #for VELA runs
                vx, vy, vz = [gc_sphere[('gas', 'momentum_%s'%s)] for s in 'xyz'] # momentum density
                #metals = gc_sphere[('gas', 'metal_ia_density')] + gc_sphere[('gas', 'metal_ii_density')]
                gas_L = L_crossing(x, y, z, vx, vy, vz, metals*cell_volume**2, gc_sphere.center)


        except:
                #for enzo runs
                density=gc_sphere[('gas', 'density')]
                vx, vy, vz = [gc_sphere[('gas', 'velocity_%s'%s)] for s in 'xyz'] 
                metals=gc_sphere[('gas', 'metal_density')]
                gas_L = L_crossing(x, y, z, density*vx, density*vy, density*vz, density*cell_volume**2, gc_sphere.center)
        '''
        galaxy_props['gas_L'].append(gas_L)
        del(gc_sphere)

        return galaxy_props



if __name__ == "__main__":
	#Should read these in from an initialization file
	#gen_name, gal_name, snap_name, snaps  = 'VELA_v2', 'VELA27', 'VELA27_a0.370', '../data/VELA27_v2/a0.370/10MpcBox_csf512_a0.370.d'

	#snap_dir = '/Volumes/wd/yt_pipeline/Runs/%s/%s/%s'%(gen_name, gal_name, snap_name+'_sunrise')

	#if not os.path.isdir(snap_dir):
	#	os.system('mkdir '+'/Volumes/wd/yt_pipeline/Runs/%s/%s'%(gen_name, gal_name))		
	#	os.system('mkdir '+snap_dir)
	#	os.system('mkdir '+snap_dir+'/input')


	#assert os.path.exists(snap_dir), 'Snapshot directory %s not found'%snap_dir
	
        if len(sys.argv)==2:
            snaps = np.asarray([sys.argv[1]])
        else:
            snaps = np.sort(np.asarray(glob.glob("*.01024")))



        print( "Calculating Galaxy Props for: ", snaps)

        abssnap = os.path.abspath(snaps[0])
        assert os.path.lexists(abssnap)

        dirname = os.path.dirname(abssnap)
        simname = os.path.basename(dirname) #assumes directory name for simulation name
        print( "Simulation name:  ", simname)

        particle_headers = []
        particle_data = []
        stars_data = []
        new_snapfiles = []
        '''
        for sn in snaps:
                aname = sn.split('_')[-1].rstrip('.d')
                particle_headers.append('PMcrd'+aname+'.DAT')
                particle_data.append('PMcrs0'+aname+'.DAT')
                stars_data.append('stars_'+aname+'.dat')
                snap_dir = os.path.join(simname+'_'+aname+'_sunrise')
                yt_fig_dir = snap_dir+'/yt_projections'
                print( "Sunrise directory: ", snap_dir)
                #if not os.path.lexists(snap_dir):
                #    os.mkdir(snap_dir)        
                #if not os.path.lexists(yt_fig_dir):
                #    os.mkdir(yt_fig_dir)        

                newf = os.path.join(snap_dir,sn)
                new_snapfiles.append(newf)
                if not os.path.lexists(newf):
                        os.symlink(os.path.abspath(sn),newf)
                        os.symlink(os.path.abspath(particle_headers[-1]),os.path.join(snap_dir,particle_headers[-1]))
                        os.symlink(os.path.abspath(particle_data[-1]),os.path.join(snap_dir,particle_data[-1]))
                        os.symlink(os.path.abspath(stars_data[-1]),os.path.join(snap_dir,stars_data[-1]))
        '''
        new_snapfiles = snaps

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

        #ts = yt.DatasetSeries(new_snapfiles)

        cosmology_parameters = {'current_redshift': 0.0,
                                'omega_lambda': 0.728,
                                'omega_matter': 0.272,
                                'hubble_constant': 0.702}
        print new_snapfiles
        ds = yt.load(new_snapfiles[0])
        #for ds,snap_dir in zip(reversed(ts, cosmology_parameters=cosmology_parameters),np.flipud(new_snapfiles)):

        dd = ds.all_data()
        #ds.domain_right_edge = ds.arr(ds.domain_right_edge,'code_length')
        #ds.domain_left_edge  = ds.arr(ds.domain_left_edge,'code_length')
        #print( ds.index.get_smallest_dx())

        #need to exit gracefully here if there's no stars.
        try:
            stars_pos_x = dd['stars', 'Coordinates'][:,0].in_units('kpc')
            assert stars_pos_x.shape[0] > 5
        except:
            pass

        scale = round(1.0/(ds.current_redshift+1.0),3)
        galaxy_props['scale'] = np.append(galaxy_props['scale'], scale)



        print( 'Determining center...')
        #max_ndens_arr = find_center(dd, ds, cen_pos = ds.domain_center.in_units('kpc')[0].value[()], units = 'kpc')
        max_ndens_arr = find_center(dd, ds, cen_pos = 0., units = 'kpc')
        print( '\tCenter = ', max_ndens_arr)


        #Generate Sphere Selection
        print( 'Determining virial radius...')
        rvir = find_rvirial(dd, ds, max_ndens_arr)
        print( '\tRvir = ', rvir)


        hc_sphere = ds.sphere(max_ndens_arr, rvir)
        

        galaxy_props['stars_maxndens'].append(max_ndens_arr.value)
        galaxy_props['rvir'] = np.append(galaxy_props['rvir'], rvir.value[()])
        galaxy_props['Mvir_dm'] = np.append(galaxy_props['Mvir_dm'], hc_sphere[('DarkMatter', 'Mass')].in_units('Msun').sum().value[()])


        #Find Galaxy Properties
        galaxy_props = find_galaxyprops(galaxy_props, ds, hc_sphere, max_ndens_arr)

        del (hc_sphere)
        sys.stdout.flush()


        # Save galaxy props file
        galaxy_props_file = '/u/rcsimons/gasoline_testing/galprops.npy'
        print( '\nSuccessfully computed galaxy properties')
        print( 'Saving galaxy properties to ', galaxy_props_file)
        np.save(galaxy_props_file, galaxy_props)  






























