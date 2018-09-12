import numpy as np
from numpy import *
import os, sys, argparse
from collections import OrderedDict
import time
#from blist import blist
import glob

import genSunriseInput as gSI


def _stars(pfilter, data):
    return data[(pfilter.filtered_type, "particle_type")] == 2

#this gets dark matter particles in zoom region only
def _darkmatter(pfilter, data):
    return data[(pfilter.filtered_type, "particle_type")] == 4



def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''\
                                Generate the cameras to use in Sunrise and make projection plots
                                of the data for some of these cameras. Then export the data within
                                the fov to a FITS file in a format that Sunrise understands.
                                ''')
 
    parser.add_argument('snap_files', nargs='?', default=None, help='Snapshot files to be analyzed.')
    
    #parser.add_argument('-s', '--snap_base', default='10MpcBox_csf512_',
    #                    help='Base of the snapshots file names.') 

    parser.add_argument('-d', '--distance', default=100000, type=float,
                        help='Distance between cameras and the center of the galaxy (in [kpc]).')

    parser.add_argument('-f', '--fov', default=50, type=float,
                        help='Field of view of the cameras at the image plane (in [kpc]).')
    
    #I think many of these aren't needed for now.
    '''
    parser.add_argument('--star_particles', default='stars',
                        help='The name given to the star particles in the yt dataset '\
                        '(as printed in ds.field_list or ds.derived_field_list).') 

    parser.add_argument('--dm_particles', default='darkmatter',
                        help='The name given to the dark matter particles in the yt dataset '\
                        '(as printed in ds.field_list or ds.derived_field_list).') 

    parser.add_argument( '--galprops_file', default='sim_dir/analysis/catalogs/*_galaxy_props.npy',
                        help='File containing the galaxy properties. A python dictionary is expected '\
                             'as generated by findGalaxyProps.py.')

    parser.add_argument('--cams_to_plot', nargs='+', default=['face','edge','45'],
                        help='Cameras for which to make slice and projection plots ')

    parser.add_argument( '--max_level', default=None, type=int,
                         help='Max level to refine when exporting the oct-tree structure.')

    parser.add_argument('--out_dir',default='sim_dir/analysis/sunrise_analysis/',
                        help='Directory where the output will be placed. A sub directory will be created '\
                            'for each snapshot') 

    parser.add_argument('--rockstar_out_dir', default='sim_dir/analysis/rockstar_output/',
                        help='Directory where to find the rockstar output, used to annotate halos '\
                            'on plots. If not found or set to None halos will not be annotated.')
    '''
    #parser.add_argument('--no_plots',action='store_true',
    #                    help='Do not generate projection plots.') 

    parser.add_argument('--format',default='ENZO',type=str,
                        help='Simulation type (ENZO or VELA)') 
    
    parser.add_argument('--no_export',action='store_true',
                        help='Do not export data to fits for Sunrise.') 

    parser.add_argument('--no_gas_p',action='store_true',
                        help='Shut off momentum of gas grid.')
    
    #parser.add_argument('--email',default='gsnyder@stsci.edu',type=str,
    #                    help='email address for job notifications')
    parser.add_argument('--email',default='rsimons@jhu.edu',type=str,
                    help='email address for job notifications')

    parser.add_argument('--group',default='s1698',type=str,
                        help='group for compute charge')
    

    args = vars(parser.parse_args())
    return args



if __name__ == "__main__":

    args = parse()

    import yt
    import sunrise_octree_exporter
    #reload(sunrise_octree_exporter)
    
    
    print( args)
    print( args['no_export'], args['no_gas_p'], args['group'])

    if args['snap_files'] is not None:
        snaps = np.asarray( [args['snap_files']])
        form = args['format'] #VELA or ENZO
        if form=='ENZO' and snaps.shape[0]==1:
            snaps=np.asarray([snaps[0]+'/'+snaps[0]]) #handle the case of single-file input generated by scripts
    else:
        try:
            snaps = np.sort(np.asarray(glob.glob("*.d")))  #VELA format
            form='VELA'
            assert snaps.shape[0] > 0
        except AssertionError as e:
            snaps = np.sort(np.asarray(glob.glob("RD????/RD????")))  #ENZO format a list of snapshots in separate directories
            form='ENZO'

    args['format']=form

    assert snaps.shape[0] > 0

    print( "Generating Sunrise Inputs for "+form+": ", snaps)

    abssnap = os.path.abspath(snaps[0])
    assert os.path.lexists(abssnap)

    if form=='VELA':
        dirname = os.path.dirname(abssnap)
    elif form=='ENZO':
        dirname = os.path.dirname(os.path.dirname(abssnap))
        
    simname = os.path.basename(dirname) #assumes directory name for simulation name
    
    print( "Simulation name:  ", simname)
    
    particle_headers = []
    particle_data = []
    stars_data = []
    new_snapfiles = []
    anames=[]
    snapdirs=[]
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

        anames.append(aname)
        snapdirs.append(snap_dir)
        
        print( "Sunrise directory: ", snap_dir)
        assert os.path.lexists(snap_dir)

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


    #exit()

        
        
    a = time.time()



    cam_dist = 100000
    if args['fov'] is not None:
        cam_fov = float(args['fov'])
    else:
        cam_fov  = 50.0

    max_level = None
    seed = 0

    export_all_fn = 'export_all.sh'
    eaf = open(export_all_fn,'w')


        
    # Loop over snapshot to generate cameras and projection plots, 
    # parallelization happens while generating the plots.
    for snapfile,aname,snap_dir in zip(new_snapfiles,anames,snapdirs):
        
        print( "Timestep name: ", aname)

        print( "Sunrise directory: ", snap_dir)
        assert os.path.lexists(snap_dir)


                    
        galprops_file = simname+'_galprops.npy'

        out_dir = os.path.join(snap_dir,'input')
        print( os.path.lexists(out_dir))
        if not os.path.lexists(out_dir):
            os.mkdir(out_dir)

        galprops = np.load(galprops_file)[()]


        if os.path.abspath(snapfile) not in galprops['snap_files']: continue
        idx = np.argwhere(galprops['snap_files']==os.path.abspath(snapfile))[0][0]

        #scale = round(1.0/(ds.current_redshift+1.0),4)
        #if scale not in galprops['scale']: continue
        #idx = np.argwhere(galprops['scale'] == scale)[0][0]

        stars_L = galprops['stars_L'][idx]
        gas_L = galprops['gas_L'][idx]

        
        try:
            L_sum = stars_L + gas_L
        except TypeError:
            L_sum = gas_L

        L = L_sum/np.sqrt(np.sum(L_sum*L_sum))
        

        #L_temp = array([0.229307690083501, 0.973325655982054, 0.00742635009091421]) #to Match with C Moody
        #This function is important for generating the cameras that we will be using
        try:
            cameras = gSI.generate_cameras(L, seed = seed, distance = cam_dist, fov = cam_fov, segments_fixed=1, segments_random=1)
        except np.linalg.linalg.LinAlgError:
            print( "Error in camera linear algebra: skipping")
            continue

        prefix = os.path.join(out_dir,simname+'_'+aname)
        gSI.write_cameras(prefix, cameras)
        sys.stdout.flush()

        qsubfn = 'export_'+aname+'.qsub'
        gSI.write_qsub_exporters(snapfile,qsubfn,aname,args)
        submitline = 'qsub '+qsubfn
        eaf.write(submitline+'\n')
        

    eaf.close()
    
    if args['no_export'] is True:
        print( "Skipping export stage, per command argument.")
        exit()

    print( "Continuing to export grids.")

    if form=='VELA':
        starfield='stars'
    elif form=='ENZO':
        starfield='stars'

    
    yt.enable_parallelism()
    
    ts = yt.DatasetSeries(new_snapfiles)

    # Send one snapshots to each processor to export 
    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])

    
    for ds in ts.piter():
        if form=='VELA':
            aname = (os.path.basename(ds._file_amr)).split('_')[-1].rstrip('.d')
            print( "Timestep name: ", aname)
            snap_dir = os.path.dirname(ds._file_amr)
            assert os.path.lexists(snap_dir)
            
            out_dir = os.path.join(snap_dir, 'input')
            assert os.path.lexists(out_dir)
            
            prefix = os.path.join(out_dir,simname+'_'+aname)
        
            snapfile = ds._file_amr
        elif form=='ENZO':
            aname=ds.basename
            snap_dir=os.path.join(ds.fullpath,simname+'_'+aname+'_sunrise')
            out_dir=os.path.join(snap_dir,'input')
            prefix=os.path.join(out_dir,simname+'_'+aname)
            snapfile=os.path.join(ds.fullpath,aname)
            ds.add_particle_filter('stars')

            
        if os.path.abspath(snapfile) not in galprops['snap_files']: continue
        idx = np.argwhere(galprops['snap_files']==os.path.abspath(snapfile))[0][0]
        
        scale = round(1.0/(ds.current_redshift+1.0),4)
        #idx = np.argwhere(galprops['scale'] == scale)[0][0]
	
        gal_center = galprops['stars_center'][idx]
        gal_center = ds.arr(gal_center, 'kpc')
        
        #export_radius = ds.arr(max(1.2*cam_dist, 1.2*cam_fov), 'kpc')
        export_radius = ds.arr(1.2*cam_fov, 'kpc')
        
        print( export_radius)
        
        export_info = gSI.export_fits(ds, gal_center, export_radius, 
                                      prefix, star_particles = starfield, 
                                      max_level=max_level, no_gas_p = args['no_gas_p'], form=form)
        

        
        export_info['sim_name'] = simname
        export_info['scale'] = scale
        export_info['snap_file'] = snapfile

        export_info_file = prefix + '_export_info.npy' #galprops_file.replace('galprops', 'export_info')
        np.save(export_info_file, export_info)
        sys.stdout.flush()

    b = time.time()
    print( 'Final time in seconds: ', b - a)


