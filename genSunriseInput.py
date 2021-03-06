'''
Working from Miguel Rocha's script: genSunriseInput.py. This script
generates the cameras to use in Sunrise and makes some projection plots for
a few of the cameras. 

'''
import numpy as np
from numpy import *
import os, sys, argparse
from collections import OrderedDict
import time
#from blist import blist
import glob

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

    parser.add_argument('segments_random', nargs='?', default=None, help='Number of random cameras.')

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

    parser.add_argument('--no_export',action='store_true',
                        help='Do not export data to fits for Sunrise.') 

    parser.add_argument('--no_gas_p',action='store_true',
                        help='Shut off momentum of gas grid.')
    
    parser.add_argument('--email',default='gsnyder@stsci.edu',type=str,
                        help='email address for job notifications')
    

    args = vars(parser.parse_args())
    return args



def generate_cameras(normal_vector, seed = 0, distance=100.0, fov=50.0, mov_ang = 0., movie = False, segments_random=7, segments_fixed=4):
    '''
    Set camera positions and orientations
    '''
    from yt.utilities.orientation import Orientation

    print( "\nGenerating cameras")
    
    north = np.array([0.,1.,0.])
    orient = Orientation(normal_vector=normal_vector, north_vector=north)
    R=np.linalg.inv(orient.inv_mat)
    camera_set = OrderedDict([
            ['face',([0.,0.,1.],[0.,-1.,0],True)], #up is north=+y
            ['edge',([0.,1.,0.],[0.,0.,-1.],True)],#up is along z
            ['backface',([0.,0.,-1.],[0.,-1.,0],True)], #up is north=+y
            ['backedge',([0.,-1.,0.],[0.,0.,-1.],True)],#up is along z
            ['45',([0.,0.7071,0.7071],[0., 0., -1.],True)],
            ['Z-axis',([0.,0.,-1.],[0.,-1.,0],False)], #up is north=+y
            ['Y-axis',([0.,1.,0.],[0.,0.,-1.],False)],#up is along z
            ['X-axis',([1.,0.,0.],[0.,0.,-1.],False)],#up is along z
            ])  


    np.random.seed()
    ts_random = np.random.random(segments_random)*np.pi*2
    #ps_random = np.random.random(segments_random)*np.pi-np.pi/2.0
    #ps_random = np.random.random(segments_random)*np.pi
    ps_random = array([np.math.acos(2*np.random.random()-1) for i in arange(segments_random)])

    np.random.seed(seed)
    ts_fixed = np.random.random(segments_fixed)*np.pi*2
    #ps_fixed = np.random.random(segments_fixed)*np.pi-np.pi/2.0
    #ps_fixed = np.random.random(segments_fixed)*np.pi
    ps_fixed = array([np.math.acos(2*np.random.random()-1) for i in arange(segments_fixed)])




    ts = np.concatenate([ts_fixed, ts_random])
    ps = np.concatenate([ps_fixed, ps_random])


    for i,(theta, phi) in enumerate(zip(ts,ps)):
        print( theta, phi)
        pos = [np.cos(theta)*np.sin(phi),np.sin(theta)*np.sin(phi),np.cos(phi)]
        vc = [np.cos(np.pi/2 - theta)*np.sin(np.pi/2-phi),np.sin(np.pi/2 - theta)*np.sin(np.pi/2-phi),np.cos(np.pi/2-phi)]
        
        if i < segments_fixed:
            camera_set['Fixed_%03i'%(i)]=(pos,vc,False)
        else:
            camera_set['Random_%03i'%(i-segments_fixed)]=(pos,vc,False)

    if movie: camera_set['Movie_angle'] = ([0.,math.cos(mov_ang*180./pi),math.sin(2*180./pi)],[0.,0.,-1.],False)

    i=0  
    cameras = OrderedDict()
    for name,(normal,north,do_rot)  in camera_set.items():
        print( name, normal, north, do_rot)
        
        orient = Orientation(normal_vector=normal, north_vector=north)
        if do_rot:
            drot = R.copy()
        else:
            drot = np.identity(3)
        sunrise_pos = np.dot(orient.normal_vector, drot)
        if name == 'Z-axis':
            sunrise_up = np.asarray([0.0,1.0,0.0])
        elif name == 'X-axis':
            sunrise_up = np.asarray([0.0,1.0,0.0])
        elif name =='Y-axis':
            sunrise_up = np.asarray([0.0,0.0,1.0])
        else:
            sunrise_up  = normal_vector.copy()

            
        if np.all(np.abs(sunrise_up-sunrise_pos)<1e-3):
            sunrise_up[0] *= 0.5 
        sunrise_direction = -1.0*sunrise_pos
        sunrise_afov = 2.0*np.arctan((fov/2.0)/distance)
        norm = lambda x: x/np.sqrt(np.sum(x*x))
        if np.all(np.abs(norm(sunrise_up)-norm(sunrise_pos))<1e-3):
            sunrise_up[0]*=0.5
            sunrise_up = norm(sunrise_up)
        line = (distance*sunrise_pos, distance*sunrise_direction, sunrise_up,
                sunrise_afov, fov, distance) 
        cameras[name] = line
        i+=1

    print( "Successfully generated cameras\n")
    return cameras

def write_cameras(prefix, cameras):
    print( "Writing cameras to ",  prefix+'.cameras')
    fn = prefix + '.cameras'
    campos = ()
    for name,row in cameras.items():
        campos += (tuple(row[1])+tuple(row[0])+tuple(row[2])+tuple([row[3]]),)
    campos = np.array(campos)
    np.savetxt(fn, campos)   
    fn =  prefix+'.camnames'
    fh = open(fn,'w')
    fh.write('\n'.join([c for c in cameras.keys()]))
    fh.close()

def export_fits(ds, center, export_radius, prefix, star_particles, max_level=None, no_gas_p = False, form='VELA'):
    '''
    Convert the contents of a dataset to a FITS file format that Sunrise
    understands.
    '''
    import sunrise_octree_exporter

    print( "\nExporting data in %s to FITS for Sunrise"%ds.parameter_filename.split('/')[-1])
    
    filename = prefix+'.fits'
    center = center.in_units('kpc')
    width = export_radius.in_units('kpc')
    info = {}
    
    fle, fre, ile, ire, nrefined, nleafs, nstars, output, output_array = \
	                                                                 sunrise_octree_exporter.export_to_sunrise(ds, filename, star_particles,  center, width, max_level=max_level, grid_structure_fn = prefix+'_grid_struct.npy', no_gas_p = no_gas_p, form=form)
    
    info['export_ile']=ile
    info['export_ire']=ire
    info['export_fle']=fle
    info['export_fre']=fre
    
    
    info['export_center']=center.value
    info['export_radius']=width.value
    info['export_max_level']=max_level
    info['export_nstars']=nstars
    info['export_nrefined']=nrefined
    info['export_nleafs']=nleafs
    info['input_filename']=filename
    
    print( "Successfully generated FITS for snapshot %s"%ds.parameter_filename.split('/')[-1])
    print( info,'\n')
    #return info,  output, output_array
    return info  #output arrays not actually used later


def write_qsub_exporters(snapname,qsubfn,aname,args):

    if args['email'] is not None:
        en = args['email']
    else:
        en = 'gsnyder@stsci.edu'
        

    if args['format']=='ENZO':
        code='enzoToSunrise.py'
        group=args['group']
    elif args['format']=='VELA':
        code='genSunriseInput.py'
        group=None
        
    qsfo = open(qsubfn,'w')
    qsfo.write('#!/bin/bash\n')
    qsfo.write('#PBS -S /bin/bash\n')
    if group is not None:
        qsfo.write('#PBS -W group_list='+group+'\n')
    qsfo.write('#PBS -l select=1:ncpus=1:model=has\n')
    qsfo.write('#PBS -l walltime=04:00:00\n')
    qsfo.write('#PBS -q normal\n')
    qsfo.write('#PBS -N sunrise_export\n')
    qsfo.write('#PBS -M '+en+'\n')
    qsfo.write('#PBS -m abe\n')
    qsfo.write('#PBS -o sunrise_export_'+aname+'pbs.out\n')
    qsfo.write('#PBS -e sunrise_export_'+aname+'pbs.err\n')
    qsfo.write('#PBS -V\n\n')
    
    qsfo.write('python $VELAYTSUNRISE_CODE/'+code+' '+os.path.basename(snapname)+' --fov='+str(args['fov'])+' --format='+str(args['format'])+' > export_test_'+aname+'.out 2> export_test_'+aname+'.err\n')
    qsfo.close()
    

    return



if __name__ == "__main__":

    args = parse()

    import yt
    import sunrise_octree_exporter
    reload(sunrise_octree_exporter)
    
    
    print( args)
    print( args['no_export'], args['no_gas_p'])
    
    if args['snap_files'] is not None:
        snaps = [args['snap_files']]
    else:
        snaps = np.asarray(glob.glob("*.d"))

    if args['segments_random'] is not None:
        segments_random=int(args['segments_random'])
    else:
        segments_random=7
        
    print( "Generating Sunrise Input for: ", snaps)

    abssnap = os.path.abspath(snaps[0])

    
    assert os.path.lexists(abssnap)
    
    dirname = os.path.dirname(abssnap)
    simname = os.path.basename(dirname) #assumes directory name for simulation name
    print( "Simulation name:  ", simname)
    
    particle_headers = []
    particle_data = []
    stars_data = []
    new_snapfiles = []
    for sn in snaps:
        aname = sn.split('_')[-1].rstrip('.d')
        particle_headers.append('PMcrd'+aname+'.DAT')
        particle_data.append('PMcrs0'+aname+'.DAT')
        stars_data.append('stars_'+aname+'.dat')
        snap_dir = os.path.join(simname+'_'+aname+'_sunrise')
        
        print( "Sunrise directory: ", snap_dir)
        if not os.path.lexists(snap_dir):
            os.mkdir(snap_dir)        

        newf = os.path.join(snap_dir,sn)
        new_snapfiles.append(newf)
        if not os.path.lexists(newf):
            os.symlink(os.path.abspath(sn),newf)
            os.symlink(os.path.abspath(particle_headers[-1]),os.path.join(snap_dir,particle_headers[-1]))
            os.symlink(os.path.abspath(particle_data[-1]),os.path.join(snap_dir,particle_data[-1]))
            os.symlink(os.path.abspath(stars_data[-1]),os.path.join(snap_dir,stars_data[-1]))


    new_snapfiles = np.asarray(new_snapfiles)


    #exit()



    #gen_name, gal_name, snap_name, snaps  = 'VELA_v2', 'VELA27', 'VELA27_a0.370', '../data/VELA27_v2/a0.370/10MpcBox_csf512_a0.370.d'
        
    #snap_dir = '/Volumes/wd/yt_pipeline/Runs/%s/%s/%s'%(gen_name, gal_name, snap_name+'_sunrise')

    #assert os.path.exists(snap_dir), 'Snapshot directory %s not found'%snap_dir

        
        
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


    #ts = yt.DatasetSeries(new_snapfiles)
        
    # Loop over snapshot to generate cameras and projection plots, 
    # parallelization happens while generating the plots.
    for snapfile in new_snapfiles:

        aname = (os.path.basename(snapfile)).split('_')[-1].rstrip('.d')
        
        print( "Timestep name: ", aname)

        snap_dir = os.path.dirname(snapfile) #os.path.join(simname+'_'+aname+'_sunrise')

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
        gas_L 	= galprops['gas_L'][idx]

        
        try:
            L_sum = stars_L + gas_L
        except TypeError:
            L_sum = gas_L

        L = L_sum/np.sqrt(np.sum(L_sum*L_sum))
        

        #L_temp = array([0.229307690083501, 0.973325655982054, 0.00742635009091421]) #to Match with C Moody
        #This function is important for generating the cameras that we will be using
        try:
            cameras = generate_cameras(L, seed = seed, distance = cam_dist, fov = cam_fov, segments_random=segments_random)
        except np.linalg.linalg.LinAlgError:
            print( "Error in camera linear algebra: skipping")
            continue

        prefix = os.path.join(out_dir,simname+'_'+aname)
        write_cameras(prefix, cameras)
        sys.stdout.flush()

        qsubfn = 'export_'+aname+'.qsub'
        write_qsub_exporters(snapfile,qsubfn,aname,args)
        submitline = 'qsub '+qsubfn
        eaf.write(submitline+'\n')


    eaf.close()
    
    if args['no_export'] is True:
        print( "Skipping export stage, per command argument.")
        exit()

    print( "Continuing to export grids.")

    yt.enable_parallelism()
    
    ts = yt.DatasetSeries(new_snapfiles)

    # Send one snapshots to each processor to export 

    for ds in ts.piter():
        aname = (os.path.basename(ds._file_amr)).split('_')[-1].rstrip('.d')
        print( "Timestep name: ", aname)
        snap_dir = os.path.dirname(ds._file_amr)
        assert os.path.lexists(snap_dir)
        
        out_dir = os.path.join(snap_dir, 'input')
        assert os.path.lexists(out_dir)
        
        prefix = os.path.join(out_dir,simname+'_'+aname)
        
        snapfile = ds._file_amr
        
        if os.path.abspath(snapfile) not in galprops['snap_files']: continue
        idx = np.argwhere(galprops['snap_files']==os.path.abspath(snapfile))[0][0]
        
        scale = round(1.0/(ds.current_redshift+1.0),4)
        #idx = np.argwhere(galprops['scale'] == scale)[0][0]
	
        gal_center = galprops['stars_center'][idx]
        gal_center = ds.arr(gal_center, 'kpc')
        
        #export_radius = ds.arr(max(1.2*cam_dist, 1.2*cam_fov), 'kpc')
        export_radius = ds.arr(1.2*cam_fov, 'kpc')
        
        print( export_radius)
        
        export_info = export_fits(ds, gal_center, export_radius, 
                                  prefix, star_particles = 'stars', 
                                  max_level=max_level, no_gas_p = args['no_gas_p'])

        export_info['sim_name'] = simname
        export_info['scale'] = scale
        export_info['snap_file'] = snapfile
        
        export_info_file = prefix + '_export_info.npy' #galprops_file.replace('galprops', 'export_info')
        np.save(export_info_file, export_info)
        sys.stdout.flush()

    b = time.time()
    print( 'Final time in seconds: ', b - a)

































