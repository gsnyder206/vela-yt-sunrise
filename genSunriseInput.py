'''
Working from Miguel Rocha's script: genSunriseInput.py. This script
generates the cameras to use in Sunrise and makes some projection plots for
a few of the cameras. 

'''
import numpy as np
from numpy import *
import os, sys
from collections import OrderedDict
import time
from blist import blist
import glob


def generate_cameras(normal_vector, distance=100.0, fov=50.0, mov_ang = 0.):
    '''
    Set camera positions and orientations
    '''
    from yt.utilities.orientation import Orientation

    print "\nGenerating cameras"
    
    north = np.array([0.,1.,0.])
    orient = Orientation(normal_vector=normal_vector, north_vector=north)
    R=np.linalg.inv(orient.inv_mat)
    camera_set = OrderedDict([
            ['face',([0.,0.,1.],[0.,-1.,0],True)], #up is north=+y
            ['edge',([0.,1.,0.],[0.,0.,-1.],True)],#up is along z
            ['45',([0.,0.7071,0.7071],[0., 0., -1.],True)],
            ['Z-axis',([0.,0.,1.],[0.,-1.,0],False)], #up is north=+y
            ['Y-axis',([0.,1.,0.],[0.,0.,-1.],False)],#up is along z
            ['Movie_angle',([0.,math.cos(mov_ang*180./pi),math.sin(2*180./pi)],[0.,0.,-1.],False)],#up is along z                        
            ])  

    segments = 10
    #np.random.seed(0), we don't want seed here
    ts = np.random.random(segments)*np.pi*2
    ps = np.random.random(segments)*np.pi-np.pi/2.0
    for i,(theta, phi) in enumerate(zip(ts,ps)):
    	print theta, phi
    	#print theta*180./pi, phi*180./pi
        pos = [np.cos(theta),0.,np.sin(phi)]
        vc  = [np.cos(np.pi/2.-theta),0.,np.sin(np.pi/2.-phi)] 
        camera_set['Random_%03i'%(i)]=(pos,vc,False)
	ts = array([-0.00742641835397925])
	ps = array([-1.80216918439022])
	for i,(theta, phi) in enumerate(zip(ts,ps)):
		#print theta*180./pi, phi*180./pi
		pos = [np.cos(theta),0.,np.sin(phi)]
		vc  = [np.cos(np.pi/2.-theta),0.,np.sin(np.pi/2.-phi)] 
		camera_set['Random_match']=(pos,vc,False)


    i=0    
    cameras = OrderedDict()
    for name,(normal,north,do_rot)  in camera_set.iteritems():
        orient = Orientation(normal_vector=normal, north_vector=north)
        if do_rot:
            drot = R.copy()
        else:
            drot = np.identity(3)
        sunrise_pos = np.dot(orient.normal_vector, drot)
        sunrise_up  = L.copy()
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

    print "Successfully generated cameras\n"
    return cameras

def write_cameras(prefix, cameras):
    print "Writing cameras to ",  prefix+'.cameras'
    fn = prefix + '.cameras'
    campos = ()
    for name,row in cameras.iteritems():
        campos += (tuple(row[1])+tuple(row[0])+tuple(row[2])+tuple([row[3]]),)
    campos = np.array(campos)
    np.savetxt(fn, campos)   
    fn =  prefix+'.camnames'
    fh = open(fn,'w')
    fh.write('\n'.join([c for c in cameras.keys()]))
    fh.close()

def export_fits(ds, center, export_radius, prefix, star_particles, max_level=None):
	'''
	Convert the contents of a dataset to a FITS file format that Sunrise
	understands.
	'''

	print "\nExporting data in %s to FITS for Sunrise"%ds.parameter_filename.split('/')[-1]

	filename = prefix+'.fits'
	center = center.in_units('kpc')
	width = export_radius.in_units('kpc')
	info = {}

	fle, fre, ile, ire, nrefined, nleafs, nstars, output, output_array = \
	sunrise_octree_exporter.export_to_sunrise(ds, filename, star_particles,  center, width, max_level=max_level)

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

	print "Successfully generated FITS for snapshot %s"%ds.parameter_filename.split('/')[-1]
	print info,'\n'
	#return info,  output, output_array
        return info  #output arrays not actually used later


if __name__ == "__main__":
	#Should read these in from an initialization file
	#gen_name, gal_name, snap_name, snaps  = 'VELA_v2.1', 'VELA10', 'VELA10_a0.330', '../data/VELA10_v2.1/10MpcBox_csf512_a0.330.d'
	#gen_name, gal_name, snap_name, snaps  = 'VELA_v2', 'VELA27', 'VELA27_a0.560', '../data/VELA27_v2/a0.560/10MpcBox_csf512_a0.560.d'	
	#gen_name, gal_name, snap_name, snaps  = 'VELA_v2', 'VELA27', 'VELA27_a0.500', '../data/VELA27_v2/a0.500/10MpcBox_csf512_a0.500.d'
	
        if len(sys.argv)==2:
            snaps = np.asarray([sys.argv[1]])
        else:
            snaps = np.asarray(glob.glob("*.d"))



        print "Generating Sunrise Input for: ", snaps

        abssnap = os.path.abspath(snaps[0])
        assert os.path.lexists(abssnap)

        dirname = os.path.dirname(abssnap)
        simname = os.path.basename(dirname) #assumes directory name for simulation name
        print "Simulation name:  ", simname

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

                print "Sunrise directory: ", snap_dir
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
	import yt
	import sunrise_octree_exporter
	reload(sunrise_octree_exporter)


	cam_dist = 100000
	cam_fov  = 50
        max_level = None




	#ts = yt.DatasetSeries(new_snapfiles)

	# Loop over snapshot to generate cameras and projection plots, 
        # parallelization happens while generating the plots.
	for snapfile in new_snapfiles:

                aname = (os.path.basename(snapfile)).split('_')[-1].rstrip('.d')
        
                print "Timestep name: ", aname

                snap_dir = os.path.dirname(snapfile) #os.path.join(simname+'_'+aname+'_sunrise')

                print "Sunrise directory: ", snap_dir
                assert os.path.lexists(snap_dir)


                    
                galprops_file = simname+'_galprops.npy'

                out_dir = os.path.join(snap_dir,'input')
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
		cameras = generate_cameras(L, distance = cam_dist, fov = cam_fov)
                prefix = os.path.join(out_dir,simname+'_'+aname)
		write_cameras(prefix, cameras)
                sys.stdout.flush()



	ts = yt.DatasetSeries(new_snapfiles,limit_level=2)

	# Send one snapshots to each processor to export 

	for ds in ts.piter():
                aname = (os.path.basename(ds._file_amr)).split('_')[-1].rstrip('.d')
                print "Timestep name: ", aname
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

		print export_radius

		export_info = export_fits(ds, gal_center, export_radius, 
                                          prefix, star_particles = 'stars', 
                                          max_level=max_level)



		export_info['sim_name'] = simname
		export_info['scale'] = scale
                export_info['snap_file'] = snapfile

		export_info_file = prefix + '_export_info.npy' #galprops_file.replace('galprops', 'export_info')
		np.save(export_info_file, export_info)
                sys.stdout.flush()

	b = time.time()
	print 'Final time in seconds: ', b - a

































