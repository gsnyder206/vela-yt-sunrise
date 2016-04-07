"""
Code to export from yt to Sunrise


"""


import time
import numpy as np
from numpy import *
import astropy
import astropy.io.fits as pyfits
import yt
from yt.funcs import get_pbar
from blist import blist
from yt.funcs import *
#import sunrise_export_gfs.octree_to_depthFirstHilbert_GFS as depthFirstHilbert


def export_to_sunrise(ds, fn, star_particle_type, fc, fwidth, nocts_wide=None, 
    debug=False,ad=None,max_level=None,  **kwargs):


	r"""Convert the contents of a dataset to a FITS file format that Sunrise
	understands.

	This function will accept a dataset, and from that dataset
	construct a depth-first octree containing all of the data in the parameter
	file.  This octree will be written to a FITS file.  It will probably be
	quite big, so use this function with caution!  Sunrise is a tool for
	generating synthetic spectra, available at
	http://sunrise.googlecode.com/ .

	Parameters
	----------
	ds : `Dataset`
	   The dataset to convert.
	fn : string
	   The filename of the output FITS file.
	fc : array
	   The center of the extraction region
	fwidth  : array  
	   Ensure this radius around the center is enclosed
	   Array format is (nx,ny,nz) where each element is floating point
	   in unitary position units where 0 is leftmost edge and 1
	   the rightmost. 

	Notes
	-----

	Note that the process of generating simulated images from Sunrise will
	require substantial user input; see the Sunrise wiki at
	http://sunrise.googlecode.com/ for more information.

	"""
	'''
	fc = fc.in_units('code_length').value
	fwidth = fwidth.in_units('code_length').value
	Nocts_root = ds.domain_dimensions/2
	'''

	fc = fc.in_units('code_length').value
	fwidth = fwidth.in_units('code_length').value
	Nocts_root = ds.domain_dimensions/2

	#we must round the dle,dre to the nearest root grid cells

	ile,ire,super_level,nocts_wide=  round_nocts_wide(Nocts_root,fc-fwidth,fc+fwidth,nwide=nocts_wide)
	assert np.all((ile-ire)==(ile-ire)[0])
	print "rounding specified region:"
	print "from [%1.5f %1.5f %1.5f]-[%1.5f %1.5f %1.5f]"%(tuple(fc-fwidth)+tuple(fc+fwidth))
	print "to (integer)   [%07i %07i %07i]-[%07i %07i %07i]"%(tuple(ile)+tuple(ire))
	assert(len(np.unique(ds.domain_width)) == 1)
	domain_width = ds.domain_width[0]
	fle,fre = ile*domain_width/Nocts_root, ire*domain_width/Nocts_root
	print "to (float)  [%1.5f %1.5f %1.5f]-[%1.5f %1.5f %1.5f]"%(tuple(fle)+tuple(fre))

	#Create a list of the star particle properties in PARTICLE_DATA
    #Include ID, parent-ID, position, velocity, creation_mass, 
    #formation_time, mass, age_m, age_l, metallicity, L_bol
	particle_data,nstars = prepare_star_particles(ds,star_particle_type,fle=fle,fre=fre, ad=ad,**kwargs)

	#Create the refinement depth-first hilbert octree structure
    #For every leaf (not-refined) oct we have a column n OCTDATA
    #Include mass_gas, mass_metals, gas_temp_m, gas_teff_m, cell_volume, SFR
    #since our 0-level mesh may have many octs,
    #we must create the octree region sitting 
    #ontop of the first mesh by providing a negative level
	ad = ds.all_data()

	output, grid_structure, nrefined, nleafs = prepare_octree(ds,ile,fle=fle,fre=fre,
													ad=ad,start_level=super_level,
													max_level=max_level, debug=debug)




	output_array = zeros((len(output[output.keys()[0]]), len(output.keys())))
	for i in arange(len(output_array[0])):
		output_array[:,i] = output[output.keys()[i]]
	#grid_structure['level']+=6
	refined = grid_structure['refined']

	np.savez('grid_structure.npz',grid_structure)


	create_fits_file(ds,fn,output,refined,particle_data,fle = ds.domain_left_edge,fre = ds.domain_right_edge)

	return fle, fre, ile, ire, nrefined, nleafs, nstars, output, output_array

def prepare_octree(ds, ile, fle=[0.,0.,0.], fre=[1.,1.,1.], ad=None, start_level=0, max_level=None, debug=True):

	if True: 
		def _MetalMass(field, data):
			return (data['metal_ia_density']*data['cell_volume']).in_units('Msun')
		ad.ds.add_field('MetalMassMsun', function=_MetalMass, units='Msun')
		
		def _TempTimesMass(field, data):
			te = data['thermal_energy']
			hd = data['H_nuclei_density']
			temp = (2.0*te/(3.0*hd*yt.physical_constants.kb)).in_units('K')
			return temp*data["cell_mass"].in_units('Msun')
		ad.ds.add_field('TemperatureTimesCellMassMsun', function=_TempTimesMass, units='K*Msun')
		
		def _cellMassMsun(field, data):
			return data["cell_mass"].in_units('Msun')
		ad.ds.add_field('CellMassMsun', function=_cellMassMsun, units='Msun')

		def _cellVolumeKpc(field, data):
			return data["cell_volume"].in_units('kpc**3')
		ad.ds.add_field('CellVolumeKpc', function=_cellVolumeKpc, units='kpc**3')


        def _pgascgsx(field, data):
            return data['momentum_x'].in_units('Msun/(kpc**2*yr)')*data['cell_volume'].in_units('kpc**3')
        ad.ds.add_field('Cellpgascgsx', function=_pgascgsx, units = 'Msun*kpc/yr')


        def _pgascgsy(field, data):
            return data['momentum_y'].in_units('Msun/(kpc**2*yr)')*data['cell_volume'].in_units('kpc**3')
        ad.ds.add_field('Cellpgascgsy', function=_pgascgsy, units = 'Msun*kpc/yr')

        def _pgascgsz(field, data):    
            return data['momentum_z'].in_units('Msun/(kpc**2*yr)')*data['cell_volume'].in_units('kpc**3')
        ad.ds.add_field('Cellpgascgsz', function=_pgascgsz, units = 'Msun*kpc/yr')




        def _cellSFRtau(field, data):
            min_dens = 0.035 #Msun/pc^3 Ceverino et al. 2009
            density = data["density"].in_units('Msun/pc**3')
            temperature = data["temperature"].in_units('K')
            volume = data["cell_volume"].in_units('pc**3')
            sfr_times_tau = np.where(np.logical_and(density >= min_dens, temperature <= 1.0e4),density*volume,np.zeros_like(density))
            return ds.arr(sfr_times_tau,'Msun')
        ad.ds.add_field('CellSFRtau', function=_cellSFRtau,units='Msun')

        #Tau_SFR = 12 Myr for VELA_v2  Ceverino et al. 2015
        #Not sure about VELA_v2.1 or VELA_v1
        #Using this general version should be applicable for any values used across resolutions
        #Must post-process SFR projections by dividing by Tau.




	fields = ["CellMassMsun","TemperatureTimesCellMassMsun","MetalMassMsun","CellVolumeKpc", "CellSFRtau", "Cellpgascgsx", "Cellpgascgsy", "Cellpgascgsz"]

	#gather the field data from octs 

	print "Retrieving field data"
	field_data = [] 
	for fi,f in enumerate(fields):
		print fi, f
		field_data = ad[f]

	del field_data


	#Initialize dicitionary with arrays containig the needed
	#properites of all octs
	total_octs = ad.index.total_octs

	mask_arr = np.zeros((2,2,2,total_octs), dtype='bool')

	block_iter = ad.blocks.__iter__()  

	for i in np.arange(total_octs):
		oct, mask = block_iter.next()
		mask_arr[:,:,:,i] = mask


	levels = oct._ires[:,:,:, :]
	icoords = oct._icoords[:,:,:, :]
	fcoords = oct._fcoords[:,:,:, :]
	fwidth = oct._fwidth[:,:,:, :]
	mask_arr = mask_arr[:,:,:,:]
	LeftEdge  = (fcoords[0,0,0,:,:]      - fwidth[0,0,0,:,:]*0.5)
	RightEdge = (fcoords[-1,-1,-1,:,:]   + fwidth[-1,-1,-1,:,:]*0.5)


	#RCS commented out fill_octree_arrays, replaced with the code above
	octs_dic = {}
	total_octs = ad.index.total_octs
	octs_dic['LeftEdge'] = LeftEdge[:,:]
	octs_dic['dx']       = fwidth[0,0,0,:,0]
	octs_dic['Level']    = levels[0,0,0,:]

	octs_dic['Fields']    = np.array([ad[f] for f in fields])

	#np.save('octs_dic.npy',octs_dic)

	octs_enclosed = np.argwhere(np.all(octs_dic['LeftEdge'] >= fle, axis=1) & 
								np.all(octs_dic['LeftEdge'] < fre, axis=1))[:,0]
	nocts_enclosed = len(octs_enclosed) 
	print 'Total_octs = %d , Nocts_enclosed = %d' % (total_octs, nocts_enclosed)
	output   = np.zeros((8*nocts_enclosed, len(fields)), dtype='float64')
	output = {}
	for field in fields:
		output[field] = []


	#Location of all octrees, at a given level, and a counter
	oct_loc = {}
	for i in np.arange(ad.index.max_level+1):
		oct_loc[str(i)] = [0,where(levels[0,0,0,:] == i)[0], 0]

	grid_structure 				  = {}
	grid_structure['level'] 	  = []
	grid_structure['refined'] 	  = []
	grid_structure['coords'] 	  = []
	grid_structure['level_index'] = []
	grid_structure['nleafs']      = 0.
	grid_structure['nrefined']    = 0.


	#Initialize, starting at the first high-level octree

	count = [0]
	current_oct_id = 0 
	current_level  = 0
	cell_id = 0
	first_visit = True
	contnue = True
	index_pointer = []
	a = time.time()

	while contnue == True:
		contnue, current_oct_id, current_level, cell_id, first_visit = OctreeDepthFirstHilbert(current_oct_id, current_level, cell_id, 
																									first_visit, oct_loc, mask_arr, 
																											grid_structure, output, octs_dic, fields)
		if oct_loc['0'][0]%5000 == 0: print oct_loc['0'][0]

	b = time.time()
	print b-a

	grid_structure = add_preamble(grid_structure)


	return output, grid_structure, grid_structure['nrefined'], grid_structure['nleafs']

def create_fits_file(ds, fn, output, refined, particle_data, fle, fre):
    #first create the grid structure
    structure = pyfits.Column("structure", format="B", array=array(refined).astype("bool"))
    cols = pyfits.ColDefs([structure])
    st_table = pyfits.BinTableHDU.from_columns(cols)
    st_table.name = "GRIDSTRUCTURE"
    st_table.header.set("hierarch lengthunit", "kpc", comment="Length unit for grid")
    fre = ds.arr(fre, 'code_length').in_units('kpc').value
    fle = ds.arr(fle, 'code_length').in_units('kpc').value
    fdx = fre-fle



    for i,a in enumerate('xyz'):
        st_table.header.set("min%s" % a, fle[i])
        st_table.header.set("max%s" % a, fre[i])
        st_table.header.set("n%s" % a, fdx[i])
        st_table.header.set("subdiv%s" % a, 2)
    st_table.header.set("subdivtp", "OCTREE", "Type of grid subdivision")

    #not the hydro grid data
    fields = ["CellMassMsun","TemperatureTimesCellMassMsun", "MetalMassMsun", "CellVolumeKpc", "CellSFRtau","Cellpgascgsx", "Cellpgascgsy", "Cellpgascgsz"]

    fd = {}
    for i,f in enumerate(fields): 
        fd[f]=array(output[f][:])
    del output

    col_list = []
    cell_mass = fd["CellMassMsun"]
    size = cell_mass.size
    tm = cell_mass.sum()
    col_list.append(pyfits.Column("mass_gas", format='D',
                    array=cell_mass, unit="Msun"))
    col_list.append(pyfits.Column("mass_metals", format='D',
                    array=fd['MetalMassMsun'], unit="Msun"))
    col_list.append(pyfits.Column("gas_temp_m", format='D',
                    array=fd['TemperatureTimesCellMassMsun'], unit="K*Msun"))
    col_list.append(pyfits.Column("gas_teff_m", format='D',
                    array=fd['TemperatureTimesCellMassMsun'], unit="K*Msun"))
    col_list.append(pyfits.Column("cell_volume", format='D',
                    array=fd['CellVolumeKpc'], unit="kpc^3"))
    col_list.append(pyfits.Column("SFR", format='D',
                    array=fd['CellSFRtau'],  unit = 'Msun'))
    p_gas_zipped = zip(fd['Cellpgascgsx'],fd['Cellpgascgsy'],fd['Cellpgascgsz'])
    col_list.append(pyfits.Column("p_gas", format='3D',
                    array=p_gas_zipped , unit = 'Msun*kpc/yr'))


    cols = pyfits.ColDefs(col_list)
    mg_table = pyfits.BinTableHDU.from_columns(cols)
    #mg_table = pyfits.new_table(cols)
    mg_table.header.set("M_g_tot", tm)
    mg_table.header.set("timeunit", "yr")
    mg_table.header.set("tempunit", "K")
    mg_table.name = "GRIDDATA"

    # Add a dummy Primary; might be a better way to do this!
    col_list = [pyfits.Column("dummy", format="E", array=np.zeros(1, dtype='float32'))]
    cols = pyfits.ColDefs(col_list)
    md_table = pyfits.BinTableHDU.from_columns(cols, nrows = len(fd['CellSFRtau']))
    #md_table = pyfits.new_table(cols)
    md_table.header.set("snaptime", ds.current_time.in_units('yr').value[()])
    md_table.name = "YT"

    phdu = pyfits.PrimaryHDU()
    phdu.header.set('nbodycod','yt')
    hls = [phdu, st_table, mg_table,md_table]
    hls.append(particle_data)
    hdus = pyfits.HDUList(hls)
    hdus.writeto(fn, clobber=True)

def round_nocts_wide(dds,fle,fre,nwide=None):
    fc = (fle+fre)/2.0

    assert np.all(fle < fc)
    assert np.all(fre > fc)
    ic = np.rint(fc*dds) #nearest vertex to the center
    ile,ire = ic.astype('int32'),ic.astype('int32')
    cfle,cfre = fc.copy(),fc.copy()
    idx = np.array([0,0,0]) #just a random non-equal array
    width = 0.0
    if nwide is None:
        #expand until borders are included and
        #we have an equaly-sized, non-zero box
        idxq,out=False,True
        while not out or not idxq:
            cfle,cfre = fc-width, fc+width
            #These .ceil and floors were rints (commented by rcs)
            ile = np.floor(cfle*dds).astype('int32')
            ire = np.ceil(cfre*dds).astype('int32')
            idx = ire-ile
            width += 0.1/dds
            #quit if idxq is true:
            idxq = idx[0]>0 and np.all(idx==idx[0])
            out  = np.all(fle>cfle) and np.all(fre<cfre) 
            out &= abs(np.log2(idx[0])-np.rint(np.log2(idx[0])))<1e-5 #nwide should be a power of 2
            assert width[0] < 1.1 #can't go larger than the simulation volume
        nwide = idx[0]
    else:
        #expand until we are nwide cells span
        while not np.all(idx==nwide):
            assert np.any(idx<=nwide)
            cfle,cfre = fc-width, fc+width
            #These .ceil and floors were rints (commented by rcs)
            ile = np.floor(cfle*dds).astype('int32')
            ire = np.ceil(cfre*dds).astype('int32')
            idx = ire-ile
            width += 1e-2*1.0/dds
    assert np.all(idx==nwide)
    assert idx[0]>0
    maxlevel = -np.rint(np.log2(nwide)).astype('int32')
    assert abs(np.log2(nwide)-np.rint(np.log2(nwide)))<1e-5 #nwide should be a power of 2
    return ile,ire,maxlevel,nwide

def prepare_star_particles(ds,star_type,pos=None,vel=None, age=None, creation_time=None,
    initial_mass=None, current_mass=None,metallicity=None, radius = None, 
    fle=[0.,0.,0.],fre=[1.,1.,1.], ad=None):

    if ad is None:
        ad = ds.all_data()
    nump = ad[star_type,"particle_ones"]
    assert nump.sum()>1 #make sure we select more than a single particle
    
    if pos is None:
        pos = yt.YTArray([ad[star_type,"particle_position_%s" % ax]
                        for ax in 'xyz']).transpose()

    idx = np.all(pos > fle, axis=1) & np.all(pos < fre, axis=1)
    assert np.sum(idx)>0 #make sure we select more than a single particle
    pos = pos[idx].in_units('kpc') #unitary units -> kpc

    if creation_time is None:
        formation_time = ad[star_type,"particle_creation_time"][idx].in_units('yr')

    if age is None:
        age = (ds.current_time - formation_time).in_units('yr')

    if vel is None:
        vel = yt.YTArray([ad[star_type,"particle_velocity_%s" % ax]
                        for ax in 'xyz']).transpose()
        # Velocity is cm/s, we want it to be kpc/yr
        #vel *= (ds["kpc"]/ds["cm"]) / (365*24*3600.)
        vel = vel[idx].in_units('kpc/yr')
    
    if initial_mass is None:
        #in solar masses
        initial_mass = ad[star_type,"particle_mass_initial"][idx].in_units('Msun')
    
    if current_mass is None:
        #in solar masses
        current_mass = ad[star_type,"particle_mass"][idx].in_units('Msun')
    
    if metallicity is None:
        #this should be in dimensionless units, metals mass / particle mass
        metallicity = ad[star_type,"particle_metallicity1"][idx]
    
    if radius is None:
        radius = ds.arr(metallicity*0.0 + 10.0/1000.0, 'kpc') #10pc radius
    
    #create every column
    col_list = []
    col_list.append(pyfits.Column("ID", format="J", array=np.arange(current_mass.size).astype('int32')))
    col_list.append(pyfits.Column("parent_ID", format="J", array=np.arange(current_mass.size).astype('int32')))
    col_list.append(pyfits.Column("position", format="3D", array=pos, unit="kpc"))
    col_list.append(pyfits.Column("velocity", format="3D", array=vel, unit="kpc/yr"))
    col_list.append(pyfits.Column("creation_mass", format="D", array=initial_mass, unit="Msun"))
    col_list.append(pyfits.Column("formation_time", format="D", array=formation_time, unit="yr"))
    col_list.append(pyfits.Column("radius", format="D", array=radius, unit="kpc"))
    col_list.append(pyfits.Column("mass", format="D", array=current_mass, unit="Msun"))
    col_list.append(pyfits.Column("age", format="D", array=age,unit='yr'))
    #For particles, Sunrise takes 
    #the dimensionless metallicity, not the mass of the metals
    col_list.append(pyfits.Column("metallicity", format="D",
        array=metallicity,unit="dimensionless")) 
    
    #make the table
    cols = pyfits.ColDefs(col_list)
    pd_table = pyfits.BinTableHDU.from_columns(cols)
    #pd_table = pyfits.new_table(cols)
    pd_table.name = "PARTICLEDATA"
    
    #make sure we have nonzero particle number
    assert pd_table.data.shape[0]>0
    return pd_table, np.sum(idx)

def save_to_gridstructure(grid_structure, level, refined, leaf):
	'''
	Function to save grid information
	'''
	grid_structure['level'].append(level)
	grid_structure['refined'].append(refined)
	if leaf:
		grid_structure['nleafs']+=1
	if refined:
		grid_structure['nrefined']+=1
	return

def add_preamble(grid_structure):
	'''
	Add the higher level oct information
	'''

	level_list = blist(grid_structure['level'])
	level_arr = array(grid_structure['level'])
	list_structure = blist(grid_structure['refined'])

	print 'Adding preamble oct structure...'
	i = 0
	while True:
		good = where(level_arr == i)[0]
		print i, len(good)
		if len(good) == 1:
			break
		else:
			i -= 1
		good_red = good[0:len(good):8]

		good = where(level_arr == i)[0]

		for n, back_good in enumerate(good_red[::-1]):
			level_list.insert(back_good, i)
			list_structure.insert(back_good, True)
			grid_structure['nrefined']+=1


		level_arr = array(level_list)



	grid_structure['refined'] = list_structure
	grid_structure['level'] = level_list
	return grid_structure

def OctreeDepthFirstHilbert(current_oct_id, current_level, cell_id, first_visit, oct_loc, mask_arr, grid_structure,  output, octs_dic, field_names, debug = False):
	
	#This is the oct mask, a 2x2x2 array that will identify leafs (True) and refined cells (False)
	mask_i = mask_arr[:,:,:, current_oct_id]
	#Flatten the mask to 8x1
	flat_mask = mask_i.ravel(order = 'F')
	fields = octs_dic['Fields'][:,:,:,:,current_oct_id]
	#p_gas_field = octs_dic['p_gas_field']
	fields_all = zeros((fields.shape[0], 8))
	for field_index in range(fields.shape[0]):
		fields_all[field_index] = fields[field_index,:,:,:].ravel(order = 'F')



	if first_visit == True:
		#It's the first time visiting this oct, so let's save 
		#the oct information here in our grid structure dict
		save_to_gridstructure(grid_structure, current_level, refined = True, leaf = False)
		if debug: print '\t'*current_level, current_level
		
	#If this is our first visit, then cell_id should be 0. If we are returning to a parent 
	#from a child, then cell_id will be larger than 0 (i.e. we don't want to return to the previous cells we visited)
	for i in arange(cell_id, len(flat_mask)): 
		is_leaf = flat_mask[i]

		if is_leaf:
			#This cell is a leaf, save the grid information and the physical properties
			save_to_gridstructure(grid_structure, current_level+1, refined = False, leaf = True)
			
			for field_index in range(fields.shape[0]):
				output[field_names[field_index]].append(fields_all[field_index,i])


			if debug:  print '\t'*current_level,current_level, 'Found a leaf in cell '+ str(int(i))+'/'+str(int(len(mask_i.ravel())))
		else:
			#This cell is not a leaf, we'll now advance in to this cell
			oct_loc[str(current_level)][2] = i + 1    #Store the cell_id that we are currently in (+1, since we'll be returning to the next cell)
			child_level 	= current_level + 1 	  #We'll be going one level down
			#Identify the index of the next available octree at child_level
			child_argument 	= oct_loc[str(child_level)][0] 
			child_oct_id 	= oct_loc[str(child_level)][1][child_argument]   

			if debug: print '\t'*current_level, current_level, 'Found a cell to refine in cell ' + str(int(i))+'/'+str(int(len(flat_mask)))

			#Visiting the next cell, going back to the for-loop first
			return True, child_oct_id, child_level, 0., True


	#Done with that oct, now we need to return to the parent level (or a new high level oct)
	#But first, let's set the cell_id quantity for this level back to 0 (starting fresh)

	oct_loc[str(current_level)][0] 	+= 1  	    #Next time we draw from this level, we will select the next index
	oct_loc[str(current_level)][2] 	= 0
	

	#ok let's see what our path is now
	if current_level > 0:
		#We're going to move up to a parent
		parent_level	= current_level - 1
		#The octree index for our current parent
		parent_argument = oct_loc[str(parent_level)][0]
		parent_oct_id 	= oct_loc[str(parent_level)][1][parent_argument]
		parent_cell_id 	= oct_loc[str(parent_level)][2]
		#Don't reset the counter just yet, need to check if we are at cell 8 in the parent
		first_visit = False       #If it turns out we need to return to parent
		while parent_cell_id == 8:
			#Done with this parent oct, reset the cells at this level and try the next level
			oct_loc[str(parent_level)][2] 	= 0
			oct_loc[str(parent_level)][0] 	+= 1
			parent_level -= 1
			if parent_level == -1: #Did we reach the high-level octs?
				parent_level = 0
				first_visit = True    #Done with this parent, on to a new high-level oct
			#Identify next index and check if we've reached the end of the high-level octs
			parent_argument = oct_loc[str(parent_level)][0]
			if (parent_level == 0) & (parent_argument == len(oct_loc[str(parent_level)][1])): 
				return False, -1, -1, -1, -1    #Done!
			else:
				#Nope, still more octs...
				parent_oct_id 	= oct_loc[str(parent_level)][1][parent_argument]
				parent_cell_id 	= oct_loc[str(parent_level)][2]


		return True, parent_oct_id, parent_level, parent_cell_id, first_visit



	elif current_level == 0:
		parent_level = current_level    #Moving to next high-level oct...
		first_visit = True       		#...which will be our first visit
		#Find next oct
		parent_argument = oct_loc[str(parent_level)][0]
		#Check if we've reached the end of the high-level oct
		if (parent_argument == len(oct_loc[str(parent_level)][1])): 
			return False, -1, -1, -1, -1   #Done!
		parent_oct_id 	= oct_loc[str(parent_level)][1][parent_argument]
		parent_cell_id 	= oct_loc[str(parent_level)][2]
		return True, parent_oct_id, parent_level, parent_cell_id, first_visit





































