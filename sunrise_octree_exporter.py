"""
Code to export from yt to Sunrise


"""


import time
import numpy as np
from numpy import *
import pyfits
import yt
from yt.funcs import get_pbar
from blist import blist
from yt.funcs import *
#import sunrise_export_gfs.octree_to_depthFirstHilbert_GFS as depthFirstHilbert



class hilbert_state():
	def __init__(self,dim=None,sgn=None,octant=None):
		if dim is None: dim = [0,1,2]
		if sgn is None: sgn = [1,1,1]
		if octant is None: octant = 5
		self.dim = dim
		self.sgn = sgn
		self.octant = octant

	def flip(self,i):
	    self.sgn[i]*=-1
	
	def swap(self,i,j):
	    temp = self.dim[i]
	    self.dim[i]=self.dim[j]
	    self.dim[j]=temp
	    axis = self.sgn[i]
	    self.sgn[i] = self.sgn[j]
	    self.sgn[j] = axis
	
	def reorder(self,i,j,k):
	    ndim = [self.dim[i],self.dim[j],self.dim[k]] 
	    nsgn = [self.sgn[i],self.sgn[j],self.sgn[k]]
	    self.dim = ndim
	    self.sgn = nsgn

	def copy(self):
	    return hilbert_state([self.dim[0],self.dim[1],self.dim[2]],
	                         [self.sgn[0],self.sgn[1],self.sgn[2]],
	                         self.octant)

	def descend(self,o):
		child = self.copy()
		child.octant = o
		if o==0:
		    child.swap(0,2)
		elif o==1:
		    child.swap(1,2)
		elif o==2:
		    pass
		elif o==3:
		    child.flip(0)
		    child.flip(2)
		    child.reorder(2,0,1)
		elif o==4:
		    child.flip(0)
		    child.flip(1)
		    child.reorder(2,0,1)
		elif o==5:
		    pass
		elif o==6:
			child.flip(1)
			child.flip(2)
			child.swap(1,2)
		elif o==7:
			child.flip(0)
			child.flip(2)
			child.swap(0,2)
		return child



	def __iter__(self):
		vertex = np.array([0,0,0]).astype('int32')
		j = 0
		for i in range(3):
		    vertex[self.dim[i]] = 0 if self.sgn[i]>0 else 1
		yield vertex, self.descend(j)
		vertex[self.dim[0]] += self.sgn[0]
		j+=1
		yield vertex, self.descend(j)
		vertex[self.dim[1]] += self.sgn[1] 
		j+=1
		yield vertex, self.descend(j)
		vertex[self.dim[0]] -= self.sgn[0] 
		j+=1
		yield vertex, self.descend(j)
		vertex[self.dim[2]] += self.sgn[2] 
		j+=1
		yield vertex, self.descend(j)
		vertex[self.dim[0]] += self.sgn[0] 
		j+=1
		yield vertex, self.descend(j)
		vertex[self.dim[1]] -= self.sgn[1] 
		j+=1
		yield vertex, self.descend(j)
		vertex[self.dim[0]] -= self.sgn[0] 
		j+=1
		yield vertex, self.descend(j)



class oct_object():
	def __init__(self, is_leaf, fcoords, fwidth, level, oct_id, child_oct_ids):
		self.is_leaf = is_leaf
		self.fcoords = fcoords
		self.fwidth = fwidth
		self.le = fcoords - 0.5*fwidth
		self.re = fcoords + 0.5*fwidth
		self.child_oct_ids = child_oct_ids
		self.n_refined_visited = 0
		self.level = level
		self.child_level = self.level + 1
		self.oct_id = oct_id		


def OctreeDepthFirstHilbert(current_oct_id, current_level, mask_arr, hilbert, \
							fcoords, fwidth, grid_structure,\
							output, octs_dic, oct_loc, field_names, debug = False, f  = 'out.out', preamble_end = None):


	if current_oct_id%10000 == 0 : print str(current_oct_id) + '/' + str(shape(mask_arr)[3])
	save_to_gridstructure(grid_structure, current_level, refined = True, leaf = False)
			

	mask_i = mask_arr[:,:,:, current_oct_id]
	fcoords_ix, fcoords_iy, fcoords_iz = fcoords[:,:,:, current_oct_id, 0],  fcoords[:,:,:, current_oct_id, 1], fcoords[:,:,:, current_oct_id, 2]
	fwidth_ix, fwidth_iy, fwidth_iz = fwidth[:,:,:, current_oct_id, 0],  fwidth[:,:,:, current_oct_id, 1], fwidth[:,:,:, current_oct_id, 2]

	flat_mask = mask_i.ravel(order = 'F')
	flat_fcoords = array(zip(fcoords_ix.ravel(order = 'F').value[()], fcoords_iy.ravel(order = 'F').value[()], fcoords_iz.ravel(order = 'F').value[()]))
	flat_fwidth = array(zip(fwidth_ix.ravel(order = 'F').value[()], fwidth_iy.ravel(order = 'F').value[()], fwidth_iz.ravel(order = 'F').value[()]))
	refined_locations = where(flat_mask == False)[0]


	fields = octs_dic['Fields'][:,:,:,:, current_oct_id-preamble_end]
	fields_all = zeros((fields.shape[0], 8))
	for field_index in range(fields.shape[0]):
		fields_all[field_index] = fields[field_index,:,:,:].ravel(order = 'F')

	#It's the first time visiting this oct, so let's save 
	#the oct information here in our grid structure dictionary
	if debug: f.write('\t'*current_level+'Entering level %i oct: found %i refined cells and %i leaf cells\n'%(current_level, len(refined_locations), 8-len(refined_locations)))

	child_level	= current_level	+ 1

	if len(refined_locations) > 0:
		child_oct_ids = oct_loc[str(child_level)][1][oct_loc[str(child_level)][0]:oct_loc[str(child_level)][0]+len(refined_locations)]
		oct_loc[str(child_level)][0] += len(refined_locations)
	else:
		child_oct_ids = None

	oct_obj = oct_object(flat_mask, flat_fcoords, flat_fwidth, current_level, current_oct_id, child_oct_ids)

	for (vertex, hilbert_child) in hilbert:
		parent_oct_le = oct_obj.le[0]
		vertex_new = vertex*oct_obj.fwidth[0]
		next_child_le = parent_oct_le + vertex_new
		i = where((oct_obj.le[:,0] == next_child_le[0]) & (oct_obj.le[:,1] == next_child_le[1]) & (oct_obj.le[:,2] == next_child_le[2]))[0][0]

		if oct_obj.is_leaf[i]:		
				if debug:  f.write('\t'*oct_obj.child_level+str(oct_obj.child_level) + '\tFound a leaf in cell %i/%i \t (x,y,z) = (%.8f, %.8f, %.8f) \n'%(i, 8, oct_obj.le[i][0], oct_obj.le[i][1], oct_obj.le[i][2]))
				#assert(current_oct_id > preamble_end)
				#This cell is a leaf, save the grid information and the physical properties
				save_to_gridstructure(grid_structure, current_level, refined = False, leaf = True)				
				assert(current_oct_id >= preamble_end)
				for field_index in range(fields.shape[0]):
					output[field_names[field_index]].append(fields_all[field_index,i])


		else:
			#This cell is not a leaf, we'll now advance in to this cell
			if debug:  f.write('\t'*child_level+str(child_level) + '\tFound a refinement in cell %i/%i \t (x,y,z) = (%.8f, %.8f, %.8f) \n'%(i, 8, oct_obj.le[i][0], oct_obj.le[i][1], oct_obj.le[i][2]))
			OctreeDepthFirstHilbert(oct_obj.child_oct_ids[oct_obj.n_refined_visited], oct_obj.child_level, mask_arr, 
									hilbert_child, fcoords, fwidth, grid_structure, output, octs_dic, oct_loc, 
									field_names, debug, f, preamble_end)

			oct_obj.n_refined_visited += 1



def add_preamble(levels, fwidth, fcoords, LeftEdge, RightEdge, mask_arr):

	i = 0
	while True:
		good = where(levels[0,0,0,:] == i)[0]
		print i, len(good)
		if len(good) == 1:
			return levels, fwidth, fcoords, mask_arr
		else:
			i -= 1


		temp_levels  = zeros((2,2,2,len(good)/8))
		temp_fwidth  = zeros((2,2,2,len(good)/8, 3))
		temp_fcoords = zeros((2,2,2,len(good)/8, 3))
		temp_maskarr = np.zeros((2,2,2,len(good)/8)).astype('bool')


		for p, gd in enumerate(good):
			print p, '\t', LeftEdge[gd]
			raw_input('')


		good_red = good[0:len(good):8]


		for n, back_good in enumerate(good_red):
			temp_levels[:,:,:,n] += i
			temp_fwidth[:,:,:,n,:] += fwidth[:,:,:,back_good,:]*2.



			child_cen = fcoords[0,0,0,back_good,:].value[()]+fwidth[0,0,0,back_good,:].value[()]			

			for ii in arange(2):
				for jj in arange(2):
					for kk in arange(2):
						temp_fcoords[ii,jj,kk, n,:] = child_cen + 0.5*temp_fwidth[ii,jj,kk,n,:]*[[(-1)**(1-ii), (-1)**(1-jj), (-1)**(1-kk)]]


			#ce = (temp_fcoords[0,0,0,n,:] + temp_fwidth[0,0,0,n,:]*0.5)
			#print n, '\t', ce
			#print raw_input('')



		levels = append(temp_levels, levels, axis = 3)
		fwidth = append(temp_fwidth, fwidth, axis = 3)
		fcoords = append(temp_fcoords, fcoords, axis = 3)
		mask_arr = append(temp_maskarr, mask_arr, axis = 3)


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

	if False: 
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




        if False:
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



	if False:
		output = {}
		for field in fields:
		    output[field] = []

	aa = shape(levels)[3]
	levels, fwidth, fcoords,  mask_arr = add_preamble(levels, fwidth, fcoords, LeftEdge, RightEdge, mask_arr)
	bb = shape(levels)[3]
	preamble_end = bb - aa
	print bb, aa

	levels = levels - min(levels[0,0,0,:])
	levels = levels.astype('int')
	

	#RCS commented out fill_octree_arrays, replaced with the code above
	octs_dic = {}
	total_octs = ad.index.total_octs
	octs_dic['LeftEdge'] = LeftEdge[:,:]
	octs_dic['dx']       = fwidth[0,0,0,:,0]
	octs_dic['Level']    = levels[0,0,0,:]

	octs_dic['Fields']    = np.array([ad[f] for f in fields])



	grid_structure 				  = {}
	grid_structure['level'] 	  = []
	grid_structure['refined'] 	  = []
	grid_structure['coords'] 	  = []
	grid_structure['level_index'] = []
	grid_structure['nleafs']      = 0.
	grid_structure['nrefined']    = 0.

	hs = hilbert_state()
	#Location of all octrees, at a given level, and a counter
	oct_loc = {}
	for i in np.arange(max(levels[0,0,0,:])+1):
		oct_loc[str(i)] = [0,where(levels[0,0,0,:] == i)[0]]

	outfile = open('debug_hilbert.out', 'w+')
	a = time.time()
	debug = True
	OctreeDepthFirstHilbert(current_oct_id = 0, current_level = 0, mask_arr = mask_arr,  hilbert = hs, 
							fcoords = fcoords, fwidth = fwidth, grid_structure = grid_structure, 
							output = output, octs_dic = octs_dic, oct_loc = oct_loc, field_names = fields, 
							debug = debug, f  = outfile, preamble_end = preamble_end)
	b = time.time()


	print b-a


	outfile.close()




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













































