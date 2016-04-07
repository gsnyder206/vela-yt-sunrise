import yt
import time



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

def OctreeDepthFirstHilbert(current_oct_id, current_level, mask_arr, hilbert, fcoords, fwidth, debug = False, f  = 'out.out'):

	mask_i = mask_arr[:,:,:, current_oct_id]
	fcoords_ix, fcoords_iy, fcoords_iz = fcoords[:,:,:, current_oct_id, 0],  fcoords[:,:,:, current_oct_id, 1], fcoords[:,:,:, current_oct_id, 2]
	fwidth_ix, fwidth_iy, fwidth_iz = fwidth[:,:,:, current_oct_id, 0],  fwidth[:,:,:, current_oct_id, 1], fwidth[:,:,:, current_oct_id, 2]

	flat_mask = mask_i.ravel(order = 'F')
	flat_fcoords = array(zip(fcoords_ix.ravel(order = 'F').value[()], fcoords_iy.ravel(order = 'F').value[()], fcoords_iz.ravel(order = 'F').value[()]))
	flat_fwidth = array(zip(fwidth_ix.ravel(order = 'F').value[()], fwidth_iy.ravel(order = 'F').value[()], fwidth_iz.ravel(order = 'F').value[()]))
	refined_locations = where(flat_mask == False)[0]

	#It's the first time visiting this oct, so let's save 
	#the oct information here in our grid structure dictionary
	if debug: f.write('\t'*current_level+'Entering level %i oct: found %i refined cells and %i leaf cells\n'%(current_level, len(refined_locations), 8-len(refined_locations)))

	child_level	= current_level	+ 1

	hilbert_order = arange(8)
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
				#This cell is a leaf, save the grid information and the physical properties
				if debug:  f.write('\t'*oct_obj.child_level+str(oct_obj.child_level) + '\tFound a leaf in cell %i/%i \t (x,y,z) = (%.8f, %.8f, %.8f) \n'%(i, 8, oct_obj.fcoords[i][0], oct_obj.fcoords[i][1], oct_obj.fcoords[i][2]))
		else:
			#This cell is not a leaf, we'll now advance in to this cell
			if debug:  f.write('\t'*child_level+str(child_level) + '\tFound a refinement in cell %i/%i \t (x,y,z) = (%.8f, %.8f, %.8f) \n'%(i, 8, oct_obj.fcoords[i][0], oct_obj.fcoords[i][1], oct_obj.fcoords[i][2]))
			OctreeDepthFirstHilbert(oct_obj.child_oct_ids[oct_obj.n_refined_visited], oct_obj.child_level, mask_arr, hilbert_child, fcoords, fwidth, debug = debug, f  = outfile)
			oct_obj.n_refined_visited += 1







if __name__ == '__main__':
	gen_name, gal_name, snap_name, snaps  = 'VELA_v2', 'VELA27', 'VELA27_a0.370', '../data/VELA27_v2/a0.370/10MpcBox_csf512_a0.370.d'
	ds = yt.load(snaps, limit_level = 2)
	snap_dir = '/Volumes/wd/yt_pipeline/Runs/%s/%s/%s'%(gen_name, gal_name, snap_name+'_sunrise')



	hs = hilbert_state()


	if True:
		ad = ds.all_data()

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




	if True:
		#Location of all octrees, at a given level, and a counter
		oct_loc = {}
		for i in np.arange(ad.index.max_level+1):
			oct_loc[str(i)] = [0,where(levels[0,0,0,:] == i)[0]]

		outfile = open('debug_hilbert.out', 'w+')


		a = time.time()
		for i in arange(len(oct_loc['0'][1])):

			if (oct_loc['0'][0] < 120000) & (oct_loc['0'][0] < 125000):	
				debug = True
			else:
				debug = False

			current_oct_id = oct_loc['0'][1][i]


			OctreeDepthFirstHilbert(current_oct_id = current_oct_id, current_level = 0, mask_arr = mask_arr,  hilbert = hs, fcoords = fcoords, fwidth = fwidth, debug = debug, f  = outfile)
			oct_loc['0'][0] += 1

			if i%10000 == 0: print str(i)+'/'+str(len(oct_loc['0'][1]))

		b = time.time()


		print 'Final time in seconds: ', b - a
		outfile.close()

























