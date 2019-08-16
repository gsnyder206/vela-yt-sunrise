import astropy.io.fits as fits
import numpy as np



def update_star_radius(sf):


    sfo = fits.open(sf)

    pd=sfo['PARTICLEDATA']

    pos=pd.data['position']  #(N,3) array

    N=pos.shape[0] #number of star particles

    #test values only!  future: 
    #actually, need to verify that we get ALL particles here

    mins=np.min(pos,axis=0)
    maxs=np.max(pos,axis=0)
    
    limits=[(mins[0]-1.0,maxs[0]+1.0),(mins[1]-1.0,maxs[1]+1.0),(mins[2]-1.0,maxs[2]+1.0)]

    alm = np.argmax([limits[0][1]-limits[0][0], limits[1][1]-limits[1][0], limits[2][1]-limits[2][0]])

    target_size= 0.2 #in kpc
    
    bins=np.int32((limits[alm][1]-limits[alm][0])/target_size)

    binsize=(limits[alm][1]-limits[alm][0])/bins
    print('star position limits: ', limits)
    print('star radius histogram bins: ', bins)
    print('star radius histogram binsize: ', binsize)
    
    H,edges = np.histogramdd(pos,bins=bins,range=limits)

    xi=np.digitize(pos[:,0],edges[0],right=True)-1 #subtract 1 because 0 is reserved for below left bound -- can never happen
    yi=np.digitize(pos[:,1],edges[1],right=True)-1
    zi=np.digitize(pos[:,2],edges[2],right=True)-1

    assert(np.min(xi) >= 0)
    assert(np.min(yi) >= 0)
    assert(np.min(zi) >= 0)

    Nvals=H[xi,yi,zi]

    Q = binsize*(Nvals/64.0)**(-1.0/3.0)
    
    uval=np.min([Q,Q*0.0+10.0],axis=0)

    finalval=np.max([0.01+0.0*uval,uval],axis=0)

    return H, edges, Nvals, Q, uval, finalval
