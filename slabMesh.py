# set of functions for reading a slab2.0 grid file and converting to a triangular mesh

#default python packages
from __future__ import division
from __future__ import absolute_import
from six.moves import range
# required anaconda packages (conda install): numpy, scipy, netcdf4, matplotlib
import numpy as np
import numpy.linalg as la
from scipy import interpolate
import scipy.io as scio
from netCDF4 import Dataset
# available from conda forge (conda install -c conda-forge meshpy)
import meshpy.triangle as triangle

# slab 2.0 import and manipulation functions
######

def load_slab2_grid(fname):
    # load slab2.0 grid file and return X,Y,Z points as 2D matrices.
    # note that the slab2.0 format has lots of NaN values where the slab is not defined
    with Dataset(fname, mode='r') as fh:
        lons = fh.variables['x'][:]
        lats = fh.variables['y'][:]
        Z    = fh.variables['z'][:]
    [X,Y]=np.meshgrid(lons,lats)
    return X,Y,Z

def crop_rectangle(X,Y,Z,lonmin=0,lonmax=360,latmin=-90,latmax=90):
    # crop a rectangle out of a 3D point dataset
    # crop a section of the fault
    # we need lats,lons which are column vectors.
    lons=X[0,:]
    lats=Y[:,0]
    # need a zero because this outputs a tuple...
    rowkeep=np.where( (lats>latmin) & (lats<=latmax) )[0] 
    colkeep=np.where( (lons>lonmin) & (lons<=lonmax) )[0]
    Ikeep=np.ix_(rowkeep,colkeep) # virtually unknown function to do a very simple task
    Xcrop=X[Ikeep]
    Ycrop=Y[Ikeep]
    Zcrop=Z[Ikeep]
    return Xcrop,Ycrop,Zcrop

def slabtop_at_zero(Z):
    # offset Z values so that 0 is maximum value
    Zoffset=Z-np.max(Z)
    return Zoffset

# interpolation that can handle nans
######
def my_griddata(xin,yin,zin,xout,yout):
    ''' fix a problem with griddata returning nans - use nearest neighbor interp to fill'''
    
    # construct (n,2) points array for griddata - remove nans first
    I_no_innan=np.where(~np.isnan(zin))
    points=np.column_stack((xin[I_no_innan].ravel(), yin[I_no_innan].ravel()))
    print(np.shape(points))
    # first use of griddata
    zout = interpolate.griddata(points, zin[I_no_innan].ravel(), (xout,yout))
    
    #now, find any output nans and use method='nearest' for them (generally points outside the convex hull)
    I_outnan=np.where(np.isnan(zout))
    zout[I_outnan] = interpolate.griddata(points,zin[I_no_innan].ravel(),
                                          (xout[I_outnan],yout[I_outnan]), method='nearest')
    return zout

# meshing functions
######
def round_trip_connect(start, end):
    return [(i, i+1) for i in range(start, end)] + [(end, start)]

def needs_refinement(vertices, area):
    bary = np.sum(np.array(vertices), axis=0)/3
    max_area = 0.005 + (la.norm(bary, np.inf)-1)*0.0005
    return bool(area > max_area)

def make_mesh(Xpoly,Ypoly):
    # create the mesh using meshpy.triangle
    points=[pt for pt in zip(Xpoly,Ypoly)]
    facets = round_trip_connect(0, len(points)-1)
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets)
    mesh = triangle.build(info, refinement_func=needs_refinement)
    mesh_points = np.array(mesh.points)
    mesh_tris = np.array(mesh.elements)
    return mesh_points,mesh_tris

# output functions
######
def save_mesh_for_blocks(pts,tris,fname):
    nc=len(pts)
    nEl=len(tris)
    #set indexing relative to 1 for matlab
    tris_out = tris + 1
    outdict={'c':pts,'nc':nc,'nEl':nEl,'v':tris_out}
    scio.savemat(fname,outdict)
    
def save_mesh_for_unicycle(pts,tris,fname):
    #nc=len(pts)
    #nEl=len(tris)
    #set indexing relative to 1 for matlab
    tris_out = tris + 1
    nedfname=fname+'.ned'
    trifname=fname+'.tri'
    with open(nedfname, 'w') as f:
        for index, point in enumerate(pts):
            #asterisk unpacks the points list into 3 items
            f.write("{} {} {} {}\n".format(index + 1, *point))
    with open(trifname, 'w') as f:
        for index,tri in enumerate(tris_out):
            f.write("{} {} {} {} 90.0\n".format(index + 1, *tri))