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

# Define a function that returns a refinement function with the threshold captured:
def get_simple_refinement(threshold, max_refinements):

    # Use a mutable container to track the count
    refinement_count = [0]

    def refinement(vertices, area):
        # Increase the counter for each triangle evaluated.
        refinement_count[0] += 1
        if refinement_count[0] > max_refinements:
            # If we have reached the maximum, stop refining.
            print("Maximum refinements reached.")
            return False
        return area > threshold
    
    return refinement

def get_depth_based_refinement(depth_interp, base_threshold=0.03, factor=0.00007, max_refinements=1e5):
    """
    Returns a refinement function that uses depth information.
    
    Parameters:
      depth_interp: callable, a function f(x, y) that returns depth z.
      base_threshold: float, the baseline area threshold.
      factor: float, scaling factor to adjust threshold based on depth.
      max_refinements: int, maximum number of refinement calls allowed.

    The threshold is computed as:
      threshold = base_threshold + factor * (depth - 1)
    so deeper (larger depth) triangles get a higher threshold;
    this assumes positive depth values.

    """
    # Use a mutable container to track the count
    refinement_count = [0]

    def refinement(vertices, area):
        # Increase the counter for each triangle evaluated.
        refinement_count[0] += 1
        if refinement_count[0] > max_refinements:
            # If we have reached the maximum, stop refining.
            print("Maximum refinements reached.")
            return False
        # Compute the barycenter (mean of the vertices).
        bary = np.mean(vertices, axis=0)
        # Lookup the depth at the barycenter. Handle missing values.
        depth = depth_interp(bary[0], bary[1])
        if depth is None or np.isnan(depth):
            depth = 0  # Or choose a default, safe value.
        threshold = base_threshold + factor * (depth - 1)
        # Return True if the triangle area exceeds the threshold.
        return area > threshold
    
    return refinement


def fix_vertical_edges(X, Y, tolerance=1e-8, offset=1e-6):
    """
    Adjusts vertices of a polygon to avoid vertical (or nearly vertical) segments.
    
    Parameters:
      X, Y : array-like
          Coordinates of the polygon vertices.
      tolerance : float
          Tolerance for considering a segment as vertical (difference in X).
      offset : float
          The fixed offset added to the x-coordinate to break verticality.
    
    Returns:
      X_fixed : np.ndarray
          The adjusted x-coordinates.
      Y : np.ndarray
          The y-coordinates (unchanged).
    """
    X_fixed = np.array(X, dtype=float)
    n = len(X_fixed)
    for i in range(n):
        # Wrap-around index for closed polygon
        j = (i + 1) % n
        dx = X_fixed[j] - X_fixed[i]
        dy = Y[j] - Y[i]
        # Calculate the angle in radians using arctan2; vertical if angle near pi/2 or -pi/2.
        angle = np.arctan2(dy, dx)
        if np.abs(np.abs(angle) - np.pi/2) < 1e-4 or np.abs(dx) < tolerance:
            # Deterministically adjust the second vertex's x value.
            # Here we always add a small offset. Alternatively, you could decide based on the vertex's position.
            X_fixed[j] = X_fixed[i] + offset
            print(f"Adjusted vertex {j} from x={X[j]} to x={X_fixed[j]} to avoid vertical segment.")
    return X_fixed, np.array(Y)

def make_mesh(Xpoly,Ypoly,threshold=0.03,max_refinements=100000):
    # create the mesh using meshpy.triangle
    try:
        Xpoly_fixed, Ypoly_fixed = fix_vertical_edges(Xpoly, Ypoly)
        points=[pt for pt in zip(Xpoly_fixed,Ypoly_fixed)]
        facets = round_trip_connect(0, len(points)-1)
        info = triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(facets)
        refinement_func = get_simple_refinement(threshold, max_refinements)
        mesh = triangle.build(info, refinement_func=refinement_func)
        mesh_points = np.array(mesh.points)
        mesh_tris = np.array(mesh.elements)
    except Exception as e:
        print("error:", e)
    return mesh_points,mesh_tris


def make_depth_variable_mesh(Xpoly,Ypoly,depth_interp,threshold=0.03,factor=0.00007,max_refinements=100000):
    # create the mesh using meshpy.triangle
    try:
        Xpoly_fixed, Ypoly_fixed = fix_vertical_edges(Xpoly, Ypoly)
        points=[pt for pt in zip(Xpoly_fixed,Ypoly_fixed)]
        facets = round_trip_connect(0, len(points)-1)
        info = triangle.MeshInfo()
        info.set_points(points)
        info.set_facets(facets)
        
        # Create the depth-based refinement function.
        refinement_func = get_depth_based_refinement(depth_interp, threshold, factor, max_refinements)

        mesh = triangle.build(info, refinement_func=refinement_func)
        mesh_points = np.array(mesh.points)
        mesh_tris = np.array(mesh.elements)
    except Exception as e:
        print("error:", e)
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