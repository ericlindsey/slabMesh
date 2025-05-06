slabMesh: Triangular meshing for megathrust slabs and other faults
-----
Eric Lindsey

Last updated May 2025 - include Sagaing fault mesh

This code contains a very simplified set of routines for creating a triangular fault mesh for use with your favorite slip-inversion program. Currently there are a number of manual steps and many hardcoded parameters - this is research-grade code, not commercial! Use at your own risk.

Steps to make a mesh, as demonstrated in the included jupyter notebook:

1. Read a Slab2.0 NetCDF grid file and crop to a region of interest, or a fault trace .csv file
2. Extend the slab up to a surface trace, if needed. 
3. Determine the polygon outline of the region you wish to mesh - from the surface down to a specified depth.
4. Create the mesh, working in 2D (one dimension is ignored)
5. Extrude the mesh to include the correct depth (or, for the Sagaing fault, longitude) at each point. This will stretch out your triangles in the down-dip direction, but it is a small effect and the alternative (fully 3D meshing) is much too complicated.
6. Output triangles and vertices in desired format.

The code also includes a modified version of Python's `griddata` method that avoids the introduction of NaN values. 
