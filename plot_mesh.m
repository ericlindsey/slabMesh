function plot_mesh(fname)
    % plot blocks mesh file: assumed .mat file format, contains variables:
    % nc: int, number of vertices
    % nEl: int, number of triangles
    % c: (nc x 3) double, contains lon,lat,depth(negative) for each vertex
    % v: (nEl x 3) int, contains triples of vertex IDs for each triangle
    %
    % Eric Lindsey, Aug 2020
    
    % load the data
    load(fname,'nc','nEl','c','v');
    x=c(:,1);
    y=c(:,2);
    z=c(:,3);
    
    % get centroid depths for each triangle
    depcol=zeros(nEl,1);
    for i=1:nEl
        depths=z(v(i,:));
        depcol(i)=sum(depths)/3;
    end
    % plot triangles
    trisurf(v,x,y,z,depcol)

end
