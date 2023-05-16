function plot_mesh(fname,varargin)
    % plot blocks mesh file: assumed .mat file format, contains variables:
    % nc: int, number of vertices
    % nEl: int, number of triangles
    % c: (nc x 3) double, contains lon,lat,depth(negative) for each vertex
    % v: (nEl x 3) int, contains triples of vertex IDs for each triangle
    %
    % optional: vector containing values to use for triangle colors
    % (default: fault depth)
    %
    % Eric Lindsey, Aug 2020
    
    % load the data
    load(fname,'nc','nEl','c','v');
    x=c(:,1);
    y=c(:,2);
    z=c(:,3);
    
    if nargin == 1
        % use centroid depths for triangle colors
        colorvec=zeros(nEl,1);
        for i=1:nEl
            depths=z(v(i,:));
            colorvec(i)=sum(depths)/3;
        end
    else
        colorvec = varargin{1};
    end
        
    % plot triangles
    trisurf(v,x,y,z,colorvec)

end
