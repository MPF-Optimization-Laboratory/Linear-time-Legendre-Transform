function [ ] = plot3d( fpath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

output = load(fpath);
x = output(:,1);
y = output(:,2);
z = output(:,3);

tri = delaunay(x,y);
trisurf(tri, x,y,z);

end

