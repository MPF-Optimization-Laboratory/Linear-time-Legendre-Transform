function [ ] = plot2d( fpath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

output = load(fpath);
x = output(:,1);
y = output(:,2);

plot(x,y);

end

