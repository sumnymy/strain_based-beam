function [h_new,pos_new] = configRecovery(strain,h_old,dl)
%--------------------------------------------------------------------------
% Purpose:
%         To recover the nodal vector from strain data
% Synopsis :
%         [h,pos] = configRecovery(strain)
% Variable Description:
% Input:
%           strain - strains of the beam element
%           -----> strain = [e1;e2;...;en]   dimension: nelem*4
%           h_old  - nodal vector at ith step
%           -----> h_old = [h1;h2;...]   dimension:nnode*12
%           pos_old  - position vector at ith step
%           -----> pos_old = [p1;p2;...]   dimension:nnode*3
%           dl  -  length of each element 
% Output:
%           h_new  - nodal vector at (i+1)th step
%           -----> h_new = [h1;h2;...]   dimension:nnode*12
%           pos_new  - position vector at (i+1)th step
%           -----> pos_new = [p1;p2;...]   dimension:nnode*3           
%--------------------------------------------------------------------------
% obtain number of elements and nodes
nelem = size(strain,1);
nnode = size(h_old,1);
%for the first element:i = 1
expG1 = elementKinematic(strain(1,:),dl/2);
h_new(1,:) = h_old(1,:);
h_new(2,:) = (expG1*h_new(1,:)')';
h_new(3,:) = (expG1*h_new(2,:)')';
%when i>1
for i = 2:nelem
    expGi = elementKinematic(strain(i,:),dl/2);
    for j = 1:3
        node = 3*(i-1)+j;
        if j == 1
            h_new(node,:) = h_new(node-1,:);
        else
            h_new(node,:) = (expGi*h_new(node-1,:)')';  %row-->colum--->row
        end
    end
end
%pos_new
pos_new = zeros(nnode,3);
for i = 1:nnode
    pos_new(i,:) = h_new(i,1:3);
end
        
        
        
        
        
        
        
        
        
        