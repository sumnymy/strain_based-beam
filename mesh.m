%2017.11.8 by MY
function [h,pos,nnode,dl]=mesh(strain0,nelem,L)
%To mesh a beam using 3 noded constant-strain elements
%----------------------------------------------------
%----------------------------------------------------
%Purpose:
%      To mesh a beam using 3 noded constant-strain elements
%Synopsis:
%      [h,pos,nnode,dl] = mesh(strain0,nelem,L)
%Variables definition:
%Input:
%      nelem  -  number of elements 
%      L      -  length of the beam
%      strain0 - initial strain
%Output:
%      h      - nodal vectors of the whole beam
%      -------->h = [h1;h2;h3;h4....]   dimension:nnode*12
%      pos    -  nodal position vector of the beam,
%      -------->pos = [node:x,y,z]  dimension:nnode*3
%      -------->pos_col = [x1;y1;z1;x2;y2;z2;...] dimention:(3*nnode)*1
%      baseVector - nodal base vectors of the beam
%      -------->baseVector =[node:wx,wy,wz]  dimension:nnode*9
%      -------->baseVector_col =[wx1;wy1;wz1;wx2;wy2;wy3;...]
%                                             dimension:(9*nnode)*1
%      nnode  -  number of nodes
%      dl     -  length of each element
%-------------------------------------------------------------------------
dl = L/nelem;
nnode = 3*nelem;
h =zeros(nnode,12);
pos = zeros(nnode,3);
baseVector = zeros(nnode,9);
for i = 1:nelem
    pos(3*i-2,:) = [(i-1)*dl,0,0];
    pos(3*i-1,:) = pos(3*i-2,:)+[dl/2,0,0];
    pos(3*i,:) = pos(3*i-1,:)+[dl/2,0,0];
end
for i = 1:nnode
    baseVector(i,:) = [1,0,0,0,1,0,0,0,1];
    h(i,:) = [pos(i,:),baseVector(i,:)];
end
[h,pos] = configRecovery(strain0,h,dl);

