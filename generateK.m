%2017.11.8 by MY
function K = generateK(dl,ex,nelem,E,G,A,Ixx,Iyy,Izz)
%--------------------------------------------------------------------------
%To generate the element stiffness matrix and assemble them into total
%stiffness matrix.
%--------------------------------------------------------------------------
%Purpose:
%     To generate the element stiffness matrix and assemble them into total
%     stiffness matrix. 
%Synposis:
%     K = generateK(dl,nelem,E,G,Ixx,Iyy,Izz)
%Variables definition:
%Input:
%      dl  -  length of each element 
%      ex  -  extensional strain, dimenstion:nelem*1
%      nelem   -  number of elements
%      E   -  elastic modulus
%      G   -  shearing modulus
%      ----->G = 0.5*E/(1+mu)
%      A   -  cross-section aera
%      Ixx -  polar moment of inertia
%      Iyy -  y-bending moment of inertia
%      Izz -  z-bending moment of inertia
%      
%Output:
%      K   -  total stiffness matrix
%      ----->K = diag([Ke,Ke,...,Ke]), dimension:(4*nelem)*(4*nelem)
%-------------------------------------------------------------------------
Ke = zeros(4,4);
K = zeros(4*nelem,4*nelem);
% Ke
Ke(1,1) = E*A;
Ke(2,2) = G*Ixx;
Ke(3,3) = E*Iyy;
Ke(4,4) = E*Izz;
% deltaSe
deltaSe = dl*(1+ex);
% K
for i = 1:nelem
    K(4*i-3:4*i,4*i-3:4*i) = deltaSe(i)*Ke;
end
