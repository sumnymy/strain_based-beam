function [expG1,dexpGdex1,dexpGdkx1,dexpGdky1,dexpGdkz1] = elementKinematic(straini,ds)
%--------------------------------------------------------------------------
% Purpose:
%         To calculate expontential expansion of (s-s0)Ke and its partial
%         derivation to strains in an element.
% Synopsis :
%         [exp,dexpde] = expG(strain,s,s0)
% Variable Description:
% Input:
%           straini - strains of ith element
%           -----> strain = [ex;kx;ky;kz]   dimension: 4*1
%           ds  -  original length of the integrated line  
% Output:
%          exp  -  expontential expansion of (s-s0)Ke  12*12
%          dexpde - partial derivation of exp to strain 12*12               
%--------------------------------------------------------------------------
% K
Ke = zeros(4,4);
Ke(1,2)= 1+straini(1);
kk = -skew(straini(2:4));
Ke(2:4,2:4) = kk;
%Éú³Étransformation matrix Th
Th = zeros(12,12);
Th(1,1)=1;Th(2,5)=1;Th(3,9)=1;
Th(4,2)=1;Th(5,6)=1;Th(6,10)=1;
Th(7,3)=1;Th(8,7)=1;Th(9,11)=1;
Th(10,4)=1;Th(11,8)=1;Th(12,12)=1;
% coefficients
deltas = ds*(1+straini(1));
lamda=norm(kk);
a0=1;
a1=deltas;
a2=(1-cos(deltas*lamda))/lamda^2;
a3=deltas/lamda^2-sin(deltas*lamda)/lamda^3;
% derivatives
% dkede
dKedex = zeros(4,4); dKedex(1,2) = 1;
dKedkx = zeros(4,4); dKedkx(3,4) = 1; dKedkx(4,3) = -1;
dKedky = zeros(4,4); dKedky(2,4) = -1;dKedky(4,2) = 1;
dKedkz = zeros(4,4); dKedkz(2,3) = 1; dKedkz(3,2) = -1;
% da0de=0
% da1de=0
% da2de
da2dex = 0;
da2Coeff = (lamda*deltas*sin(lamda*deltas)-2*(1-cos(lamda*deltas)))/lamda^4;
da2dkx = da2Coeff*straini(2);
da2dky = da2Coeff*straini(3);
da2dkz = da2Coeff*straini(4);
% da3de
da3dex = 0;
da3Coeff = 3*sin(lamda*deltas)/lamda^5-(2*deltas+deltas*cos(lamda*deltas))/lamda^4;
da3dkx = da3Coeff*straini(2);
da3dky = da3Coeff*straini(3);
da3dkz = da3Coeff*straini(4);
%divided into two situation to solve dexpKede
if lamda == 0
   % expKe
   expKe = a0*eye(4)+a1*Ke; 
   %dexpKedex
   dexpKedex = a1*dKedex;
   %dexpKedkx
   dexpKedkx = a1*dKedkx;
   %dexpKedky
   dexpKedky = a1*dKedky;
   %dexpKedkz
   dexpKedkz = a1*dKedkz;
else 
    % expKe
    expKe = a0*eye(4)+a1*Ke+a2*Ke*Ke+a3*Ke*Ke*Ke;
    %dexpKedex
    dexpKedex = da3dex*Ke*Ke*Ke+a3*(dKedex*Ke*Ke+Ke*dKedex*Ke+Ke*Ke*dKedex)+...
        da2dex*Ke*Ke+a2*(dKedex*Ke+Ke*dKedex)+a1*dKedex;
    %dexpGKedkx
    dexpKedkx = da3dkx*Ke*Ke*Ke+a3*(dKedkx*Ke*Ke+Ke*dKedkx*Ke+Ke*Ke*dKedkx)+...
        da2dkx*Ke*Ke+a2*(dKedkx*Ke+Ke*dKedkx)+a1*dKedkx;
    %dexpKedky
    dexpKedky = da3dky*Ke*Ke*Ke+a3*(dKedky*Ke*Ke+Ke*dKedky*Ke+Ke*Ke*dKedky)+...
        da2dky*Ke*Ke+a2*(dKedky*Ke+Ke*dKedky)+a1*dKedky;
    %dexpKedkz
    dexpKedkz = da3dkz*Ke*Ke*Ke+a3*(dKedkz*Ke*Ke+Ke*dKedkz*Ke+Ke*Ke*dKedkz)+...
        da2dkz*Ke*Ke+a2*(dKedkz*Ke+Ke*dKedkz)+a1*dKedkz;
end
%get expG & dexpGde
expG_bar = blkdiag(expKe,expKe,expKe);
expG = Th*expG_bar*Th';
dexpGdex_bar = blkdiag(dexpKedex,dexpKedex,dexpKedex);
dexpGdex = Th*dexpGdex_bar*Th';
dexpGdkx_bar = blkdiag(dexpKedkx,dexpKedkx,dexpKedkx);
dexpGdkx = Th*dexpGdkx_bar*Th';
dexpGdky_bar = blkdiag(dexpKedky,dexpKedky,dexpKedky);
dexpGdky = Th*dexpGdky_bar*Th';
dexpGdkz_bar = blkdiag(dexpKedkz,dexpKedkz,dexpKedkz);
dexpGdkz = Th*dexpGdkz_bar*Th';
% control the outputs
if nargout == 1
   expG1 = expG; 
else
    expG1 = expG;
    dexpGdex1 = dexpGdex;
    dexpGdkx1 = dexpGdkx;
    dexpGdky1 = dexpGdky;
    dexpGdkz1 = dexpGdkz;
end
