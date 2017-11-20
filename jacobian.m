%2017.11.9 by MY
%function [dhde,dpde,dthetade]=jacobian(strain,nelem,h,pos) 
function varargout=jacobian(strain,nelem,h,dl) 
%--------------------------------------------------------------------------
% Purpose:
%         To calculate Jacobian matrice
% Synopsis :
%         varargout=jacobian(strain,nelem,h,dl) 
%         [dhde,dpde,dthetade]=jacobian(strain,nelem,h,dl)
% Variable Description:
% Input:
%           strain - strains of the elements at this configuration
%           -----> strain = [ex1,kx1,ky1,kz1;ex2,kx2,ky2,kz2;...] 
%                  dimension: nelem*4
%           nelem - number of the elements       
%           nnode - number of the nodes
%           h     - nodal vector of the beam
%           -------->h = [h1;h2;h3;h4....]   dimension:nnode*12
%           dl  -  length of each element 
% Output:
%          dhde  -  jacobian matrix of dh/de  
%                  dimension:(36*nelem)*(4*nelem)
%          dpde  -  jacobian matrix of dp/de
%                  dimension:(9*nelem)*(4*nelem)
%          dthetade  -  jacobian matrix of dtheta/de
%                  dimension:(9*nelem)*(4*nelem)
%--------------------------------------------------------------------------
%boundary
hbc = zeros(12,nelem);
for i = 1:nelem
    hbc(:,i)=h(3*i-2,:)';
end
%dexpGde
expG = zeros(12,12,nelem);
dexpGdex = zeros(12,12,nelem);
dexpGdkx = zeros(12,12,nelem);
dexpGdky = zeros(12,12,nelem);
dexpGdkz = zeros(12,12,nelem);
dexp2Gdex = zeros(12,12,nelem);
dexp2Gdkx = zeros(12,12,nelem);
dexp2Gdky = zeros(12,12,nelem);
dexp2Gdkz = zeros(12,12,nelem);
for ii = 1:nelem

    [G,Gex,Gkx,Gky,Gkz] = elementKinematic(strain(ii,:),dl/2);
    expG(:,:,ii) = G;
    dexpGdex(:,:,ii)=Gex;
    dexpGdkx(:,:,ii)=Gkx;
    dexpGdky(:,:,ii)=Gky;
    dexpGdkz(:,:,ii)=Gkz;
    dexp2Gdex(:,:,ii) = G*Gex+Gex*G;
    dexp2Gdkx(:,:,ii) = G*Gkx+Gkx*G;
    dexp2Gdky(:,:,ii) = G*Gky+Gky*G;
    dexp2Gdkz(:,:,ii) = G*Gkz+Gkz*G;
end
%dhde
dhde = zeros(36*nelem,4*nelem);
for k = 1:nelem
  for i = 1:nelem
     for j = 1:3
         node = 3*(i-1)+j;
         if k>i
            dhde(12*node-11:12*node,4*k-3) = zeros(12,1);
            dhde(12*node-11:12*node,4*k-2) = zeros(12,1);
            dhde(12*node-11:12*node,4*k-1) = zeros(12,1);
            dhde(12*node-11:12*node,4*k) = zeros(12,1);
         elseif k == i
            if j == 1
            dhde(12*node-11:12*node,4*k-3) = zeros(12,1);
            dhde(12*node-11:12*node,4*k-2) = zeros(12,1);
            dhde(12*node-11:12*node,4*k-1) = zeros(12,1);
            dhde(12*node-11:12*node,4*k) = zeros(12,1);
            elseif j==2
            dhde(12*node-11:12*node,4*k-3) = dexpGdex(:,:,i)*hbc(:,i);
            dhde(12*node-11:12*node,4*k-2) = dexpGdkx(:,:,i)*hbc(:,i);
            dhde(12*node-11:12*node,4*k-1) = dexpGdky(:,:,i)*hbc(:,i);
            dhde(12*node-11:12*node,4*k) = dexpGdkz(:,:,i)*hbc(:,i);
            elseif j==3
            dhde(12*node-11:12*node,4*k-3) = dexp2Gdex(:,:,i)*hbc(:,i);
            dhde(12*node-11:12*node,4*k-2) = dexp2Gdkx(:,:,i)*hbc(:,i);
            dhde(12*node-11:12*node,4*k-1) = dexp2Gdky(:,:,i)*hbc(:,i);
            dhde(12*node-11:12*node,4*k) = dexp2Gdkz(:,:,i)*hbc(:,i); 
            end
         elseif k<i
            
            if j == 1
            dhde(12*node-11:12*node,4*k-3) = dhde(12*(node-1)-11:12*(node-1),4*k-3);
            dhde(12*node-11:12*node,4*k-2) = dhde(12*(node-1)-11:12*(node-1),4*k-2);
            dhde(12*node-11:12*node,4*k-1) = dhde(12*(node-1)-11:12*(node-1),4*k-1);
            dhde(12*node-11:12*node,4*k) = dhde(12*(node-1)-11:12*(node-1),4*k);
            elseif j==2
            dhde(12*node-11:12*node,4*k-3) = expG(:,:,i)*dhde(12*(node-1)-11:12*(node-1),4*k-3);
            dhde(12*node-11:12*node,4*k-2) = expG(:,:,i)*dhde(12*(node-1)-11:12*(node-1),4*k-2);
            dhde(12*node-11:12*node,4*k-1) = expG(:,:,i)*dhde(12*(node-1)-11:12*(node-1),4*k-1);
            dhde(12*node-11:12*node,4*k) = expG(:,:,i)*dhde(12*(node-1)-11:12*(node-1),4*k);
            elseif j==3
            dhde(12*node-11:12*node,4*k-3) = expG(:,:,i)*dhde(12*(node-1)-11:12*(node-1),4*k-3);
            dhde(12*node-11:12*node,4*k-2) = expG(:,:,i)*dhde(12*(node-1)-11:12*(node-1),4*k-2);
            dhde(12*node-11:12*node,4*k-1) = expG(:,:,i)*dhde(12*(node-1)-11:12*(node-1),4*k-1);
            dhde(12*node-11:12*node,4*k) = expG(:,:,i)*dhde(12*(node-1)-11:12*(node-1),4*k);
            end                
         end
     end
 end
end
%dpde
nnode = 3*nelem;
dpde = zeros(9*nelem,4*nelem);
for i = 1:nnode
    dpde(3*i-2:3*i,:)=dhde(12*i-11:12*i-9,:);
end
%dthetade
dthetade = zeros(0,0);   %initial value
for k = 1:nelem      %column
    dthetadei = zeros(0,0);
    for i = 1:nelem      %row
        dthetadeik = zeros(0,0);
        for j = 1:3
            node = 3*(i-1)+j;
            Gn = [h(node,4:6);h(node,7:9);h(node,10:12)];
            C_GE = Gn;   %rotation matrix
            node4 = 12*(node-1)+4;
            node5 = 12*(node-1)+5;
            node6 = 12*(node-1)+6;
            node7 = 12*(node-1)+7;
            node8 = 12*(node-1)+8;
            node9 = 12*(node-1)+9;
            node10 = 12*(node-1)+10;
            node11 = 12*(node-1)+11;
            node12 = 12*(node-1)+12; 
            exNode = 4*k-3;
            kxNode = 4*k-2;
            kyNode = 4*k-1;
            kzNode = 4*k;
            dGndex = [dhde(node4,exNode),dhde(node5,exNode),dhde(node6,exNode);...
                      dhde(node7,exNode),dhde(node8,exNode),dhde(node9,exNode);...
                     dhde(node10,exNode),dhde(node11,exNode),dhde(node12,exNode)];
            dGndkx = [dhde(node4,kxNode),dhde(node5,kxNode),dhde(node6,kxNode);...
                      dhde(node7,kxNode),dhde(node8,kxNode),dhde(node9,kxNode);...
                     dhde(node10,kxNode),dhde(node11,kxNode),dhde(node12,kxNode)];
            dGndky = [dhde(node4,kyNode),dhde(node5,kyNode),dhde(node6,kyNode);...
                      dhde(node7,kyNode),dhde(node8,kyNode),dhde(node9,kyNode);...
                      dhde(node10,kyNode),dhde(node11,kyNode),dhde(node12,kyNode)];
            dGndkz = [dhde(node4,kzNode),dhde(node5,kzNode),dhde(node6,kzNode);...
                      dhde(node7,kzNode),dhde(node8,kzNode),dhde(node9,kzNode);...
                      dhde(node10,kzNode),dhde(node11,kzNode),dhde(node12,kzNode)];     
            %dthetadex_E
            dthetadexSkew_G = -dGndex*C_GE';
            dthetadex_G = anti_skew(dthetadexSkew_G);
            dthetadex = C_GE'*dthetadex_G;
            %dthetadkx_E
            dthetadkxSkew_G = -dGndkx*C_GE';
            dthetadkx_G = anti_skew(dthetadkxSkew_G);
            dthetadkx = C_GE'*dthetadkx_G; 
            %dthetadky_E
            dthetadkySkew_G = -dGndky*C_GE';
            dthetadky_G = anti_skew(dthetadkySkew_G);
            dthetadky = C_GE'*dthetadky_G;
            %dthetadkz_E
            dthetadkzSkew_G = -dGndkz*C_GE';
            dthetadkz_G = anti_skew(dthetadkzSkew_G);
            dthetadkz = C_GE'*dthetadkz_G;
            %dthetade_E(j)
            dthetadeikj = [dthetadex,dthetadkx,dthetadky,dthetadkz]; %行—组装单元内部各应变分量
            dthetadeik = [dthetadeik;dthetadeikj];  %列—组装j
        end
        dthetadei = [dthetadei;dthetadeik];  %列—组装i
    end
    dthetade = [dthetade,dthetadei]; %行—组装k
end
varargout{1} = dhde;
varargout{2} = dpde;
varargout{3} = dthetade;
% if nargout == 1
%     dhde1 = dhde;
% elseif nargout ==2
%     dhde1 = dhde;
%     nnode = 3*nelem;
%     dpde = zeros(9*nelem,4*nelem);
%     for i = 1:nnode
%         dpde(3*i-2:3*i,:)=dhde(12*i-11:12*i-9,:);
%     end
%     dpde1 = dpde;
% end

