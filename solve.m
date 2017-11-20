function [h,pos,strain,iterations] = solve(h0,pos0,strain0,fPoint,mPoint,hBC,posBC,dl,nelem,E,G,aera,Ixx,Iyy,Izz,maxiteration,eps,alphan)
%--------------------------------------------------------------------------
% Purpose:
%         To solve the static equations iteratively
% Synopsis :
%         [h,pos,strain] = solve(h0,pos0,strain0,hBC,maxiteration,eps)
% Variable Description:
% Input:
%           strain0 - strains of the undeformed beam
%           -----> strain0 = [e1;e2;...;en]   dimension: nelem*4
%           h0    - nodal vector of the undefored beam
%           -----> h0 = [h1;h2;...]   dimension:nnode*12
%           pos0   - position vector of the undefomed beam
%           -----> pos0 = [p1;p2;...]   dimension:nnode*3
%           fPoint  -  point force vector, dimension: (9*nelem)*1
%           mPoint  -  point moment vector, dimension: (9*nelem)*1
%           hBC  -  boundary condition defined in main program
%           -----> if hBC = 0,then, the beam is fixed in root
%           posBC - postion in which the boundary condition occured
%           -----> posBC = [ielem1,jnode1;ielem2,jnode2;...]
%           dl  -  length of each element 
%           nelem   -  number of elements
%           E   -  elastic modulus
%           G   -  shearing modulus
%           ----->G = 0.5*E/(1+mu)
%           aera -  cross-section aera
%           Ixx -  polar moment of inertia
%           Iyy -  y-bending moment of inertia
%           Izz -  z-bending moment of inertia
%           maxiteration - max number of iterations
%           eps  -   relative error to judge convergence
%           alphan  -  numerical damping coefficient
% Output:
%           h  - nodal vector of the deformed beam
%           -----> h = [h1;h2;...]   dimension:nnode*12
%           pos  - position vector of the deformed beam
%           -----> pos = [p1;p2;...]   dimension:nnode*3  
%           strain - strains of the deformed beam
%           -----> strain = [e1;e2;...;en]   dimension: nelem*4
%           iterations - iteration numbers up to convergence
%--------------------------------------------------------------------------
h = h0;
pos = pos0;
strain = strain0;
iterations = 0;
if hBC ==0
    unknown = matrix2col(strain);
else
    ielem = posBC(:,1);
    jnode = posBC(:,2);
    nodeid = 3*(ielem-1)+jnode;
    unknown = zeros(4*nelem+3*length(nodeid),1);
end

for iteration =1:maxiteration
    iterations = iterations+1;
    posi = pos;
    straini = strain;
    straini_col = matrix2col(straini);
    unknowni = unknown;
    exi = straini(:,1);
% calculate total stiffness matrix K at ith step
    KF = generateK(dl,exi,nelem,E,G,aera,Ixx,Iyy,Izz);
% calculate the generalized force Re at ith step
    [dhde,dpde,dthetade]=jacobian(strain,nelem,h,dl);
    Re = dpde'*fPoint+dthetade'*mPoint;
%boundary condition
    if hBC == 0
        K = KF;
        R = Re;
        row = 0;
    else
        Kce = dpde(3*nodeid-2:3*nodeid,:);
        Rce = Kce*straini_col-(posi(nodeid,:)'-hBC');
        row = size(Kce,1);
        R = [Re;Rce];
        K = [KF,Kce';Kce,zeros(row,row)];
    end   
% calculate the strains at (i+1)th step
    unknown = alphan*unknowni+(1-alphan)*K\R;
    strain_col = unknown(1:end-row);
% calculate new nodal vector (i+1)th
    strain = col2matrix(strain_col,4);
    [h,pos] = configRecovery(strain,h,dl);
    
% judging convergence
    residual = norm(unknown-unknowni)/norm(unknowni);
    if residual < eps 
        break;
    end
end
  