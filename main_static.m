%%Meng Yang, 2017.11
%This program is to solve the nonlinear beam problem using strain-based
%formulation. The varialbes are strains and lagrange multipliers
%representing the boundary conditions.
%

clear 
clc
disp('Please wait Programme is under Run')
%--------------------------------------------------------------------------
% % Geometrical and material properties of beam
%--------------------------------------------------------------------------
L = 1.04;                       %Length
E = 206*10^9;                %Elastic modulus
mu = 0.33;                   %poison ratio
G = 0.5*E/(1+mu);            %shearing modulus
width = 0.038;              %width of cross section
height = 0.001;            %height of cross-section
aera = width*height;         %aero of crosss-section
phu = 2700;                  %material density
mass = phu*aera;             %mass per unit length
Iyy = 1/12*width*height^3;   %y-bend moment of inertia 
Izz = 1/12*height*width^3;   %z-bend moment of inertia
Ixx = Iyy+Izz;               %polar moment of inertia

%%-----Element definition-----%%
nelem = 50;                   %total number of elements
strain = zeros(nelem,4);
%--------------------------------------------------------------------------
% Input data for nodal connectivity for each element
%--------------------------------------------------------------------------
[h,pos,nnode,dl]=mesh(strain,nelem,L);
%--------------------------------------------------------------------------
% Applied loads
%--------------------------------------------------------------------------
    %point force
tipforce = 2;
fPoint = zeros(9*nelem,1);
fPoint(end) = tipforce;
    %point moment
tipmoment = 0;
mPoint = zeros(9*nelem,1);
mPoint(end-1) = tipmoment;
%--------------------------------------------------------------------------
% Initial conditions,initial curvatures.
%--------------------------------------------------------------------------
h0 = h;
pos0 = pos;
strain0 = strain;
%--------------------------------------------------------------------------
% Boundary condition. Nodal displacement
%--------------------------------------------------------------------------
ielem = nelem;
jnode = 3;
posBC = [ielem,jnode];
nodeBC = 3*(ielem-1)+jnode;
pBC = pos0(nodeBC,:);   %row vector
%--------------------------------------------------------------------------
% Solve strains iteratively.
%--------------------------------------------------------------------------
maxiteration = 200;
eps = 10^-3;
damping = 0;

maxloadstep = 10;
P = 0;
deltaP = tipforce/maxloadstep;
fPointi = zeros(9*nelem,1);

for loadstep = 1:maxloadstep
    P = P+deltaP;
    fPointi(end) = P; 
[h,pos,strain,iteration] = solve(h,pos,strain,fPointi,mPoint,0,posBC,dl,nelem,E,G,aera,Ixx,Iyy,Izz,maxiteration,eps,damping);
end
% for iteration =1:maxiteration
%     posi = pos;
%     straini = strain;
%     straini_col = matrix2col(straini);
%     exi = straini(:,1);
% % calculate total stiffness matrix K at ith step
%     KF = generateK(dl,exi,nelem,E,G,aera,Ixx,Iyy,Izz);
% % calculate the generalized force Re at ith step
%     [dhde,dpde,dthetade]=jacobian(strain,nelem,h,dl);
%     Re = dpde'*fPoint+dthetade'*mPoint;
% % %boundary condition
% %         Kce = dpde(end-2:end,:);
% %         Rce = Kce*straini_col-(posi(end,:)'-pos0(end,:)');
% %         row = size(Kce,1);
% %         R = [Re;Rce];
% %         K = [KF,Kce';Kce,zeros(row,row)];
% % 
% % % calculate the strains at (i+1)th step
% %     unknown = K\R;
% %     strain_col = unknown(1:end-row);
% strain_col = KF\Re;
% % strain_col = damping*straini_col+(1-damping)*KF\Re;
% % calculate new nodal vector (i+1)th
%     strain = col2matrix(strain_col,4);
%     [h,pos] = configRecovery(strain,h,dl);
%     
% % judging convergence
%     residual = norm(strain_col-straini_col)/norm(straini_col);
%     if residual < eps 
%         break;
%     end
% end

%--------------------------------------------------------------------------
% Strain-configuration recovery.Post process
%--------------------------------------------------------------------------

configurationPlot(pos,pos0);























