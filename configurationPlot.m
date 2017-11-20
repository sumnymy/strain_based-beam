function configurationPlot(pos,pos0)
%--------------------------------------------------------------------------
% Purpose:
%         To plot the configuration of the beam reference line
% Synopsis :
%         configurationPlot(pos)
% Variable Description:
% Input:
%           pos    -  nodal position vector of the deformed beam,
%           pos0   -  nodal position vector of the undeformed beam,
% Output:
%--------------------------------------------------------------------------
nnode = size(pos,1);
nelem = nnode/3;
% coordinate 
repeatedNode = [];
for i = 1:nelem-1
    repeatedNode = [repeatedNode,3*i+1];
end
pos(repeatedNode,:) = [];
pos0(repeatedNode,:) = [];
coordinate = pos;x = coordinate(:,1);y = coordinate(:,2);z=coordinate(:,3);
coordinate0 = pos0; x0 = coordinate0(:,1);y0 = coordinate0(:,2);z0=coordinate0(:,3);
%plot configuration
figure;
% plot3(x,y,z,'r');
% hold on;
% plot3(x0,y0,z0,'b');
plot(x,z,'r',x0,z0,'b')
axis equal
